"""
Agent Planner — The Decision Engine.

LangChain-based autonomous agent that plans and executes a drug discovery campaign.
This is the "Intelligence Layer" described in the North Star vision document.

The Agent:
  1. Evaluates the target by querying telemetry for past runs on similar targets.
  2. Dispatches modules (ingestion → generation → screening → docking).
  3. Reviews attrition rates and adapts strategy.
  4. Produces a final recommendation report.

Usage:
  python agent_planner.py --target target.pdb --mode dry_run --max_iterations 3

Environment variables:
  DISCOVERY_LLM_API_KEY — API key for the LLM backbone (Gemini, OpenAI, etc.)
  DISCOVERY_LLM_PROVIDER — "google" (default), "openai", or "anthropic"

Dependencies:
  pip install langchain langchain-core langchain-google-genai
"""

import os
import sys
import json
import argparse
import subprocess
import uuid
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.resolve()
sys.path.insert(0, str(PROJECT_ROOT))

from telemetry import TelemetryDB

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DEFAULT_MAX_ITERATIONS = 5
DATA_DIR = PROJECT_ROOT / "data"

# Map generation mode → conda env. Keeps the subprocess launcher consistent
# with the orchestrator's routing so the planner works for every mode.
_ENV_FOR_MODE = {
    "targetdiff": "targetdiff_env",
    "pocket2mol": "pocket2mol_env",
}


def _conda_python(env_name: str) -> list[str]:
    """Return the command prefix for running python in a specific conda env."""
    conda_bin = os.environ.get("CONDA_EXE", "conda")
    return [conda_bin, "run", "-n", env_name, "python"]

# Check for LangChain availability
try:
    from langchain_core.tools import tool
    from langchain_core.messages import HumanMessage, SystemMessage
    HAS_LANGCHAIN = True
except ImportError:
    HAS_LANGCHAIN = False

try:
    from langchain_google_genai import ChatGoogleGenerativeAI
    HAS_GOOGLE_AI = True
except ImportError:
    HAS_GOOGLE_AI = False


# ---------------------------------------------------------------------------
# Tool Definitions
# ---------------------------------------------------------------------------

def _run_module(cmd, desc):
    """Run a subprocess command and return result."""
    print(f"  [Tool] {desc}...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        output = result.stdout + result.stderr
        if result.returncode != 0:
            return f"FAILED (exit code {result.returncode}): {output[-500:]}"
        return f"SUCCESS: {output[-500:]}"
    except subprocess.TimeoutExpired:
        return "FAILED: Timeout after 300s"
    except Exception as e:
        return f"FAILED: {e}"


def make_tools(campaign_id: str, db_path: str, data_dir: Path, exec_mode: str = "simulation"):
    """Create the tool functions for the agent, bound to a specific campaign.

    Returns plain callables. If LangChain is available, the functions are
    decorated as @tool so they can be bound to an LLM.
    """

    gen_env = _ENV_FOR_MODE.get(exec_mode, "base")

    # Only the `simulation` mode uses stub docking; every real generation
    # backend (rdkit / pocket2mol / targetdiff / production) gets real Vina.
    dock_mode = "simulation" if exec_mode == "simulation" else "production"

    def run_ingestion(pdb_path: str) -> str:
        """Run P2Rank (default) or fpocket to identify binding pockets.
        Returns the path to the manifest JSON file."""
        output_dir = str(data_dir / "01_ingestion")
        cmd = _conda_python("base") + [
            str(PROJECT_ROOT / "modules" / "01_ingestion" / "run_pocket.py"),
            "--pdb", pdb_path,
            "--output_dir", output_dir,
            "--db_path", db_path,
            "--campaign_id", campaign_id,
        ]
        return _run_module(cmd, f"Ingesting {pdb_path}")

    def run_generation(manifest_path: str) -> str:
        """Generate candidate ligands using the configured backend
        (simulation / rdkit / pocket2mol / targetdiff). Takes a manifest JSON
        from ingestion. Returns path to generated SDF file."""
        output_dir = str(data_dir / "02_generation")
        cmd = _conda_python(gen_env) + [
            str(PROJECT_ROOT / "modules" / "02_generation" / "run_generation.py"),
            "--manifest", manifest_path,
            "--output_dir", output_dir,
            "--mode", exec_mode,
            "--db_path", db_path,
            "--campaign_id", campaign_id,
        ]
        return _run_module(cmd, f"Generating molecules (mode={exec_mode}, env={gen_env})")

    def run_screening(sdf_path: str) -> str:
        """Screen molecules through drug-likeness filters + PAINS + ADMET-AI.
        Returns path to screening report JSON."""
        output_dir = str(data_dir / "03_screening")
        cmd = _conda_python("base") + [
            str(PROJECT_ROOT / "modules" / "03_screening" / "run_screening.py"),
            "--input_sdf", sdf_path,
            "--output_dir", output_dir,
            "--db_path", db_path,
            "--campaign_id", campaign_id,
        ]
        return _run_module(cmd, "Screening molecules")

    def run_docking(manifest_path: str, candidates_dir: str) -> str:
        """Score molecules via AutoDock Vina (production) or return
        simulation scores. Returns path to docking results CSV."""
        output_dir = str(data_dir / "04_docking")
        cmd = _conda_python("base") + [
            str(PROJECT_ROOT / "modules" / "04_docking" / "run_docking.py"),
            "--manifest", manifest_path,
            "--candidates_dir", candidates_dir,
            "--output_dir", output_dir,
            "--mode", dock_mode,
            "--db_path", db_path,
            "--campaign_id", campaign_id,
        ]
        return _run_module(cmd, f"Docking (mode={dock_mode})")

    def query_telemetry(query_type: str) -> str:
        """Query the telemetry database. query_type can be:
        'campaign_summary' — summary of current campaign,
        'all_runs' — list all runs,
        'molecules' — list scored molecules from the latest screening run."""
        db = TelemetryDB(db_path)
        try:
            if query_type == "campaign_summary":
                summary = db.get_campaign_summary(campaign_id)
                return json.dumps(summary, indent=2, default=str)
            elif query_type == "all_runs":
                runs = db.query_runs(campaign_id=campaign_id)
                return json.dumps(runs[:20], indent=2, default=str)
            elif query_type == "molecules":
                runs = db.query_runs(campaign_id=campaign_id, module_name="03_screening")
                if runs:
                    mols = db.query_molecules(runs[0]["run_id"])
                    return json.dumps(mols[:50], indent=2, default=str)
                return "No screening runs found."
            else:
                return f"Unknown query type: {query_type}. Use: campaign_summary, all_runs, molecules"
        finally:
            db.close()

    def read_file(file_path: str) -> str:
        """Read the contents of a JSON or CSV file for analysis. Returns the first 2000 characters."""
        try:
            with open(file_path) as f:
                content = f.read(2000)
            return content
        except Exception as e:
            return f"Error reading file: {e}"

    all_fns = [run_ingestion, run_generation, run_screening, run_docking, query_telemetry, read_file]

    # Wrap with LangChain @tool if available (needed for LLM tool-calling mode)
    if HAS_LANGCHAIN:
        from langchain_core.tools import tool as lc_tool
        wrapped = []
        for fn in all_fns:
            wrapped.append(lc_tool(fn))
        return wrapped

    # Without LangChain: return plain functions with .name attribute for tool_map lookup
    for fn in all_fns:
        fn.name = fn.__name__  # type: ignore[attr-defined]
        fn.invoke = lambda args, _fn=fn: _fn(**args) if isinstance(args, dict) else _fn(args)  # type: ignore[attr-defined]
    return all_fns


# ---------------------------------------------------------------------------
# Agent Controller
# ---------------------------------------------------------------------------

SYSTEM_PROMPT = """You are the Adaptive Discovery Orchestrator — an autonomous decision engine for computational drug discovery.

Your mission: Given a target protein PDB file, plan and execute a discovery campaign to identify promising small-molecule drug candidates.

You have these tools available:
- run_ingestion: Identify binding pockets on the target protein (P2Rank by default)
- run_generation: Generate candidate ligands. Backend is fixed per campaign by --mode
  (simulation, rdkit, pocket2mol, targetdiff, or production)
- run_screening: Filter molecules for drug-likeness, SA, PAINS, and ADMET
- run_docking: Score surviving molecules via AutoDock Vina (or stubs in simulation mode)
- query_telemetry: Check campaign progress and historical data
- read_file: Read output files for analysis

Standard workflow:
1. Run ingestion on the PDB to get a manifest
2. Run generation using the manifest
3. Screen the generated molecules
4. Dock the survivors
5. Query telemetry for the campaign summary
6. Write a final recommendation

Important rules:
- Always check the output of each step before proceeding
- If a step fails, note the error and try to recover or adapt
- Pay attention to attrition rates. If >95% of molecules are filtered, note this as a concern
- Be decisive. You have a limited number of iterations.
"""


class AgentController:
    """LLM-driven agent loop for drug discovery campaigns."""

    def __init__(self, target_pdb: str, mode: str = "simulation",
                 max_iterations: int = DEFAULT_MAX_ITERATIONS,
                 db_path: str | None = None):

        self.target_pdb = Path(target_pdb).resolve()
        self.mode = mode
        self.max_iterations = max_iterations

        # Campaign setup
        self.campaign_id = f"campaign_{uuid.uuid4().hex[:8]}"
        self.db_path = db_path or str(PROJECT_ROOT / "data" / "telemetry.db")
        self.data_dir = DATA_DIR / self.campaign_id
        self.data_dir.mkdir(parents=True, exist_ok=True)

        # Create tools (works with or without LangChain)
        self.tools = make_tools(self.campaign_id, self.db_path, self.data_dir, mode)

        # Decision log
        self.decisions = []

    def run(self) -> str:
        """Execute the discovery campaign."""
        print(f"\n{'='*60}")
        print(f"  Adaptive Discovery Orchestrator")
        print(f"  Campaign:  {self.campaign_id}")
        print(f"  Target:    {self.target_pdb.name}")
        print(f"  Mode:      {self.mode}")
        print(f"  Max Iters: {self.max_iterations}")
        print(f"{'='*60}\n")

        # Log agent start to telemetry
        db = TelemetryDB(self.db_path)
        agent_run_id = db.start_run(
            campaign_id=self.campaign_id,
            module_name="agent_planner",
            input_path=str(self.target_pdb),
            parameters={
                "mode": self.mode,
                "max_iterations": self.max_iterations,
                "llm_available": HAS_LANGCHAIN and HAS_GOOGLE_AI,
            },
        )

        try:
            if HAS_LANGCHAIN and HAS_GOOGLE_AI and os.environ.get("DISCOVERY_LLM_API_KEY"):
                report = self._run_with_llm()
            else:
                print("[Agent] No LLM API key found. Running deterministic pipeline.")
                report = self._run_deterministic()

            # Log completion
            db.complete_run(agent_run_id, "success", str(self.data_dir / "campaign_report.md"))

            # Write report
            report_path = self.data_dir / "campaign_report.md"
            with open(report_path, "w") as f:
                f.write(report)

            print(f"\n[Agent] Campaign report: {report_path}")
            return report

        except Exception as e:
            db.complete_run(agent_run_id, "failed", error_trace=str(e))
            raise
        finally:
            db.close()

    def _run_with_llm(self) -> str:
        """Run with LLM-driven decision making."""
        llm = ChatGoogleGenerativeAI(
            model="gemini-2.0-flash",
            google_api_key=os.environ["DISCOVERY_LLM_API_KEY"],
        )

        # Bind tools to LLM
        llm_with_tools = llm.bind_tools(self.tools)

        messages = [
            SystemMessage(content=SYSTEM_PROMPT),
            HumanMessage(content=(
                f"Execute a drug discovery campaign on target: {self.target_pdb}\n"
                f"Campaign ID: {self.campaign_id}\n"
                f"Mode: {self.mode}\n"
                f"Data directory: {self.data_dir}\n"
                f"Begin by running ingestion on the PDB file."
            ))
        ]

        for iteration in range(self.max_iterations):
            print(f"\n--- Agent Iteration {iteration + 1}/{self.max_iterations} ---")

            response = llm_with_tools.invoke(messages)
            messages.append(response)

            # Process tool calls
            if response.tool_calls:
                for tc in response.tool_calls:
                    tool_name = tc["name"]
                    tool_args = tc["args"]
                    print(f"  [Agent] Calling: {tool_name}({tool_args})")

                    self.decisions.append({
                        "iteration": iteration + 1,
                        "tool": tool_name,
                        "args": tool_args,
                        "timestamp": datetime.now(timezone.utc).isoformat(),
                    })

                    # Execute tool
                    tool_fn = next((t for t in self.tools if t.name == tool_name), None)
                    if tool_fn:
                        result = tool_fn.invoke(tool_args)
                        from langchain_core.messages import ToolMessage
                        messages.append(ToolMessage(content=str(result), tool_call_id=tc["id"]))
                    else:
                        from langchain_core.messages import ToolMessage
                        messages.append(ToolMessage(content=f"Unknown tool: {tool_name}", tool_call_id=tc["id"]))
            else:
                # No more tool calls — agent is done
                return response.content or self._generate_summary()

        return self._generate_summary()

    def _run_deterministic(self) -> str:
        """Fallback: run a fixed pipeline sequence without LLM reasoning."""
        steps = [
            ("run_ingestion", {"pdb_path": str(self.target_pdb)}),
            ("run_generation", None),  # args determined from previous output
            ("run_screening", None),
            ("run_docking", None),
            ("query_telemetry", {"query_type": "campaign_summary"}),
        ]

        tool_map = {t.name: t for t in self.tools}
        outputs = {}

        # Step 1: Ingestion
        print("\n[Agent] Step 1/5: Ingestion")
        result = tool_map["run_ingestion"].invoke(steps[0][1])
        print(f"  Result: {result[:200]}")
        outputs["ingestion"] = result

        # Find manifest
        manifest_path = self.data_dir / "01_ingestion"
        manifests = list(manifest_path.glob("*_manifest.json"))
        if not manifests:
            return self._failure_report("Ingestion failed — no manifest found.", outputs)
        manifest = str(manifests[0])

        # Step 2: Generation
        print("\n[Agent] Step 2/5: Generation")
        result = tool_map["run_generation"].invoke({"manifest_path": manifest})
        print(f"  Result: {result[:200]}")
        outputs["generation"] = result

        gen_dir = self.data_dir / "02_generation"
        sdf_files = list(gen_dir.glob("*.sdf"))
        if not sdf_files:
            return self._failure_report("Generation failed — no SDF found.", outputs)
        sdf_path = str(sdf_files[0])

        # Step 3: Screening
        print("\n[Agent] Step 3/5: Screening")
        result = tool_map["run_screening"].invoke({"sdf_path": sdf_path})
        print(f"  Result: {result[:200]}")
        outputs["screening"] = result

        # Step 4: Docking
        print("\n[Agent] Step 4/5: Docking")
        screening_dir = str(self.data_dir / "03_screening")
        result = tool_map["run_docking"].invoke({
            "manifest_path": manifest,
            "candidates_dir": screening_dir,
        })
        print(f"  Result: {result[:200]}")
        outputs["docking"] = result

        # Step 5: Summary
        print("\n[Agent] Step 5/5: Campaign Summary")
        result = tool_map["query_telemetry"].invoke({"query_type": "campaign_summary"})
        outputs["summary"] = result

        return self._generate_summary(outputs)

    def _generate_summary(self, outputs=None):
        """Build the campaign report from outputs."""
        report = f"""# Campaign Report: {self.campaign_id}

**Date**: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}
**Target**: {self.target_pdb.name}
**Mode**: {self.mode}
**Iterations**: {len(self.decisions)} tool calls

## Pipeline Execution Summary

"""
        if outputs:
            for step, result in outputs.items():
                status = "✅" if "SUCCESS" in str(result) else "❌"
                report += f"### {step.replace('_', ' ').title()} {status}\n"
                report += f"```\n{str(result)[:300]}\n```\n\n"

        if self.decisions:
            report += "## Decision Log\n\n"
            report += "| # | Tool | Args |\n|---|---|---|\n"
            for d in self.decisions:
                report += f"| {d['iteration']} | {d['tool']} | {json.dumps(d['args'])[:60]} |\n"

        report += f"\n## Data Directory\n`{self.data_dir}`\n"
        return report

    def _failure_report(self, reason, outputs):
        """Generate a failure report."""
        report = f"# Campaign FAILED: {self.campaign_id}\n\n**Reason**: {reason}\n\n"
        for step, result in outputs.items():
            report += f"### {step}\n```\n{str(result)[:200]}\n```\n\n"
        return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Adaptive Discovery Orchestrator — Agent Planner")
    parser.add_argument("--target", required=True, help="Target PDB file")
    parser.add_argument("--mode",
                        choices=["simulation", "rdkit", "pocket2mol", "targetdiff", "production"],
                        default="simulation",
                        help="Generation backend (production = rdkit + Vina docking)")
    parser.add_argument("--max_iterations", type=int, default=DEFAULT_MAX_ITERATIONS,
                        help="Maximum agent iterations (dead-man switch)")
    parser.add_argument("--db_path", default=None, help="Telemetry DB path")
    args = parser.parse_args()

    controller = AgentController(
        target_pdb=args.target,
        mode=args.mode,
        max_iterations=args.max_iterations,
        db_path=args.db_path,
    )
    report = controller.run()
    print("\n" + report)


if __name__ == "__main__":
    main()
