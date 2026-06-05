# Agent Harness — the user interface to the pipeline.
#
# Every target is a one-liner over Docker + rclone. The contract: you should
# never need to remember a `docker run` invocation or a conda command.
#
#   make build                         # build the image locally
#   make run TARGET=1M17 MODE=production
#   make run TARGET=2HYY MODE=targetdiff NUM=30 GPU=1
#   make cloud-run TARGET=2HYY MODE=targetdiff NUM=30   # rent a GPU, run, sync
#   make pull / make push              # sync data/ with Cloudflare R2
#   make dashboard                     # regenerate the static dashboard
#
# No Docker? Use the *-local targets (run via conda instead). Needs the `base`
# (+ `targetdiff_env` for cascade/targetdiff) conda envs — see README
# "Running locally without Docker":
#   make run-local TARGET=1IEP MODE=cascade NUM=5     # full pipeline, no Docker
#   make dashboard-local                              # regenerate dashboard, no Docker
#
# Override any capitalised variable on the command line: `make run TARGET=6P3D`.

# --- Configuration -----------------------------------------------------------
IMAGE      ?= ghcr.io/vladsft/autonomous-drug-discovery:latest
REPO_ROOT  := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
DATA_DIR   := $(REPO_ROOT)/autonomous_drug_discovery/data
DASH_DIR   := $(REPO_ROOT)/dashboard

# Run parameters (override on the command line).
TARGET     ?= 1M17
MODE       ?= simulation
NUM        ?=
GPU        ?=
CAMPAIGN   ?=
# No-Docker (`*-local`) targets run via conda instead of Docker. DEVICE selects
# the compute device for GPU-capable generation; on a CPU-only box leave it
# `auto` (it detects no GPU and uses cpu). CONDA overrides the conda binary.
DEVICE     ?= auto
CONDA      ?= conda

# Cloudflare R2: an rclone remote named `r2` (see .env.example / make bootstrap).
R2_BUCKET  ?= agent-harness
R2_REMOTE  := r2:$(R2_BUCKET)

# Derived docker arguments.
DATA_MOUNT := -v $(DATA_DIR):/app/autonomous_drug_discovery/data
RUN_AS     := --user $(shell id -u):$(shell id -g) -e HOME=/tmp
GPU_FLAG   := $(if $(GPU),--gpus all,)
NUM_FLAG   := $(if $(NUM),--num_samples $(NUM),)
DOCKER_RUN := docker run --rm $(GPU_FLAG) $(RUN_AS) $(DATA_MOUNT)
# No-Docker invocation: run the orchestrator/dashboard directly in the `base`
# conda env (the env must exist — see "Running locally without Docker" in README).
CONDA_BASE := $(CONDA) run --no-capture-output -n base python

.DEFAULT_GOAL := help
.PHONY: help bootstrap build pull push run run-local cloud-run test dashboard dashboard-local deploy logs clean

help: ## Show this help
	@grep -hE '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk -F':.*?## ' '{printf "  \033[1m%-14s\033[0m %s\n", $$1, $$2}'

bootstrap: ## One-time setup: pull the image, check rclone + credentials
	docker pull $(IMAGE)
	@command -v rclone >/dev/null || { echo "rclone not installed — see https://rclone.org/install/"; exit 1; }
	@test -f $(REPO_ROOT)/.env || { echo "Missing .env — copy .env.example to .env and fill it in."; exit 1; }
	@echo "Bootstrap OK. Try: make run TARGET=1M17 MODE=simulation"

build: ## Build the Docker image locally (CI builds the canonical one)
	docker build -t $(IMAGE) $(REPO_ROOT)

# pull/push source .env so host rclone picks up the RCLONE_CONFIG_R2_* vars —
# one credentials file drives both `make` and the cloud pod.
pull: ## Download campaign data + telemetry from Cloudflare R2
	set -a; . $(REPO_ROOT)/.env; set +a; rclone copy --progress $(R2_REMOTE) $(DATA_DIR)

push: ## Upload local campaign data + telemetry to Cloudflare R2
	set -a; . $(REPO_ROOT)/.env; set +a; rclone copy --progress $(DATA_DIR) $(R2_REMOTE)

run: ## Run the full pipeline locally (TARGET=, MODE=, NUM=, GPU=1)
	$(DOCKER_RUN) $(IMAGE) \
		orchestrator.py run data/processed/$(TARGET).pdb --mode $(MODE) $(NUM_FLAG)

run-local: ## No-Docker run via conda (TARGET=, MODE=cascade, NUM=, DEVICE=auto)
	@test -f $(DATA_DIR)/processed/$(TARGET).pdb || { echo "Missing $(DATA_DIR)/processed/$(TARGET).pdb"; exit 1; }
	$(CONDA_BASE) $(REPO_ROOT)/autonomous_drug_discovery/orchestrator.py run \
		$(DATA_DIR)/processed/$(TARGET).pdb --mode $(MODE) --device $(DEVICE) $(NUM_FLAG)

cloud-run: ## Provision a RunPod GPU, run the pipeline, sync to R2, tear down
	TARGET=$(TARGET) MODE=$(MODE) NUM=$(NUM) IMAGE=$(IMAGE) \
		bash $(REPO_ROOT)/scripts/cloud_run.sh

test: ## Run the test suite inside the image
	docker run --rm $(RUN_AS) $(IMAGE) \
		bash -lc 'cd /app && conda run --no-capture-output -n base python -m pytest -q autonomous_drug_discovery/tests'

dashboard: ## Regenerate the static dashboard from local telemetry
	$(DOCKER_RUN) -v $(DASH_DIR):/app/dashboard $(IMAGE) \
		/app/scripts/regenerate_dashboard.py --data-dir data --out /app/dashboard

dashboard-local: ## No-Docker dashboard regen via conda (reads telemetry → dashboard/)
	$(CONDA_BASE) $(REPO_ROOT)/scripts/regenerate_dashboard.py \
		--data-dir $(DATA_DIR) --out $(DASH_DIR)

deploy: dashboard ## Regenerate the dashboard, then hand off to CI for Pages
	@echo "Dashboard regenerated. Commit dashboard/ and push — CI deploys to GitHub Pages."

logs: ## Tail a campaign's logs from R2 (make logs CAMPAIGN=campaign_xxxx)
	@test -n "$(CAMPAIGN)" || { echo "Set CAMPAIGN=campaign_xxxx"; exit 1; }
	set -a; . $(REPO_ROOT)/.env; set +a; rclone cat $(R2_REMOTE)/$(CAMPAIGN)/run.log

clean: ## Remove local campaign outputs and Python cruft (mirrors .gitignore)
	rm -rf $(DATA_DIR)/campaign_* \
	       $(DATA_DIR)/candidates $(DATA_DIR)/screened $(DATA_DIR)/results \
	       $(DATA_DIR)/logs $(DATA_DIR)/outputs
	find $(DATA_DIR)/processed -maxdepth 1 -name '*_p2rank' -type d -exec rm -rf {} + 2>/dev/null || true
	find $(DATA_DIR)/processed -maxdepth 1 -name '*_manifest.json' -delete 2>/dev/null || true
	find $(REPO_ROOT) -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
	find $(REPO_ROOT) -name '.pytest_cache' -type d -exec rm -rf {} + 2>/dev/null || true
	find $(REPO_ROOT) -name '.ruff_cache' -type d -exec rm -rf {} + 2>/dev/null || true
