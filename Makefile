# Scop3P-Toolkit image builds.
#
# The app images are one line each on top of a shared base, so the base has to exist
# first. That ordering is the whole reason this file exists -- `docker compose build`
# cannot express it.
#
#   make base                    build the shared runtime
#   make apps                    build every app image
#   make scop3p-toolkit          build just the published portal image
#   make scan                    vulnerability scan of the base
#   make scan-app APP=rinalign   vulnerability scan of one app image
#   make sizes                   report image sizes
#   make clean                   remove the locally built images
#
# Overridable: NAMESPACE, VERSION, BASE_TAG, PLATFORM.

NAMESPACE ?= bio2byte
VERSION   ?= local
BASE_TAG  ?= $(NAMESPACE)/scop3p-base:$(VERSION)
PLATFORM  ?= linux/amd64

# TM-align is committed as a prebuilt amd64 ELF, so the images are amd64-only.
export DOCKER_DEFAULT_PLATFORM = $(PLATFORM)

BUILD_DATE := $(shell date -u +%Y-%m-%dT%H:%M:%SZ)
VCS_REF    := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)

APPS := peptide-mapper structure-viz mutation-effect topology-viewer rinalign scop3p-toolkit

BUILD_ARGS := \
	--build-arg BUILD_DATE=$(BUILD_DATE) \
	--build-arg VCS_REF=$(VCS_REF) \
	--build-arg VERSION=$(VERSION)

.PHONY: all base apps $(APPS) scan scan-app sizes clean help

all: apps

help:
	@sed -n 's/^# \{0,2\}//p' $(MAKEFILE_LIST) | sed -n '1,20p'

# The base carries the venv, the binaries and apps/, and has no CMD, so it is not
# runnable and cannot be confused for an app image.
base:
	docker build -f docker/Dockerfile.base $(BUILD_ARGS) -t $(BASE_TAG) .

# Each app image is FROM the base, adding only SCOP3P_APP_NAME and a CMD. The build
# context is docker/apps, so these take about a second once the base exists.
$(APPS): base
	docker build -f docker/apps/$@.Dockerfile \
		--build-arg BASE_IMAGE=$(BASE_TAG) \
		-t $(NAMESPACE)/$@:$(VERSION) docker/apps

apps: $(APPS)

# Vulnerability scan. Reports CRITICAL and HIGH; the base is the interesting target
# because every app image inherits from it unchanged.
scan: base
	docker run --rm -v /var/run/docker.sock:/var/run/docker.sock \
		aquasec/trivy:latest image --scanners vuln \
		--severity CRITICAL,HIGH $(BASE_TAG)

scan-app:
	@test -n "$(APP)" || { echo "usage: make scan-app APP=<name>"; exit 2; }
	docker run --rm -v /var/run/docker.sock:/var/run/docker.sock \
		aquasec/trivy:latest image --scanners vuln \
		--severity CRITICAL,HIGH $(NAMESPACE)/$(APP):$(VERSION)

sizes:
	@docker images --format '{{.Repository}}:{{.Tag}}\t{{.Size}}' \
		| grep -E '^$(NAMESPACE)/(scop3p-base|$(shell echo $(APPS) | tr ' ' '|'))' || true

clean:
	-docker rmi -f $(BASE_TAG) $(foreach app,$(APPS),$(NAMESPACE)/$(app):$(VERSION))
