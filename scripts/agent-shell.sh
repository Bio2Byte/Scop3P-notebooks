#!/usr/bin/env bash
#
# Run a command against this repo inside a container, with the source mounted READ-ONLY.
#
# This is the default way for an automated agent to do shell work here: it keeps the
# host's zsh profile and brew tooling out of the picture, so results do not depend on
# what happens to be installed locally, and it makes it impossible for an exploratory
# command to modify the tree by accident.
#
#   scripts/agent-shell.sh                        # interactive bash, read-only source
#   scripts/agent-shell.sh grep -rn TODO apps     # one command, read-only source
#   scripts/agent-shell.sh --tests                # the full suite, read-only source
#   MODE=rw scripts/agent-shell.sh <cmd>          # writable, only when changing files
#   IMAGE=... scripts/agent-shell.sh <cmd>        # a different image
#
# Mounts:
#   /src      this repo. Read-only unless MODE=rw.
#   /scratch  writable scratch space, for anything that needs to write.
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Default image: native to the host, so investigation is fast. The project's own image is
# linux/amd64 (TM-align ships as an amd64 ELF) and therefore emulated on Apple Silicon;
# use it only when the real dependency set is actually needed -- see --tests below.
IMAGE="${IMAGE:-python:3.12-slim}"
TEST_IMAGE="${TEST_IMAGE:-bio2byte/scop3p-base:local}"
MODE="${MODE:-ro}"
SCRATCH="${SCRATCH:-${TMPDIR:-/tmp}/scop3p-agent-scratch}"

if [[ "${MODE}" == "rw" ]]; then
  MOUNT_FLAG=""
  echo "agent-shell: source mounted READ-WRITE (MODE=rw)" >&2
else
  MOUNT_FLAG=":ro"
fi

mkdir -p "${SCRATCH}"

run_in_container() {
  local image="$1"; shift
  # Pin the platform to the image's own, so an amd64 project image on Apple Silicon runs
  # under emulation quietly instead of warning on every invocation.
  local platform
  platform="$(docker image inspect "${image}" --format '{{.Os}}/{{.Architecture}}' 2>/dev/null || true)"
  docker run --rm -i ${TTY_FLAG:-} ${platform:+--platform "${platform}"} \
    -v "${REPO_ROOT}:/src${MOUNT_FLAG}" \
    -v "${SCRATCH}:/scratch" \
    -w /src \
    -e PYTHONPATH=/src/apps \
    -e PYTHONDONTWRITEBYTECODE=1 \
    -e SCOP3P_LOG_DIR=/scratch/logs \
    -e HOME=/scratch \
    "${image}" \
    "$@"
}

# Interactive shell when given nothing to run. bash, not zsh: the slim images do not ship
# zsh, and bash avoids the zsh-only quoting traps (glob expansion of an unquoted
# --include=*.py, ":local" being read as a history modifier).
if [[ $# -eq 0 ]]; then
  TTY_FLAG="-t" run_in_container "${IMAGE}" bash
  exit $?
fi

if [[ "$1" == "--tests" ]]; then
  shift
  # Default to the whole suite, but if the caller names any path, run only those. Passing
  # both a directory and a file made pytest collect a narrowed set while still *looking*
  # like a full run -- a passing count that means less than it appears is worse than no
  # count at all. The resolved command is echoed so the scope is never in doubt.
  # Split the arguments into paths to run and flags to pass through. An argument counts
  # as a path only if it actually exists, so a flag's value (the "foo" in -k foo) is not
  # mistaken for one. Targets are then passed exactly once: passing them again inside the
  # flag list made pytest run every selected test twice, which doubled the reported pass
  # count -- a number that looks like more coverage while being the same tests.
  targets=()
  flags=()
  for arg in "$@"; do
    if [[ "${arg}" != -* && -e "${REPO_ROOT}/${arg}" ]]; then
      targets+=("${arg}")
    else
      flags+=("${arg}")
    fi
  done
  if [[ ${#targets[@]} -eq 0 ]]; then
    targets=(tests/unit tests/integration)
  fi
  echo "agent-shell: pytest ${targets[*]} ${flags[*]:-}" >&2
  # -p no:cacheprovider: pytest would otherwise try to write .pytest_cache into the
  # read-only source mount and fail before collecting anything.
  run_in_container "${TEST_IMAGE}" \
    python -m pytest "${targets[@]}" -q -p no:cacheprovider ${flags[@]+"${flags[@]}"}
  exit $?
fi

run_in_container "${IMAGE}" "$@"
