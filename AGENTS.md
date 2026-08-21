# Working agreements for automated agents

Conventions for coding agents working in this repository. Written for agents, but the
reasoning applies to anyone.

## Run shell work in a container, with the source read-only

**Do not run investigative shell commands on the host.** Run them inside a container with
this repo mounted, and mount it **read-only** by default.

```bash
scripts/agent-shell.sh grep -rn "TODO" apps      # one command, source read-only
scripts/agent-shell.sh                           # interactive bash, source read-only
scripts/agent-shell.sh --tests                   # full suite, source read-only
scripts/agent-shell.sh --tests -k sifts          # arguments pass through to pytest
MODE=rw scripts/agent-shell.sh <cmd>             # writable -- only to change files
IMAGE=<image> scripts/agent-shell.sh <cmd>       # a different image
```

Two reasons, both practical:

- **The host is someone's working machine.** Running there picks up their `zsh` profile and
  their `brew` tooling, so a result may depend on whatever they happen to have installed,
  and any dev tool installed to test an idea pollutes their environment. A container is
  clean, fresh and reproducible.
- **Read-only makes accidents impossible.** Searching, reading, and testing a hypothesis
  never need write access, so they should not have it. Mount read-write only when the task
  is actually to change files.

### Which image

| Purpose | Image | Notes |
|---|---|---|
| Searching, reading, ad-hoc Python | `python:3.12-slim` | native `arm64` on Apple Silicon, fast |
| Real dependency set: b2bTools, TM-align, the test suite | `bio2byte/scop3p-base:local` | `linux/amd64`, so emulated on Apple Silicon |

The project images are `linux/amd64` because `TM-align` is committed as an amd64 ELF
(`ELF 64-bit LSB pie executable, x86-64`), which is also why every compose service inherits
`platform: linux/amd64` from the shared `x-app` anchor in `docker-compose.yml`. Emulation costs roughly 2x:
the full suite takes about 15s in the container against 7.5s natively. That is a fine price
for a clean environment, but it is why the default image for investigation is the native one.

### Use bash, not zsh

The containers run **bash**. The slim images do not ship zsh, and bash avoids the zsh-only
quoting traps that have actually caused wrong results here:

- `grep -rn "x" --include=*.py .` — zsh expands the unquoted `*.py` and the grep fails with
  "no matches found", which reads exactly like "no occurrences in the codebase".
- `"$var:local"` — zsh treats `:local` as a history modifier.

### What a read-only mount breaks, and the fix

Anything that writes into the tree fails. The wrapper already handles the common cases:

- `PYTHONDONTWRITEBYTECODE=1`, so no `__pycache__` writes are attempted.
- `-p no:cacheprovider` for pytest, which would otherwise try to create `.pytest_cache`
  in the source tree and fail before collecting a single test.
- `SCOP3P_LOG_DIR=/scratch/logs` and `HOME=/scratch`, so logs and any tool state land on
  the writable `/scratch` mount.

Write anything else that needs a filesystem to `/scratch` too.

### Two things that will catch you out

**Importing anything from `common` needs the project image.** `apps/common/__init__.py`
imports `services`, which imports pandas, so `python:3.12-slim` cannot import even a
dependency-free module like `common.structure_labels`. Use the slim image for text work
(grep, AST, reading) and `IMAGE=bio2byte/scop3p-base:local` the moment you need to import
project code.

**The hardened base image has no pip.** It was stripped deliberately. To use a linter or
another dev tool, install it in the slim image instead:

```bash
scripts/agent-shell.sh bash -c 'pip install --quiet pyflakes && python -m pyflakes apps/'
```

That works because pyflakes only parses source -- it does not import it, so it does not
need the project dependencies.

## Run the apps in containers too

The same applies to running a protocol for verification: mount the working tree read-only
into the project image and invoke shiny directly, so what runs is the code you just edited
rather than whatever is baked into the image.

```bash
docker run -d --name scop3p-agent-<app> --platform linux/amd64 \
  -p 8057:8057 \
  -v "$PWD:/src:ro" -v /tmp/agent-scratch:/scratch -w /src \
  -e PYTHONPATH=/src/apps -e SCOP3P_LOG_DIR=/scratch -e HOME=/scratch \
  bio2byte/scop3p-base:local \
  python -m shiny run --host 0.0.0.0 --port 8057 apps/<app>/app.py
```

Name your containers `scop3p-agent-*` so it is unambiguous which ones are yours to remove,
and give logs a writable `/scratch` mount -- the default log directory is not writable and
the app would otherwise fall back to a temp path you then have to hunt for.

## Do not disturb a running container

Someone may be manually testing a protocol while an agent works. Never run `docker stop`,
`docker rm`, `docker compose down`, or `docker rmi` against a container you did not start,
and keep off port `8000` — that is where the toolkit is usually published. Start your own
dev server on a free port instead, and stop only what you started.

## Logging

There are no `print` statements in `apps/`, and tests enforce that. Everything goes through
`common.logging_utils`, which also provides the per-session experiment trail and the
`quiet_third_party` wrapper for libraries that print. See the "Logs" section of
[`apps/README.md`](apps/README.md) before adding any output.
