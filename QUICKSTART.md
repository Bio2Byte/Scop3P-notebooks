# Quickstart

Every protocol in this repository can be opened in four ways. Pick the row that matches
what you want to do:

| I want to… | Use | Needs |
|---|---|---|
| Just use the tools, all of them, in one place | **Docker launcher** | Docker |
| Run one tool in isolation | **Docker, one app per container** | Docker |
| Develop or modify an app | **Shiny app locally** | Python 3.12 |
| Read and modify the science | **Jupyter notebook** | Python 3.12, or Binder |
| Use a notebook as an app, without seeing code | **Voilà** | Python 3.12, or Binder |

---

## Fastest path: the launcher in Docker

The launcher (`scop3p-toolkit`) is a single container holding every app, with a navbar
to switch between them. Build it:

```bash
make scop3p-toolkit
```

Run it:

```bash
docker run --rm -p 8000:8000 bio2byte/scop3p-toolkit:local
```

Then open **<http://localhost:8000>**.

`make` is the entry point because the images come in two layers: a shared base
(`docker/Dockerfile.base`) carrying the Python environment, the binaries and `apps/`,
and a one-line image per app (`docker/apps/<app>.Dockerfile`) adding just its name and
its command. The base has to exist first, and `make` handles that ordering. Which app
you get is decided by which Dockerfile you build, so there is no `--target` and no way
for stage ordering to hand you the wrong app.

Other targets: `make base`, `make apps` (all six), `make sizes`, `make scan`.

`TM-align` is committed as a prebuilt amd64 ELF binary, so the images are amd64-only;
the Makefile sets `DOCKER_DEFAULT_PLATFORM=linux/amd64` for you. On Apple Silicon this
runs under emulation, so the first base build takes a few minutes and later ones are
cached. App images build in about a second each once the base exists.

### Or pull the published image

The published image is `bio2byte/scop3p-toolkit`; `v0.3.2` was the first release with
all five apps:

```bash
docker run --rm -p 8000:8000 bio2byte/scop3p-toolkit:v0.3.3
```

- **Pin an explicit version, `v0.3.2` or newer.** Older tags predate the Topology
  Viewer, RIN Alignment and the peptide upload mode — `v0.2.5` gives you three apps,
  not five.
- **`latest` tracks recent releases.** It is pushed by hand for the `v0.3.x` releases
  and moves automatically with CI from `v0.4.0` onward. Pinning an explicit version is
  still the reproducible choice — `latest` tells you nothing about which release you
  got.

### Jump straight to one app

The launcher accepts a preselect parameter, which is also what the navbar links use:

| App | URL |
|---|---|
| Structure Visualisation | <http://localhost:8000/?app=structure-viz> |
| RIN Alignment | <http://localhost:8000/?app=rinalign> |
| Mutation Effect | <http://localhost:8000/?app=mutation-effect> |
| Peptide Mapper | <http://localhost:8000/?app=peptide-mapper> |
| Topology Viewer | <http://localhost:8000/?app=topology-viewer> |
| Help | <http://localhost:8000/?app=help> |

Your choice is remembered in a cookie, so the next plain visit to `/` reopens it.

**Not sure which tool you need?** The navbar's **Help?** entry lists the mission, scope
and use cases of every protocol, with a question-first table at the top — pick the
question you arrived with and it points you at the tool.

Every accession field has a **Load example** button, so you can see any protocol work
without knowing an accession to type.

---

## What each protocol is called, in each format

| Protocol | Notebook | Shiny app | Docker target | Port |
|---|---|---|---|---|
| Peptide Mapper | `notebooks/Peptide_mapper_scop3p_voila.ipynb`<br>`notebooks/Peptide_mapper_fileupload_voila.ipynb` | `apps/peptide_mapper/app.py` | `peptide-mapper` | 8001 |
| Structure Visualisation | `notebooks/Scop3P_PTM_structure_viz_voila_app.ipynb` | `apps/structure_viz/app.py` | `structure-viz` | 8002 |
| Mutation Effect | `notebooks/Scop3P_b2b_mutation_effect_voila_app.ipynb` | `apps/mutation_effect/app.py` | `mutation-effect` | 8003 |
| Topology Viewer | `notebooks/topology_viewer/topology_viewer.ipynb` | `apps/topology_viewer/app.py` | `topology-viewer` | 8004 |
| RIN Alignment | `notebooks/RINAlign_align_and compare_networks.ipynb` | `apps/rinalign/app.py` | `rinalign` | 8005 |
| py3Dmol peptide mapper | `notebooks/scop3p_py3dmol_mapper.ipynb` | — | — | — |
| *(the launcher)* | — | `apps/portal/main.py` | `scop3p-toolkit` | 8000 |

The two Peptide Mapper notebooks are one Shiny app with two source tabs: *Scop3P
peptides* and *Upload your own*. For migration status and known gaps per protocol, see
[`docs/PROTOCOL_STATUS.md`](docs/PROTOCOL_STATUS.md).

---

## 1. Docker, one app per container

Use this to run a single tool on its own port, without the launcher. Build the shared
base once, then bring up whichever app you want:

```bash
make base
```

```bash
docker compose up -d topology-viewer
```

Open <http://localhost:8004>. Substitute any target from the table above — the host port
differs per app, the container always listens on 8000.

Stop it:

```bash
docker compose stop topology-viewer
```

Everything at once, launcher plus all five single-app containers on ports 8000–8005:

```bash
make apps && docker compose up -d
```

Note that CI only publishes `scop3p-toolkit`. The per-app images are a local
convenience; there is nothing to pull for them.

### How large is it, and how safe?

```bash
make sizes
```

```bash
make scan
```

`make scan` runs Trivy against the base image, which is where essentially all of the
surface lives — the app images add only an `ENV` and a `CMD`. See the security notes in
[`apps/README.md`](apps/README.md) for the current posture and the two remaining
fixable findings.

---

## 2. Shiny apps locally, without Docker

Install the app stack once:

```bash
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-shiny.txt -r requirements-biophysics.txt
```

`requirements-biophysics.txt` pins CPU-only PyTorch. Measured on linux/amd64, that file
alone installs into a 1.3 GB virtualenv instead of 4.9 GB, with none of the twelve CUDA
packages the default `torch` wheel drags in. Only Structure Visualisation and Mutation
Effect actually need it — they run b2bTools — but installing it keeps all five apps
runnable from one environment.

Run one app:

```bash
PYTHONPATH=apps .venv/bin/shiny run --reload --port 8050 apps/topology_viewer/app.py
```

Open <http://localhost:8050>. Change the file and the port for other apps.

**`PYTHONPATH=apps` is not optional for `peptide_mapper`, `structure_viz` and
`mutation_effect`.** `shiny run` puts the app file's own directory on `sys.path`, not
`apps/`, so `from common import …` fails without it. `topology_viewer` and `rinalign` add
`apps/` themselves and work either way; the prefix is harmless, so just always use it.

The launcher is an ASGI app rather than a Shiny app, so it runs under uvicorn:

```bash
PYTHONPATH=apps .venv/bin/python -m uvicorn portal.main:app --port 8000
```

---

## 3. Jupyter notebooks

### On Binder, nothing to install

Open any notebook straight from the [Jupyter Notebook index in the
README](README.md#jupyter-notebook-index), which carries a Binder badge per protocol.

### Locally

The notebooks use a **different dependency set** from the Shiny apps — Voilà, ipywidgets,
nglview and py3Dmol, none of which the apps need:

```bash
python3.12 -m venv .venv-nb && .venv-nb/bin/pip install -r requirements.txt
```

```bash
.venv-nb/bin/jupyter lab
```

Then open a notebook from `notebooks/`. Two quirks worth knowing:

- **The Topology Viewer notebook must be started from its own directory.** It is a
  three-line launcher over the `topology/` package sitting next to it, so `topology` has
  to be importable from the working directory:

  ```bash
  cd notebooks/topology_viewer && ../../.venv-nb/bin/jupyter lab topology_viewer.ipynb
  ```

- **The RINAlign filename contains a space.** Quote it in any shell command:
  `"notebooks/RINAlign_align_and compare_networks.ipynb"`.

---

## 4. Voilà

Voilà renders a notebook as an app: the widgets, none of the code cells. It is the
original delivery format for these protocols, and the Shiny apps are ports of it.

### On Binder

Same index as above — each protocol has an *Interactive app (Voilà)* badge beside its
notebook badge. If a Voilà link errors, open the *Notebook (JupyterLab)* link first to
authenticate, then retry.

### Locally

Using the notebook environment from the previous section:

```bash
.venv-nb/bin/voila --port 8080 notebooks/Peptide_mapper_scop3p_voila.ipynb
```

It opens a browser at <http://localhost:8080> and executes the notebook on load, so the
first paint is slow — tens of seconds for the protocols that import b2bTools and torch,
or that call the Scop3P and AlphaFold APIs before rendering. A blank page for a while is
normal. The Topology Viewer again needs its own working directory:

```bash
cd notebooks/topology_viewer && ../../.venv-nb/bin/voila --port 8081 topology_viewer.ipynb
```

To browse and launch any notebook instead of pinning one, point Voilà at the folder —
**with an absolute path**:

```bash
.venv-nb/bin/voila --port 8080 "$PWD/notebooks"
```

A relative directory does not work with the pinned Voilà 0.5.8: it passes the path
straight through as the content manager's root, and current `jupyter_core` rejects a
non-absolute root, so the tree page loads but every notebook you click returns HTTP 500
(`abs_root=PosixPath('notebooks') is not absolute`). Relative paths to a single
*notebook file* are fine — only directory mode is affected.

---

## Which format should I use?

- **Docker launcher** — the intended way to *use* the toolkit, and what gets deployed as
  a Galaxy interactive tool. Multi-user safe.
- **Docker per app** — a single tool on a fixed port, for embedding or a focused demo.
- **Shiny locally** — developing an app. `--reload` restarts on save. Note that editing
  the `topology/` package under `notebooks/` needs a manual restart, because its browser
  assets are read into module constants at import.
- **Notebook** — reading, adapting or extending the science, and the place where new
  protocols start.
- **Voilà** — a widget UI for a protocol that has no Shiny app yet, or for comparing a
  port against the original.

---

## If something does not work

**The navbar only shows three apps.** You are running an old published image. Pull
`v0.3.2` or newer (see the top of this document); tags before it predate the Topology
Viewer, RIN Alignment and the peptide upload mode.

**`ModuleNotFoundError: No module named 'common'`.** Missing `PYTHONPATH=apps` on a local
`shiny run`.

**`ModuleNotFoundError: No module named 'topology'`.** The Topology Viewer notebook was
started from the wrong directory; `cd notebooks/topology_viewer` first. For the *Shiny*
app this cannot happen — it resolves the package through
`apps/common/topology_bridge.py`, and if resolution fails the app shows a card naming
every path it tried. `SCOP3P_TOPOLOGY_PATH` overrides the location.

**The Docker build fails on `notebooks/topology_viewer/topology`.** Something is
excluding `notebooks/` from the build context. `.dockerignore` must not.

**An app build fails with `pull access denied` or `manifest unknown` for
`bio2byte/scop3p-base:local`.** The base image has not been built yet. Run `make base`
(the `make <app>` targets do it for you; a bare `docker compose build` does not).

**The 3D panel stays blank.** All 3D viewers (NGL, Mol\*, py3Dmol, D3) load from public
CDNs at runtime. They need outbound network access from the *browser*, not the server. In
the Topology Viewer, the status line inside the view reports which engine loaded.

**An app tab shows a spinner forever.** Report it — that is the signature of a Shiny
output suspended in a hidden tab, and it needs `@output(suspend_when_hidden=False)` on
that output. See the note at the end of [`apps/README.md`](apps/README.md).

**A protocol reports zero modifications, or fails with `JSONDecodeError: Expecting
value: line 1 column 1`.** The Scop3P endpoint it is calling has moved. Scop3P serves
its single-page app from a catch-all route, so a retired endpoint answers `200 OK` with
HTML rather than a 404 — which is why the failure surfaces as a JSON decode error deep
inside `requests`, or as a silent zero if the caller swallows it. The Shiny apps target
the v1 API and report this properly, as do the notebooks. If you see it, check
<https://iomics.ugent.be/scop3p/api/v1/openapi.json> for the current paths.

**Where are the logs?** Inside a container, `/var/log/scop3p_toolkit/`, alongside a
`metadata.yml` recording the session's dependency and tool versions. `docker-compose.yml`
bind-mounts them to `./logs/<service>/`. Override with `SCOP3P_LOG_DIR`.

---

## More detail

- [`apps/README.md`](apps/README.md) — per-app features, dependencies, and how to add an app
- [`docs/PROTOCOL_STATUS.md`](docs/PROTOCOL_STATUS.md) — what is migrated, what is not, and why
- [`docs/use-cases/`](docs/use-cases/) — per-app parity specs against the source notebooks
- [`README.md`](README.md) — project background, Scop3P citation, Binder index

## What a run leaves behind

Every protocol writes a step-by-step record of what you did, so a figure can be traced
back to the inputs and actions that produced it.

```bash
ls logs/structure-viz/
# scop3p_toolkit_log_<stamp>.log     everything: the record plus diagnostics
# scop3p_toolkit_trail_<stamp>.log   the record on its own
# metadata.yml                       versions, tools, session start, log paths
```

The trail file reads in order:

```
step=3  action=input  detail="UniProtKB accession = P07949"
step=4  action=click  detail="Set protein"
step=5  action=result detail="34 PDB entries cross-referenced from P07949"
```

Each line carries a session id, so two browser tabs stay separable. Blocked actions and
failures appear in the same sequence at `WARNING` and `ERROR`, so a run that went wrong
shows *where* it went wrong.

For more detail on a misbehaving prediction, raise the level:

```bash
SCOP3P_LOG_LEVEL=DEBUG
```

That surfaces the third-party output (b2bTools progress, dependency warnings) which is
captured to `DEBUG` rather than printed, so it never floods the console at normal levels.
See "Logs" in [`apps/README.md`](apps/README.md) for the full scheme.

