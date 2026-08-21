# Scop3P Shiny Apps

This directory contains the five converted Python Shiny apps, a help page, and one
wrapper portal:

- [`apps/peptide_mapper/app.py`](peptide_mapper/app.py): Peptide Mapper
- [`apps/structure_viz/app.py`](structure_viz/app.py): Structure Visualisation
- [`apps/topology_viewer/app.py`](topology_viewer/app.py): Topology Viewer
- [`apps/mutation_effect/app.py`](mutation_effect/app.py): Mutation Effect
- [`apps/rinalign/app.py`](rinalign/app.py): RIN Alignment
- [`apps/help/app.py`](help/app.py): Help page — mission, scope and use cases per protocol
- [`apps/portal/main.py`](portal/main.py): single-root all-in-one selector

For which notebooks have and have not been converted yet, see [`docs/PROTOCOL_STATUS.md`](../docs/PROTOCOL_STATUS.md).

For how to launch any of them, see [`QUICKSTART.md`](../QUICKSTART.md).

## Apps

### Peptide Mapper
- Loads peptide and modification data from Scop3P, or from an uploaded TSV/CSV
- Auto-detects the five relevant columns in an uploaded search-engine export
- Filters and selects peptide spans
- Maps selections on AlphaFold structures
- Exports a styled HTML session, the PDB, and a mapped-residue TSV

### Structure Visualisation
- Fetches PTMs from Scop3P and, optionally, annotated PTMs from UniProt, merged so a
  site described by both is listed once
- Fetches disease-associated variants
- Renders AlphaFold structures, or a PDB entry picked from a dropdown of the structures
  cross-referenced from the accession's own UniProt record -- each option labelled with its
  method, resolution and chains, cascading into a chain picker that shows the UniProt range
  that chain covers
- Runs Bio2Byte predictions and 3D coloring
- Builds residue interaction networks, coloured either by PTM/variant status or by a
  predicted Bio2Byte property (fill carries the value, border keeps the site status)
- Translates UniProt positions onto each PDB entry's own residue numbering via SIFTS, and
  names the method on screen so the reader knows whether marks were placed authoritatively
  or inferred
- Runs TM-align comparisons

### Mutation Effect
- Fetches WT UniProt sequence and Scop3P PTMs
- Runs Bio2Byte WT and mutant predictions
- Compares WT vs mutant residue features
- Generates mutation-centric inference summaries

### Topology Viewer
- Draws a 2D secondary-structure topology diagram beside a 3D viewer
- Fetches an AlphaFold DB model or a PDBe entry by accession, or takes an upload
- Reads secondary structure from the file, else derives it from coordinates
  (`mkdssp`, then `biotite`, then a built-in P-SEA with no dependencies)
- Overlays Scop3P PTMs and UniProt disease variants, in accession mode only
- Three layouts, two 3D engines (NGL and Mol*)

### RIN Alignment
- Builds residue interaction networks from two structures
- Diffs two models of the same protein into conserved, lost and gained contacts
- Aligns two different proteins by Weisfeiler-Lehman topology plus Hungarian assignment
- Overlays Scop3P PTMs and UniProt disease variants on the network views
- Five views, including a linked network-and-3D panel

### Help
- One card per protocol: mission, scope, use cases, and what input it needs
- A question-first table at the top: pick the question, get the tool
- Static content only, so it has no reactive outputs and cannot hit the
  hidden-output suspension trap noted at the end of this file
- Not a protocol, so it has no Docker image of its own; it is reachable through the
  portal navbar and at `/?app=help`

### PTM sources and residue codes

The PTM table merges two sources. Scop3P contributes experimentally observed
modifications; UniProt contributes annotated PTMs of every kind (acetylation,
methylation, glycosylation, lipidation), which is what the "Include UniProt PTMs"
checkbox controls. Site identity is accession + three-letter residue + position: where
both sources describe the same site it is listed once, keeping the Scop3P naming and
folding in the UniProt reference.

`residue` is always a **three-letter code**, produced by
`common.structure_viz.residue_three_letter`. This matters beyond tidiness: Scop3P v1
reports the descriptive name ("Phosphoserine"), and both the 3D viewer's colour map and
the UniProt site key expect `SER`. While `residue` held the descriptive name, every PTM
rendered in the viewer's fallback colour and no Scop3P row could ever match a UniProt
one. The descriptive text is preserved in a separate `modification` column.

The viewer's colour map covers SER/THR/TYR for phospho-acceptors plus ASN/LYS/CYS for
the glycosylation, acetylation and lipidation sites the UniProt source adds. A unit test
pins that contract.

## Shared UI Vocabulary

Everything user-facing that repeats across apps lives in
[`apps/common/ui_shell.py`](common/ui_shell.py), so it can only be inconsistent in one
place:

| Helper | Purpose |
|---|---|
| `scop3p_shell(name, intro, *children)` | Page shell: hero, CSS, title |
| `scop3p_card(title, *children)` | One panel |
| `scop3p_footer()` | Shared footer and logo strip |
| `ACCESSION_LABEL` | **The** name for a UniProt identifier field |
| `scop3p_field_row(*children)` | An input and its buttons on one baseline |
| `scop3p_example_button(id)` | The "Load example" affordance |

Three conventions worth knowing:

- **One name for the accession field.** `ACCESSION_LABEL` is `"UniProtKB accession"`.
  The apps previously each invented their own — `ACC_ID (UniProt accession number)`,
  `UniProt`, `UniProt accession (AlphaFold DB / PDBe)` — which made one toolkit read as
  five unrelated tools. Import the constant; do not retype the label. Status and guard
  strings use the same wording.
- **Inputs and their buttons share a baseline.** Shiny renders a text input as a label
  above a control with a margin below, while an action button has no label, so a button
  dropped into an adjacent column floats above the field it belongs to.
  `scop3p_field_row` pins its children to the bottom edge and drops the trailing
  margin. It wraps on narrow columns, which is why the buttons stack below the input in
  the narrower sidebars. Note that the row pins button `width`/`flex`, because an app
  stylesheet setting `.btn { width: 100% }` for its own button grid would otherwise
  stretch these across the row — scope such rules to the grid that needs them.
- **Every accession field has a "Load example" button.** Each app defines an
  `EXAMPLE_ACCESSION` constant, shows it in the placeholder, and a `load_example`
  handler fills the field and says what to click next. The placeholder documents the
  example; the button makes it usable in one click.
- **Result cards match the height of their controls.** The two-column grids
  (`.scop3p-two-col`, `.pm-main-grid`, `.tv-main-grid`) use `align-items: stretch`, so
  an empty results panel is never a short box floating beside a tall controls card. Do
  not set `align-items: start` on them.

### Choosing a structure

Every protocol that lets the user pick a structure offers a **dropdown**, never a free-text
entry code, and every one of them describes an entry the same way. The wording lives in
[`apps/common/structure_labels.py`](common/structure_labels.py) so it cannot drift apart
again -- it had, into three notations for the same facts:

| Protocol | Picker | Label |
|---|---|---|
| Structure Visualisation | entry, then chain | `2IVS · X-ray · 2.00 A · chains A, B` |
| Topology Viewer | entry, then chain | `6Q2O · Electron Microscopy · 3.65 A · chains E, F · 56% cover` |
| RIN Alignment | entry-and-chain in one | `9C5S · chain A · UniProt 17-30 · X-ray · 1.01 A` |

RIN Alignment selects a chain directly because a residue interaction network is built on
exactly one chain; the other two pick an entry and then cascade into its chains, labelled
with the UniProt range that chain covers (`A (705-1013)`). The AlphaFold option is called
`AlphaFold model (full length)` everywhere.

Options are only ever what the accession actually cross-references, so a user cannot select
a structure that does not contain the protein. Placeholders distinguish the two cases that
look alike in an empty dropdown: *Set a protein first* against *No PDB entries for this
protein*.

`StructureRef.label()` in `notebooks/topology_viewer/topology/` is deliberately left alone:
that dataclass is shared with the notebook and its test suite, so the Shiny app formats its
own labels instead of changing it.

### Waiting time

These protocols do slow things -- a Bio2Byte prediction, a KD-tree over a few hundred
residues, TM-align in a subprocess, a network round trip per annotation source. Every action
that waits uses a **task button** from [`apps/common/busy.py`](common/busy.py), which:

- **disables itself the instant it is clicked**, client-side, so a second click cannot queue
  the same work again (measured: disabled 3ms after the click, and zero of a burst of clicks
  accepted across an 11-second wait);
- shows a **spinner and a "Working…" label** for the duration;
- **re-enables when the work ends -- on success and on failure alike.**

Two facts about Shiny drive the design, both measured against a running app rather than
assumed:

**A synchronous handler cannot show a live status message.** Anything an effect writes is
queued and flushed only after the effect returns, so a handler that sets "Working…" and then
works for three seconds displays nothing until it is finished. The button's own busy state
still works, because that is client-side.

**Moving work off the loop needs a thread, not just `async`.** `reactive.extended_task` runs
its coroutine on the event loop, so a blocking call inside `async def` still pins it;
`asyncio.to_thread` is what actually frees it. That is also what stops one slow upstream
freezing *every* connected session.

So there are two tiers, and which one an action uses is a judgement about its wait:

| Helper | Behaviour | Used for |
|---|---|---|
| `task_button` | disable + spinner + auto re-enable; status message appears when the work ends | most slow actions |
| `background_task_button` + `background(...)` | all of the above, plus a **live** status message during the work, and the loop stays free for other sessions | the Bio2Byte predictions, the longest wait in the suite |

`background_task_button` sets `auto_reset=False` (with the default the button re-enables on
the next flush -- measured at 177ms into a 3-second task, so it would go live while the work
was still running). That brings one obligation: **every exit path must call `finish_task`**,
early returns included. Miss one and the button stays dead for the rest of the session with
nothing on screen to explain it. A test walks each app's AST to enforce the pairing, because
that failure is invisible to every other kind of test.

Instant, local actions -- *Load example*, *Clear*, *Reset* -- stay plain buttons on purpose;
a spinner that flashes for 20ms is noise, not feedback. A test pins that too.

### Lookups are cached process-wide

Six protocols run in one process behind the portal and they ask the same upstreams for the
same things: the UniProt FASTA for an accession is fetched by Mutation Effect *and* by
Structure Visualisation, Scop3P modifications by four of them, and every browser session
repeated the lot. Against upstreams that fail intermittently, each avoidable request was
another chance to fail in front of a user.

[`apps/common/cache.py`](common/cache.py) memoizes these lookups **once per process**, so
they are shared across sessions and across protocols. Measured on `P07949`: the first fetch
takes 0.26s, and the same sequence requested by a different protocol and then by a
different session both return in under a millisecond.

| Cache | Filled by |
|---|---|
| `uniprot.sequence.fasta` | Structure Visualisation **and** Mutation Effect, sharing one entry |
| `uniprot.pdb.xrefs` | the structure pickers |
| `uniprot.ptm.features`, `uniprot.disease.variants` | annotation sources |
| `scop3p.modifications` | every protocol that reads Scop3P |
| `uniprot.entry.info` | RIN Alignment |
| `b2b.prediction.table` | **Run predictions** on the Bio2Byte tab |
| `b2b.prediction.raw` | Mutation Effect's wild-type prediction |

**Bio2Byte predictions** are the most expensive thing the toolkit computes -- about 16
seconds for 1100 residues -- and deterministic for a given sequence, so they are the
biggest single win: pressing **Run predictions** again, or opening the protocol in another
session, returns in under a millisecond. They also carry the sharpest correctness risk in
the whole cache, and two separate protections are in place:

- **The key includes the sequence, not just the accession.** Mutation Effect predicts a
  *mutant* under the wild type's accession, so an accession-only key would hand back the
  wild type's numbers for the mutant -- silently destroying the comparison that protocol
  exists to make, with no error anywhere.
- **Mutants are not cached at all.** A mutant sequence is seen once, so a slot spent on it
  never earns a hit and, worse, would evict the wild-type prediction every later comparison
  needs: exploring a dozen mutations would keep re-running the one prediction worth keeping.
  The caller declares which it has (`wild_type=False`) rather than this being guessed.

Both are wanted, not one or the other: the bypass is a call-site decision and call sites get
edited, so the sequence in the key is what still protects a caller who forgets it.

The two prediction caches are deliberately separate -- one returns a normalised DataFrame,
the other b2bTools' raw payload -- and are bounded to 32 entries rather than 256, because a
prediction frame is orders of magnitude larger than a JSON lookup.

Four properties are deliberate, and each is pinned by test because losing any of them
would be invisible:

- **Failures are never cached.** The one that matters most. These upstreams fail
  transiently, so a cached exception would make the retry button useless and the only cure
  a restart -- turning a blip into an outage.
- **The key ignores the receiver.** Each session builds its own service instance, so keying
  on it would give every session a private cache: no sharing at all, while still looking
  like a working cache. Whether to drop the first argument is decided from the *signature*
  (a first parameter named `self`), not from the call -- dropping it unconditionally would
  discard the real argument of a plain function and collapse every accession onto one entry.
- **Results are copied on the way out, on both the store and hit paths.** A DataFrame handed
  to two sessions is one object, and either could corrupt the other's view. Copying only on
  the way *in* protects the first caller and leaves every later one holding the live object.
- **The name is the cache identity.** Reusing a name is how two protocols share a result;
  building a fresh cache per function would look shared while sharing nothing, and would
  orphan all but the last from the registry so `clear_all()` silently missed them.

Entries expire after an hour and each cache is bounded to 256 entries (LRU). `cache_report()`
returns hit/miss counters per cache.

**Structure files** are handled separately, in a process-wide directory
(`SCOP3P_STRUCTURE_CACHE_DIR`, default `$TMPDIR/scop3p_structures`) rather than in memory: a
downloaded PDB entry is an immutable upstream artefact, so 2IVT is the same 226 KB for every
session. Downloads take a per-filename lock, so two tabs asking at the same moment do not
both fetch it or read a half-written file. Generated files -- trimmed chain segments,
rendered HTML, TM-align output -- stay in the per-session workdir, because they are that
session's output rather than a copy of something upstream.

**Why not memcached.** It was considered and it is the wrong shape here: the data is a few MB
per accession, does not need to outlive the process, and there is one process per container.
A cache server would add a service to every image, a network hop, serialisation, and a new
failure mode in exchange for nothing this needs. The point to revisit is if the toolkit is
ever scaled to several worker processes or replicas -- then the sharing becomes genuinely
cross-process and the trade changes.

### A failed lookup does not look like an empty result

An empty dropdown after a failed request looked exactly like a protein with no structures,
so the user was told something false about their protein and had no reason to press the
button again. The two now read differently:

| Situation | What the picker says |
|---|---|
| No accession set yet | *Set a protein first* |
| Lookup succeeded, nothing found | *No PDB entries for this protein* |
| Lookup **failed** | *Lookup failed - press Set protein to retry* (or *Fetch*, per app) |

RIN Alignment sets the same message on both structure pickers when its fetch raises. The
status card already explained what happened, but the dropdown is where the user is looking.

### Annotation lookups share one HTTP policy

[`apps/common/http_lookup.py`](common/http_lookup.py) is the single policy for the metadata
requests that populate these pickers and the annotation sources: a 5-second connect bound,
a 20-second read bound, and one retry. Both bounds were earned from observed failures, not
guessed:

- `rest.uniprot.org` timing out on connect, and separately dropping a TLS handshake with
  `UNEXPECTED_EOF_WHILE_READING`, while a direct request moments earlier answered in 0.04s.
  One retry recovers the run instead of leaving the dropdown empty -- which reads to the
  user as "this protein has no structures".
- The EBI Proteins API accepting a connection in 0.06s and then never answering. Only the
  read bound covers that, and these calls run in synchronous reactive effects, so an
  unbounded read freezes *every* connected session rather than just the one that made it.

A persistent failure still raises: retrying must never turn an outage into a silently empty
result. File downloads are deliberately excluded and keep a longer timeout, because a large
structure legitimately takes time.

## Docker Compose Services

Defined in [`docker-compose.yml`](../docker-compose.yml):

- `peptide-mapper` -> host port `8001`
- `structure-viz` -> host port `8002`
- `mutation-effect` -> host port `8003`
- `topology-viewer` -> host port `8004`
- `rinalign` -> host port `8005`
- `scop3p-toolkit` -> host port `8000`

## Run Independently

Build one app:

```bash
make peptide-mapper
make structure-viz
make mutation-effect
make topology-viewer
make rinalign
```

Use `make`, not `docker compose build`: the app images are `FROM` a shared base that has
to exist first, and Compose cannot express build ordering. `make <app>` builds the base
if needed. (`docker compose build <app>` works once `make base` has run.)

Start one app:

```bash
docker compose up -d peptide-mapper
docker compose up -d structure-viz
docker compose up -d mutation-effect
docker compose up -d topology-viewer
docker compose up -d rinalign
```

Open in browser:

- Peptide Mapper: `http://localhost:8001`
- Structure Visualisation: `http://localhost:8002`
- Mutation Effect: `http://localhost:8003`
- Topology Viewer: `http://localhost:8004`
- RIN Alignment: `http://localhost:8005`

Stop one app:

```bash
docker compose stop peptide-mapper
docker compose stop structure-viz
docker compose stop mutation-effect
docker compose stop topology-viewer
docker compose stop rinalign
```

## Run All-In-One Mode

Build the single-container toolkit:

```bash
make scop3p-toolkit
```

Start the toolkit:

```bash
docker compose up -d scop3p-toolkit
```

Open in browser:

- toolkit root: `http://localhost:8000`

The toolkit exposes a selector navbar at the root URL. You can also preselect an app with:

- `http://localhost:8000/?app=peptide-mapper`
- `http://localhost:8000/?app=structure-viz`
- `http://localhost:8000/?app=topology-viewer`
- `http://localhost:8000/?app=mutation-effect`
- `http://localhost:8000/?app=rinalign`
- `http://localhost:8000/?app=help`

Stop the toolkit:

```bash
docker compose stop scop3p-toolkit
```

## Run Everything

If you want all independent apps plus the toolkit at the same time:

```bash
docker compose up -d peptide-mapper structure-viz mutation-effect scop3p-toolkit
```

Stop all services:

```bash
docker compose down
```

## Logs

Follow logs for one service:

```bash
docker compose logs -f peptide-mapper
docker compose logs -f structure-viz
docker compose logs -f mutation-effect
docker compose logs -f scop3p-toolkit
```

Each service also writes the same Python logging records to a timestamped file inside
`/var/log/scop3p_toolkit` in the container. Docker Compose mounts that directory to
service-specific host paths:

- `logs/peptide-mapper/`
- `logs/structure-viz/`
- `logs/mutation-effect/`
- `logs/scop3p-toolkit/`

Each mounted directory holds three things:

| File | What it is |
|---|---|
| `scop3p_toolkit_log_<stamp>.log` | everything: the experiment record plus diagnostics |
| `scop3p_toolkit_trail_<stamp>.log` | **the experiment record on its own** |
| `metadata.yml` | context-only FAIR execution metadata, and the paths of the two logs |

`metadata.yml` records app name, session start time, image version/revision/build date,
Python runtime, relevant package versions and available external tools.

### The experiment record

A protocol run is an experiment, so the log is written to read as its record: which
protocol was opened, what went in, which actions were taken in what order, and what each
one produced. It is written to its own file as well as the combined log, so a run can be
handed over as a standalone document rather than grepped out of a log interleaved with
diagnostics.

```
10:20:23 session=0f8098ea step=1  action=open   detail="opened protocol Structure Visualisation"
10:20:32 session=0f8098ea step=2  action=click  detail="Load example"
10:20:32 session=0f8098ea step=3  action=input  detail="UniProtKB accession = P07949"
10:20:32 session=0f8098ea step=4  action=click  detail="Set protein"
10:20:35 session=0f8098ea step=5  action=result detail="34 PDB entries cross-referenced from P07949"
10:20:43 session=0f8098ea step=8  action=click  detail="Fetch PTMs"
10:20:44 session=0f8098ea step=9  action=result detail="23 PTM sites" scop3p=23 uniprot=off
10:20:51 session=0f8098ea step=12 action=click  detail="Run Bio2Byte predictions"
10:20:55 session=0f8098ea step=13 action=result detail="Bio2Byte predictions over 1114 residues" properties=8
10:21:53 session=0f8098ea step=18 action=click  detail="Build RIN"
10:22:01 session=0f8098ea step=19 action=result detail="RIN built for chain A" cutoff=8.0 edges=1319 nodes=283 numbering=sifts-api
```

Two properties make this a record rather than a pile of lines, and both are pinned by
test:

- **Ordering.** Every step carries an incrementing number, so the sequence is explicit
  rather than inferred from timestamps that can tie at millisecond resolution.
- **Attribution.** Every step carries a short session id. One process serves many browser
  sessions concurrently -- the same user with two tabs open is enough -- and without a
  discriminator two interleaved runs are indistinguishable.

The vocabulary is a closed set, so the trail can be parsed and so protocols do not each
invent their own wording: `open`, `input`, `select`, `click`, `result`, `blocked`,
`failed`, `export`.

The portal dispatcher sits outside any Shiny session, so it has no per-user trail of its
own; switching protocol is logged there as `event=navbar_click` with the requested and
selected app. Each protocol then opens its own record with an `open` step.

### Levels

| Level | What goes there |
|---|---|
| `DEBUG` | third-party output, captured rather than printed (see below) |
| `INFO` | the experiment record, and each handler's own start/finish detail |
| `WARNING` | degraded but continuing: a blocked action, an upstream that failed while an optional source stayed off |
| `ERROR` | the action failed. The traceback goes to the app logger; the record gets a one-line step saying where in the sequence it happened |

A failure is never recorded at `INFO`. "0 PDB entries" and "the lookup failed" must not
read the same: the first is an answer about the protein, the second is a broken run, and
recording it as a result would quietly turn an outage into a finding.

`SCOP3P_LOG_LEVEL` quiets diagnostics. It does not touch the record -- the trail logger
pins its own level, so raising the threshold to `WARNING` cannot lose the experiment.

### Third-party output

`b2bTools` reports progress with `print` -- 215 call sites -- and its dependencies warn on
every prediction: scikit-learn's `InconsistentVersionWarning` fires once per unpickled
model, so three predictions produced three copies of a three-line warning. Worse,
`warnings` go straight to stderr and so never reached the log file at all.

Every b2bTools call now runs inside `quiet_third_party`, which routes both stdout and
warnings into the log at `DEBUG`. Measured over three predictions, console output went
from 20 lines to 1. Nothing is discarded -- run with `SCOP3P_LOG_LEVEL=DEBUG` to get it
back when a prediction misbehaves.

There are no `print` statements in `apps/`, and a test walks the AST of every module to
keep it that way. Output that goes through the logger carries a level, a timestamp and a
logger name; a `print` carries none of those, cannot be silenced, and never reaches the
log file.

Override the log location inside a container with `SCOP3P_LOG_DIR` if needed:

```bash
docker run --rm \
  -e SCOP3P_LOG_DIR=/var/log/scop3p_toolkit \
  -v "$(pwd)/logs/peptide-mapper:/var/log/scop3p_toolkit" \
  -p 8001:8000 \
  bio2byte/peptide-mapper:0.1.0
```

## Run Locally Without Docker

`shiny run` puts the app file's own directory on `sys.path`, not `apps/`, so
`from common... import` needs help. `apps/topology_viewer/app.py` and
`apps/rinalign/app.py` add `apps/` themselves and run without the prefix; the three
older apps do not, so pass it:

```bash
PYTHONPATH=apps shiny run --reload --port 8050 apps/topology_viewer/app.py
```

```bash
PYTHONPATH=apps shiny run --reload --port 8052 apps/peptide_mapper/app.py
```

The all-in-one portal is an ASGI app rather than a Shiny app, so it runs under uvicorn:

```bash
PYTHONPATH=apps python -m uvicorn portal.main:app --port 8000
```

## Where the Topology Package Lives

The Topology Viewer's science lives in `notebooks/topology_viewer/topology/`, not under
`apps/`. That is deliberate: it is the single source of truth shared with the Voila
notebook and with `notebooks/topology_viewer/test_topology.py`, and moving it would
break both.

`apps/common/topology_bridge.py` locates it, trying in order:

1. `$SCOP3P_TOPOLOGY_PATH`, if set
2. `/opt/scop3p/topology_viewer`, where `docker/Dockerfile.base` copies it
3. `<repo>/notebooks/topology_viewer`, for a source checkout

The container path sits outside `/apps` so that a `-v ./apps:/apps` dev mount cannot
mask it, and it is deliberately **not** added to `PYTHONPATH`: the bridge must take the
same resolution path in the image as it does locally and in CI, otherwise a broken
bridge would ship green. If nothing resolves, the app renders a card naming every path
it tried instead of failing with an import traceback — which matters because
`apps/portal/main.py` imports every app at module scope, so an exception there would
take down all five.

Note for image builders: `.dockerignore` must not exclude `notebooks/`.

## Image Layout

Two layers, deliberately:

- **`docker/Dockerfile.base`** builds `bio2byte/scop3p-base`: the Python virtualenv,
  `TM-align`, `hmmer`, the topology package and all of `apps/`. It declares no `CMD`
  and selects no app, so it is not runnable and cannot be mistaken for an app image.
- **`docker/apps/<app>.Dockerfile`** builds one app image `FROM` that base, adding only
  `SCOP3P_APP_NAME` and a `CMD`.

This replaced a single multi-stage Dockerfile whose leaf stages were the apps. That
worked, but which app you got from `docker build` depended on which stage happened to
be last, so adding a stage in the wrong place silently published the wrong app. Now the
Dockerfile you name is the app you get.

The base must be built before any app image, which `docker compose` cannot express, so
use the Makefile:

```bash
make base
```

```bash
make apps
```

```bash
make sizes
```

Because every app image is the base plus two metadata lines, they share all layers: six
images cost essentially one image's worth of disk, and each app builds in about a second
once the base exists.

## Security Posture

`make scan` runs Trivy against the base, which is where the whole surface lives.

The base is `python:3.12-slim` **pinned by digest** so scans are reproducible and an
upstream retag cannot change what ships. It is already Debian 13 (trixie) and Python
3.12.14 — the newest Python available, because b2bTools declares `requires_python <3.13`
and its `torch~=2.2.2` pin has no cp313 wheels.

Measured with Trivy, OS packages only:

| | CRITICAL | HIGH |
|---|---|---|
| `python:3.12-slim` as shipped | 4 | 49 |
| `+ apt-get upgrade` | 4 | 13 |
| `+ perl-base removed` (what we ship) | **0** | **9** |

All four CRITICALs were in `perl-base` and none has an upstream fix, so patching alone
could not clear them. Nothing in the image executes perl — `hmmer` is a C binary,
`TM-align` is static, and no app shells out to perl. `t-coffee` was the only perl
consumer and no app invokes it either; it appeared solely in `session_metadata.py`'s
version-probe table, which degrades to a null entry when a binary is absent. So perl and
t-coffee are both gone and `metadata.yml` records `t_coffee` as absent.

`pip` is also removed from the runtime: nothing installs packages at run time, and pip
vendors its own `msgpack` and `urllib3` that cannot be patched independently. That
cleared 8 further findings. `setuptools` stays.

The 9 remaining OS HIGHs (`gzip`, `libacl1`, `ncurses`, `openssl`) have no fix available
upstream. Re-run `make scan` after any change to the apt block, and refresh the pinned
digest deliberately.

### Two findings that need b2bTools

Scanning the Python packages leaves exactly two fixable items, and both are pinned by
b2bTools rather than by us:

| Package | We ship | Fixed in | b2bTools pin | Finding |
|---|---|---|---|---|
| `torch` | 2.2.2+cpu | 2.6.0 | `torch~=2.2.2` | CVE-2025-32434 (CRITICAL) |
| `urllib3` | 1.26.20 | 2.6.0+ | `urllib3~=1.26.20` | 4 × HIGH |

Neither can be resolved without b2bTools relaxing those constraints. On exposure: the
apps never import torch and never call `torch.load` — only b2bTools does, on the model
weights bundled inside its own package, so the CVE is not reachable from user input in
the Shiny apps. Note that the *notebook* version of the mutation-effect protocol does
monkeypatch `torch.load` to `weights_only=False`, which is precisely the unsafe pattern;
that affects the Binder image, not these.

## Scop3P API

The apps talk to the **Scop3P v1 REST API** directly, through
`apps/common/services.py::Scop3PClient`:

| Resource | Endpoint |
|---|---|
| Modifications | `GET /scop3p/api/v1/proteins/{accession}/modifications` |
| Peptides | `GET /scop3p/api/v1/proteins/{accession}/peptides` |

Full spec: <https://iomics.ugent.be/scop3p/api/v1/openapi.json>

The `scop3p` PyPI client is deliberately **not** used. Version 1.1.0, the latest
release, still targets the pre-v1 query-string endpoints
(`/scop3p/api/modifications?accession=...`), which the current deployment no longer
serves.

The notebooks under `notebooks/` target v1 too, through their own small `_scop3p_get`
helper — they stay self-contained rather than importing from `apps/`, but they use the
same field maps and the same content-type guard.

Two properties of this API are worth knowing before debugging it:

- **A retired endpoint returns 200, not 404.** Scop3P has a `GET /scop3p/{catchall}`
  route that serves the single-page-app HTML, so a request to an old path comes back
  `200 OK` with `content-type: text/html`. `raise_for_status()` is therefore useless
  here, and the only symptom is `JSONDecodeError: Expecting value: line 1 column 1`
  from inside `requests`. `Scop3PClient._get_json` checks the content type for exactly
  this reason and raises `Scop3PApiError` naming the URL and pointing at the OpenAPI
  document.
- **v1 renamed every field to snake_case** and returns a bare list rather than a dict.
  `Scop3PClient` maps them back to the column names the toolkit has always used
  (`position`, `residue`, `name`, `peptideSequence`, `peptideStart`, ...), so nothing
  downstream had to change. `modifiedResidue` no longer exists in the peptides payload
  and is recovered from the sequence: the modification position is 1-based within the
  peptide, and `peptide_start + position - 1 == uniprot_position`.

An accession Scop3P does not cover returns `200` with `[]`. That is not an error, and
it is indistinguishable from an unknown accession — so an empty frame means "no
coverage", while a transport or endpoint problem raises.

## Dependencies

`requirements-shiny.txt` is the app runtime. `requirements-biophysics.txt` carries
b2bTools and the Scop3P API client, and pins **CPU-only PyTorch**:

- b2bTools runs inference on CPU, so the default PyPI `torch` wheel is waste here — it
  is CUDA-enabled and pulls eleven `nvidia-*` packages plus `triton`.
  `torch==2.2.2+cpu` from `https://download.pytorch.org/whl/cpu` pulls no CUDA runtime
  at all. Measured by installing `requirements-biophysics.txt` into a fresh venv on
  `python:3.12-slim`, linux/amd64:

  | | virtualenv | CUDA packages |
  |---|---|---|
  | default `torch` | 4.90 GB | 12 |
  | `torch==2.2.2+cpu` | 1.30 GB | 0 |

  That is 3.59 GB (73%) smaller, and about 2.4 GB less to download. The venv is the
  dominant layer of the runtime image, so the saving carries through almost
  one-to-one: the `scop3p-toolkit` image is 2.2 GB uncompressed (454 MB compressed)
  where it would otherwise be near 5.8 GB.
- The pin carries environment markers so a plain `pip install -r
  requirements-biophysics.txt` still works off linux/amd64; the PyPI wheel is already
  CPU-only on macOS.
- The torch version is constrained, not chosen: b2bTools requires `torch~=2.2.2` on
  Python 3.12, the image's base. `torch==1.13.1+cpu` has no cp312 wheel, and moving the
  base to Python 3.10 to get it would also break `apps/common/models.py`, which uses
  `enum.StrEnum` (3.11+).

## Adding Another App

1. `apps/<name>/app.py` exposing a module-level `app = App(content_ui, server)`, with
   the reusable logic in `apps/common/<name>.py` so it can be unit-tested.
2. Register it in `apps/portal/main.py`: one import plus one `APP_OPTIONS` entry
   (`key -> (label, Font Awesome class, app)`). Without this it is unreachable in the
   published all-in-one image.
3. Add `docker/apps/<name>.Dockerfile` — four lines: an `ARG BASE_IMAGE`, the `FROM`,
   `ENV SCOP3P_APP_NAME`, and the `CMD`. Copy any existing one. Add the name to `APPS`
   in the `Makefile`.
4. Add a `docker-compose.yml` service on the next free host port.
5. Add a `docs/use-cases/<name>.md` parity spec.
6. Add unit tests plus a construct assertion in `tests/integration/test_app_smoke.py`.
   The suite runs offline; monkeypatch HTTP.
7. Add a card to `apps/help/app.py`'s `PROTOCOLS` tuple. A smoke test asserts that the
   Help page and `APP_OPTIONS` describe exactly the same set of apps, so a new app
   without a Help entry fails CI.
8. Update this file, the root `README.md`, and `docs/PROTOCOL_STATUS.md`.

### One Shiny detail worth knowing

Outputs inside a `nav_panel` that is **not** the initially-active tab must be declared
`@output(suspend_when_hidden=False)`. Shiny decides suspension from the client-reported
`.clientdata_output_<id>_hidden` value, and `Session._is_hidden()` treats "never
reported" as hidden — so such an output is suspended at page load and is never woken
when the user opens its tab. It sits at "recalculating" forever with nothing logged
anywhere. This affected every tab after the first in `mutation_effect` and
`structure_viz`; both are fixed.
