# Scop3P Shiny Apps

This directory contains the three converted Python Shiny apps and one wrapper portal:

- [`apps/peptide_mapper/app.py`](/Users/adrian/workspace/vub/Scop3P-notebooks/apps/peptide_mapper/app.py): Peptide Mapper
- [`apps/structure_viz/app.py`](/Users/adrian/workspace/vub/Scop3P-notebooks/apps/structure_viz/app.py): Structure Visualisation
- [`apps/mutation_effect/app.py`](/Users/adrian/workspace/vub/Scop3P-notebooks/apps/mutation_effect/app.py): Mutation Effect
- [`apps/portal/main.py`](/Users/adrian/workspace/vub/Scop3P-notebooks/apps/portal/main.py): single-root all-in-one selector

## Apps

### Peptide Mapper
- Loads peptide and modification data from Scop3P
- Filters and selects peptide spans
- Maps selections on AlphaFold structures
- Exports a styled HTML session

### Structure Visualisation
- Fetches PTMs and disease-associated variants
- Renders PDB or AlphaFold structures
- Runs Bio2Byte predictions and 3D coloring
- Builds residue interaction networks
- Runs TM-align comparisons

### Mutation Effect
- Fetches WT UniProt sequence and Scop3P PTMs
- Runs Bio2Byte WT and mutant predictions
- Compares WT vs mutant residue features
- Generates mutation-centric inference summaries

## Docker Compose Services

Defined in [`docker-compose.yml`](/Users/adrian/workspace/vub/Scop3P-notebooks/docker-compose.yml):

- `peptide-mapper` -> host port `8001`
- `structure-viz` -> host port `8002`
- `mutation-effect` -> host port `8003`
- `scop3p-suite` -> host port `8000`

## Run Independently

Build one app:

```bash
docker compose build peptide-mapper
docker compose build structure-viz
docker compose build mutation-effect
```

Start one app:

```bash
docker compose up -d peptide-mapper
docker compose up -d structure-viz
docker compose up -d mutation-effect
```

Open in browser:

- Peptide Mapper: `http://localhost:8001`
- Structure Visualisation: `http://localhost:8002`
- Mutation Effect: `http://localhost:8003`

Stop one app:

```bash
docker compose stop peptide-mapper
docker compose stop structure-viz
docker compose stop mutation-effect
```

## Run All-In-One Mode

Build the single-container suite:

```bash
docker compose build scop3p-suite
```

Start the suite:

```bash
docker compose up -d scop3p-suite
```

Open in browser:

- Suite root: `http://localhost:8000`

The suite exposes a selector navbar at the root URL. You can also preselect an app with:

- `http://localhost:8000/?app=peptide-mapper`
- `http://localhost:8000/?app=structure-viz`
- `http://localhost:8000/?app=mutation-effect`

Stop the suite:

```bash
docker compose stop scop3p-suite
```

## Run Everything

If you want all independent apps plus the suite at the same time:

```bash
docker compose up -d peptide-mapper structure-viz mutation-effect scop3p-suite
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
docker compose logs -f scop3p-suite
```
