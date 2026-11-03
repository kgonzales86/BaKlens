# BaKlens

**An interactive bacterial genome explorer built in Python.** Upload a DNA FASTA to inspect sequence composition, contig statistics, and GC skew. Add Bakta annotations to explore genes and RNA features on circular and linear maps.

BaKlens runs locally, needs no API key, and does not upload your sequences to external services. It is a research and learning tool, not a clinical diagnostic system.

## Quick start

Python 3.11 or 3.12 is recommended. From this folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m streamlit run app.py --server.address 127.0.0.1
```

Open the local URL printed by Streamlit (normally http://localhost:8501). The app starts with a **synthetic demonstration**, including illustrative annotations. It is not a real organism or a real Bakta run.

On Windows, activate using `.venv\Scripts\activate`; the core dashboard and JSON import work there. The local Bakta runner targets Linux/macOS; on Windows use WSL or import results from another machine.

## What you can do

- Upload multi-record DNA FASTA, with strict IUPAC DNA validation and duplicate-ID checks.
- Compare contig lengths, GC%, and ambiguous-base percentages.
- Inspect adjustable sliding-window GC content, GC skew, and ambiguous-base tracks.
- View one sequence at a time as a circular composition plot.
- Import Bakta JSON, or execute an installed local Bakta pipeline.
- Browse gene/RNA tracks, select a linear coordinate range, and search annotations.
- Export CSV tables, a Markdown summary, raw annotation JSON, and self-contained interactive HTML plots.

The default limits are **50 MB per uploaded file, 25 million DNA bases and 10,000 sequences**. Each selected sequence is limited to 50,000 analysis windows. Circular composition plots are sampled to about 3,000 windows for display; CSV exports retain all calculated windows. This is intended for bacterial assemblies, not whole human genomes.

## Three ways to use it

### 1. FASTA only — no Bakta installation

Choose **DNA FASTA**, upload `.fna`, `.fa` or `.fasta`, and explore. Composition calculations do not need annotations. `examples/demo.fna` is available for a first run.

### 2. Existing Bakta results

Choose **Bakta JSON** and upload the JSON produced by Bakta. The file must contain `sequences` (including `id` and `nt`) and `features`.

Alternatively, upload FASTA first and expand **Add Bakta annotations**. The app requires identical sequence IDs and bases across the two files; it rejects a mismatch instead of silently attaching the wrong genes. If Bakta renamed your contigs, load its JSON alone, or rerun Bakta using `--keep-contig-headers`.

### 3. Run Bakta from the dashboard — optional setup

Bakta is a separate, computationally heavier dependency with its own database. It is not included in `requirements.txt`. Consult the [official installation and database instructions](https://github.com/oschwengers/bakta#installation) for the supported environment and current requirements.

For a Linux setup using Conda, one approach is a dedicated environment containing both Bakta and this application:

```bash
conda create -n baklens-bakta -c conda-forge -c bioconda python=3.11 bakta
conda activate baklens-bakta
python -m pip install -r requirements.txt
bakta_db download --output "$HOME/bakta-databases" --type light
```

Set `BAKTA_DB` to the **installed database directory containing `version.json`**, not its parent. For example, if the installer created `db-light`:

```bash
export BAKTA_DB="$HOME/bakta-databases/db-light"
python -m streamlit run app.py --server.address 127.0.0.1
```

Check the actual directory printed by the installer. The light database is an option for constrained machines; it does not provide all the annotation capabilities of the full database. Check available disk space before downloading. Database downloads and full annotation are not instantaneous.

In the app, upload FASTA, expand **Add Bakta annotations**, then click **Run local Bakta**. This runs only on the machine hosting the app. No external web service is contacted by BaKlens.

The wrapper uses an argument list (no shell), a temporary working directory, preserved sequence headers, and a one-hour timeout. Bakta runs with its own plotting disabled because BaKlens generates the views. On successful completion, JSON and the process log can be downloaded; other Bakta outputs are removed with temporary files. For complete Bakta outputs, run its CLI separately:

```bash
bakta --db "$BAKTA_DB" --keep-contig-headers --threads 2 \
  --output results/my-genome --prefix my-genome genome.fna
```

Bakta is designed for bacteria and plasmids. Do not interpret the acceptance of a FASTA file as verification of its taxonomic identity. BaKlens does not automatically determine organism or circular topology, or set Bakta's `--complete` flag.

## Methods and interpretation

| Quantity | Definition |
|---|---|
| GC percentage | `100 × (G+C)/(A+C+G+T)` |
| GC skew | `(G−C)/(G+C)` on the submitted strand |
| Ambiguous bases | Valid IUPAC DNA symbols other than A/C/G/T, including N |
| N count | Count of N only; a subset of ambiguous bases |
| N50 | Sequence length at which descending lengths reach at least half the assembly size |
| L50 | Number of those sequences needed to reach that threshold |

- Undefined ratios are missing values, not zero; plots leave gaps.
- GC percentage for the assembly is weighted by canonical bases, not averaged across contigs.
- Windows never cross contig boundaries. The final window may be partial; overlapping windows stop after the first window reaches the end. If step exceeds window size, gaps between windows are intentional.
- Displayed and exported feature coordinates are **1-based and inclusive**. Origin-spanning features marked `edge` by Bakta are split into two graphical segments while retaining the original coordinates in tables.
- A circular plot is a layout, not evidence that a sequence is circular. Contigs are not concatenated into a fictitious chromosome.
- GC/skew anomalies do not establish contamination or replication origins. N50 does not establish assembly correctness.
- Annotation labels describe predictions, not experimental validation. The hypothetical-CDS count matches the phrase `hypothetical protein`; it is not a general measure of annotation accuracy.
- The written summary is generated from statistics using templates. **No generative AI operates inside the app.** AI-assisted development should be described separately from product capabilities.

## Project structure

```text
app.py                  Streamlit interface and workflow
baklens/analysis.py     FASTA validation, statistics, windows and report
baklens/annotation.py   Bakta JSON parsing and optional local runner
baklens/plots.py        Interactive genome and assembly figures
examples/               Reproducible synthetic FASTA and annotation fixtures
tests/                  Calculation, parsing, runner and interface tests
.github/workflows/      GitHub Actions test workflow
```

## Tests

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
```

Tests cover ambiguous bases, undefined denominators, weighted GC, N50, partial windows, invalid input, annotation identity matching, coordinate bounds, origin-spanning features, the runner's command contract, and Streamlit navigation. The local runner test uses a fake subprocess, **not a real Bakta installation**.

Development verification: 33 tests passed on Python 3.12. Streamlit started successfully and its AppTest harness exercised the interface; browser screenshot verification was unavailable in this environment. Full annotation with a downloaded Bakta database has not been verified in the development environment; validate it locally before presenting that part as end-to-end tested. GitHub Actions is configured for Python 3.11 and 3.12 but is only run after you push the repository.

## Publish to GitHub

Create an empty repository named `BaKlens` in your GitHub account, then run from this folder:

```bash
git init
git add .
git commit -m "Build BaKlens genome exploration dashboard"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/BaKlens.git
git push -u origin main
```

Replace `YOUR_USERNAME`. No repository has been created or pushed by this project bundle. Keep real genomes, databases, environment files, and credentials out of the repository. The example fixtures are deliberately synthetic and reproducible with `python examples/generate_demo.py`.

## Credits and scope

BaKlens integrates existing tools; it does not introduce a new annotation algorithm. Bakta already provides circular plots and a viewer. This project's focus is a local dashboard combining adjustable composition analysis, contig comparisons, annotation exploration, and reproducible exports.

- [Bakta](https://github.com/oschwengers/bakta) — optional annotation engine and JSON format.
- [Biopython](https://biopython.org/) — FASTA parsing.
- [Streamlit](https://streamlit.io/), [Plotly](https://plotly.com/python/), [pandas](https://pandas.pydata.org/) and [NumPy](https://numpy.org/) — interface, visualizations, and analysis.

Bakta, its databases, and other dependencies retain their respective licenses and citation requirements. No Bakta database is bundled. BaKlens source is released under the MIT license below; this does not relicense its dependencies.
