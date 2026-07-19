# Genome Firewall

Research prototype for the Hack-Nation Genome Firewall challenge.

The goal is to predict whether antibiotics are likely to fail, likely to work, or require a no-call from a reconstructed, quality-checked bacterial FASTA genome. This is defensive decision support only: every result must be confirmed by standard laboratory testing, and the project must not design, modify, optimize, or suggest changes to organisms.

## Problem Scope

Input:

- One reconstructed bacterial genome in FASTA format.
- Initial supported species: `Escherichia coli`.

Output per antibiotic:

- `likely_to_fail`, `likely_to_work`, or `no_call`.
- Calibrated confidence score.
- Evidence category: known resistance gene/mutation, statistical association only, or no known resistance signal.
- Mandatory standard-lab-confirmation warning.

## Dataset Research

The organizers did not provide a fixed challenge dataset, so we need a public fallback.

### Primary Choice: BV-BRC

BV-BRC is the best match for this challenge because it provides genome-level AMR phenotype records, including resistance phenotype and/or MIC values measured by antibiogram/AMR panel assays. BV-BRC also includes computational AMR predictions, so our scripts explicitly filter for `evidence == "Laboratory Method"` and exclude computational predictions from labels.

Current local BV-BRC pull:

- Species filter: `Escherichia coli`
- Raw lab-evidence AMR rows fetched: `20,000`
- Usable cleaned labeled rows: `2,762`
- Unique genomes with cleaned labels: `1,105`
- Antibiotics with most usable labels:
  - `ampicillin`: 1,063
  - `amoxicillin`: 360
  - `ciprofloxacin`: 135
  - `cefotaxime`: 105
  - `gentamicin`: 102
  - `chloramphenicol`: 98
  - `tetracycline`: 94

Useful local files:

- `data/raw/bvbrc_ecoli_lab_amr.jsonl`
- `data/processed/bvbrc_ecoli_labels.csv`
- `reports/metrics/bvbrc_ecoli_label_counts.csv`

### Annotation Tools And Databases

These are useful for genome-to-feature conversion, but they are not by themselves complete supervised datasets because they do not provide matching lab susceptibility labels for each genome.

- AMRFinderPlus: NCBI tool/database for acquired AMR genes and resistance-associated point mutations from nucleotide and/or protein sequence. This is the default challenge annotation tool.
- ResFinder: identifies acquired genes and chromosomal mutations mediating AMR in total or partial bacterial DNA sequence.
- CARD/RGI: predicts resistomes from nucleotide or protein data using CARD homology/SNP models.
- ABRicate: screens contigs against databases such as NCBI, CARD, ARG-ANNOT, ResFinder, MEGARes, and others.
- cAMRah: workflow that runs multiple AMR tools/databases together, including AMRFinderPlus, ResFinder, RGI/CARD, ABRicate/NCBI, ABRicate/ARG-ANNOT, and BV-BRC AMR detection.

### Other Dataset Options

- NCBI Pathogen Detection / AST Browser: useful secondary source because it links antibiotic susceptibility test data with NCBI Pathogen Detection isolates and AMRFinderPlus results. It is a good fallback if BV-BRC sequence retrieval becomes limiting.
- Kaggle mirrors: acceptable only for tutorials or demos when license, original source, and version are clear. A Kaggle copy alone is not a trustworthy benchmark for this challenge.

## Current Implementation Plan

1. Build a cleaned BV-BRC `E. coli` lab-label table.
2. Select a quality-filtered starter cohort with 3-5 antibiotics and both resistant/susceptible examples.
3. Download FASTA files for selected genomes through the BV-BRC API.
4. When `amrfinder` is installed, run AMRFinderPlus on each FASTA and parse TSV output into binary features.
5. If `amrfinder` is still unavailable, use a temporary k-mer baseline only for plumbing/demo validation, clearly marked as non-final.
6. Split by genetic group using BV-BRC `cgmlst_hc100` when available; later replace/augment with `mash`, `fastANI`, or `mmseqs2`.
7. Train one calibrated logistic-regression model per antibiotic.
8. Add no-call logic for low confidence, poor genome quality, conflicting evidence, missing target gate, or out-of-distribution metadata.
9. Evaluate balanced accuracy, resistant/susceptible recall, F1, AUROC, PR-AUC, Brier score, no-call rate, and reliability plots.
10. Build a Streamlit demo for FASTA upload, prediction report, evidence display, and safety warning.

## Major Improvement Areas

The local cohort has now been synchronized against the manually downloaded FASTAs. The cleaned BV-BRC label table has `2,762` usable lab-labeled rows across `1,105` genomes; the current local cohort contains `497` quality-filtered FASTA-backed genomes and `1,136` lab-label rows across six antibiotics. The AMRFinderPlus TSV generation is still catching up, so the next model retrain should happen after those annotations finish.

Priority improvements:

- Keep the cohort synchronized from existing FASTAs using `scripts/prepare_ecoli_cohort.py --use-existing-fasta`.
- Retain auto-selected balanced antibiotics instead of hard-coding a small list.
- Run AMRFinderPlus over the larger FASTA set in WSL.
- Retrain the AMRFinderPlus models and regenerate reports.
- Improve split quality with `mash`, `fastANI`, or `mmseqs2` clustering instead of relying only on BV-BRC `cgmlst_hc100`.
- Expand target-presence verification beyond BV-BRC annotation markers by supporting Bakta/Prokka GFF inputs for uploaded genomes.
- Improve no-call calibration with a larger calibration split.
- Add more drug-specific evidence logic, so explanations emphasize resistance mechanisms relevant to each antibiotic.

## Strong-Submission Readiness

The challenge's strong-submission criteria are addressed in the design, but several still need stronger evidence before submission:

| Criterion | Current status | Remaining work |
| --- | --- | --- |
| One species and a few antibiotics | E. coli scope is enforced; three antibiotics currently have trained models and four are configured for target checks | Finish the larger cohort and retain only drugs with both reliable classes and enough lab labels |
| Calibrated confidence and no-call | Sigmoid calibration, Brier score, reliability plots, and configurable no-call thresholds are implemented | Recalibrate with a larger validation/calibration split and report calibration error; current test sets have only seven rows per drug |
| Honest generalization | `cgmlst_hc100` groups are kept on one side of the split; training now writes a split-audit CSV with overlap counts | Add per-group/generalization reporting; sequence-similarity clustering would be stronger |
| Honest explanations | AMRFinderPlus hits and target-gate matches are shown in reports; AMRFinderPlus report categories are now drug-specific | Keep expanding curated relevant AMR classes and mutation interpretation per antibiotic |
| Molecular-target compatibility | Target markers are checked from BV-BRC annotations for the current metadata genomes | Refresh target annotations after the cohort is expanded and add an annotation path for uploaded FASTAs |
| Defensive scope and oversight | E. coli-only scope, lab-confirmation warning, and no organism-design functionality are present | Keep these limits visible in the demo and presentation |

Current data-processing warning: `data/raw/fasta/` contains `500` FASTA files and the synchronized cohort metadata now contains `497` quality-filtered genomes. AMRFinderPlus TSVs are still incomplete while the WSL annotation run continues, so only genomes with completed TSVs can be used in the next AMRFinderPlus feature build.

## Model Architecture

The model layer is implemented under `src/genome_firewall/models/`.

Current behavior:

- Merges labels, genome metadata, and a feature matrix by `genome_id`.
- Trains one binary model per antibiotic.
- Uses `cgmlst_hc100` as the grouped split key when available, falling back to `genome_id`.
- Writes split-audit files such as `reports/metrics/amrfinder_split_audit.csv` with train/test group counts and any group overlap.
- Skips antibiotics with insufficient class balance.
- Uses balanced logistic regression wrapped in sigmoid calibration.
- Saves one model artifact per antibiotic plus a manifest.
- Converts resistant probability into:
  - `likely_to_fail` when `P(resistant) >= 0.70`
  - `likely_to_work` when `P(resistant) <= 0.30`
  - `no_call` otherwise
- These thresholds are configurable with `--likely-to-work-max` and `--likely-to-fail-min`.
- Emits the mandatory safety warning in prediction outputs.

### Grouped Split

Grouped split means related genomes are kept together when splitting train and test data. We use BV-BRC `cgmlst_hc100` as the current group key when available.

Why this matters:

- Bacterial datasets often contain nearly identical strains.
- A random row split can put near-duplicate genomes in both train and test.
- That inflates scores because the model can memorize strain lineage instead of learning resistance signals.
- A grouped split is stricter: if one genome from a genetic group is in test, related genomes from that group should not appear in train.

Current limitation: `cgmlst_hc100` is a useful BV-BRC grouping field, but a stronger split should be built with sequence similarity tools such as `mash`, `fastANI`, or `mmseqs2`.

### Threshold Tuning

The default decision thresholds are conservative:

```text
P(resistant) >= 0.70 -> likely_to_fail
P(resistant) <= 0.30 -> likely_to_work
otherwise            -> no_call
```

These can and should be tuned on a calibration set. Tradeoff:

- Wider no-call band, for example `0.20/0.80`: fewer wrong confident calls, more no-calls.
- Narrower no-call band, for example `0.40/0.60`: more predictions, higher risk of confident errors.

Example:

```bash
conda run -n genome python scripts/run_amrfinder_pipeline.py --likely-to-work-max 0.25 --likely-to-fail-min 0.75
```

Current model files:

- Training package: `src/genome_firewall/models/training.py`
- Decision rules: `src/genome_firewall/models/decisions.py`
- Metrics: `src/genome_firewall/models/metrics.py`
- Prediction helper: `src/genome_firewall/models/predict.py`
- Train command wrapper: `scripts/train_models.py`
- Prediction command wrapper: `scripts/predict_from_features.py`

Model outputs:

- Manifest: `models/genome_firewall/manifest.json`
- Per-antibiotic models: `models/genome_firewall/*.joblib`
- Held-out metrics: `reports/metrics/model_metrics.csv`
- Held-out prediction report: `reports/metrics/heldout_predictions.csv`
- Full cohort prediction report: `reports/metrics/cohort_predictions.csv`
- Model card: `reports/model_card.md`
- Reliability plots: `reports/figures/reliability_*.png`
- Called-decision matrix plots: `reports/figures/called_confusion_*.png`

## Demo App

The Streamlit app is implemented in `app.py`.

Current app capabilities:

- Default to the AMRFinderPlus gene/mutation model.
- Review AMRFinderPlus cohort predictions with drug-relevant supporting AMR hits.
- Switch to the temporary k-mer model only for upload-flow demonstration.
- Show model metrics.
- Show target-gate metadata for supported antibiotics.
- Show the model card.
- Always display standard-lab-confirmation language.

Run locally:

```bash
conda run -n genome streamlit run app.py --server.port 8501
```

Then open `http://localhost:8501`.

Live FASTA upload is disabled for the AMRFinderPlus model inside the Windows Streamlit app because AMRFinderPlus runs in WSL. To score a new FASTA with the AMRFinderPlus model, add it to `data/raw/fasta/`, run AMRFinderPlus from WSL, then rerun the AMRFinderPlus pipeline below.

## Drug Compatibility Gate

The problem statement requires a deterministic drug compatibility gate based on the presence or absence of the drug's molecular target.

Current implementation:

- Configured in `configs/antibiotics.yaml`.
- Uses curated, citation-backed target markers for the configured antibiotics.
- Fetches general BV-BRC genome annotations with `scripts/fetch_bvbrc_genome_features.py`.
- Builds `data/interim/target_presence.csv` with `scripts/build_target_presence.py`.
- Applied in `src/genome_firewall/models/target_gate.py`.
- Enforced in `scripts/predict_amrfinder_models.py`.
- Adds `target_gate_status`, `target_gate_reason`, `target_present`, `matched_target_markers`, and `molecular_targets` to AMRFinderPlus prediction reports.
- Forces predictions to `no_call` if the species or antibiotic is outside configured scope.
- If target presence is absent or unknown, blocks a confident `likely_to_work` call and returns `no_call`; a `likely_to_fail` call can still be reported because resistance evidence can remain clinically relevant even when target annotation is incomplete.

Configured target classes:

- `cefotaxime`, `ceftazidime`, `ampicillin`, and `amoxicillin`: penicillin-binding proteins.
- `gentamicin`: 30S ribosomal subunit / 16S rRNA A-site markers.
- `chloramphenicol`: 50S ribosomal subunit / 23S rRNA peptidyl-transferase center markers.

The current target marker YAML cites public references, including NCBI Bookshelf pages and peer-reviewed review articles. The markers are class-level compatibility checks, not resistance mechanisms and not proof of susceptibility.

Current target-module outputs for the synchronized 497-genome cohort:

- BV-BRC genome feature rows fetched: `2,657,869`
- Target-presence rows built: `2,982`
- Target status: `496` genomes have configured target markers present for each of `ampicillin`, `amoxicillin`, `cefotaxime`, `ceftazidime`, `chloramphenicol`, and `gentamicin`; `1` genome has unknown annotation status for each drug

Current limitation:

- This target implementation works for genomes that already have BV-BRC `genome_feature` annotations.
- Live uploaded FASTA still needs a local annotation step, such as Bakta or Prokka, before the same target gate can be applied.
- Product/gene names vary across annotation systems, so marker lists should be reviewed and expanded as more antibiotics are added.

## Current Environment Status

The active conda env is expected to be named `genome`.

Python package requirements are listed in `requirements.txt`.

Verified installed:

- `pandas`, `numpy`, `scipy`, `scikit-learn`
- `matplotlib`, `seaborn`, `plotly`
- `pyarrow`, `joblib`, `tqdm`, `pyyaml`, `requests`
- `streamlit`, `gradio`
- `shap`
- `pypdf`, `pymupdf`, `pdfplumber`

Missing from the current env check:

- `biopython` / import name `Bio`
- `reportlab`
- `amrfinder`
- `mmseqs`
- `mash`
- `fastANI` / `fastani`

Retry install command:

```bash
conda install -n genome -c conda-forge -c bioconda biopython reportlab ncbi-amrfinderplus mmseqs2 mash fastani -y
```

Windows remains the main Python/modeling environment. If Windows conda cannot solve/install the bio CLI tools, use WSL only for AMRFinderPlus, MMseqs2, Mash, and FastANI.

Current WSL note: `wsl --list --verbose` only showed `docker-desktop`. That is Docker Desktop's internal distro, not the right place to maintain AMRFinderPlus. Install/use a normal Ubuntu or Debian WSL distro for bio CLI tools.

## Commands

Check installed Python modules and CLI tools:

```bash
conda run -n genome python scripts/check_environment.py
```

Print exact WSL AMRFinderPlus commands for this project path:

```bash
conda run -n genome python scripts/print_wsl_amrfinder_commands.py
```

Build the cleaned BV-BRC label table:

```bash
conda run -n genome python scripts/build_bvbrc_starter_dataset.py
```

Prepare a local starter cohort and FASTA files:

```bash
conda run -n genome python scripts/prepare_ecoli_cohort.py --max-genomes 120 --target-genomes 50
```

Prepare a larger cohort for a stronger model:

```bash
conda run -n genome python scripts/prepare_ecoli_cohort.py --max-genomes 1200 --target-genomes 500 --top-antibiotics 6
```

This can download several GB of FASTA data and makes the WSL AMRFinderPlus step much longer. For a faster middle ground:

```bash
conda run -n genome python scripts/prepare_ecoli_cohort.py --max-genomes 700 --target-genomes 250 --top-antibiotics 5
```

Use already downloaded browser FASTAs without redownloading:

```bash
conda run -n genome python scripts/prepare_ecoli_cohort.py --use-existing-fasta --skip-fasta --max-genomes 1200 --target-genomes 500 --top-antibiotics 6
```

Skip FASTA download when only refreshing metadata/labels:

```bash
conda run -n genome python scripts/prepare_ecoli_cohort.py --max-genomes 120 --target-genomes 50 --skip-fasta
```

Run AMRFinderPlus annotations when `amrfinder` is available:

```bash
conda run -n genome python scripts/run_amrfinder_annotations.py
```

Recommended Windows plus WSL workflow:

1. Keep using Windows conda env `genome` for all Python scripts, model training, reports, and the app.
2. Install only AMRFinderPlus inside an Ubuntu/Debian WSL environment.
3. From WSL, write AMRFinderPlus TSVs directly into the shared Windows project folder at `data/interim/amrfinder/`.
4. Return to Windows and run `scripts/build_amrfinder_features.py`, then `scripts/train_models.py`.

Minimal WSL install/run shape:

```bash
micromamba create -n amr -c conda-forge -c bioconda ncbi-amrfinderplus -y
micromamba activate amr
amrfinder_update
cd /mnt/c/Users/ChandraSekhara/Desktop/Hacknation
mkdir -p data/interim/amrfinder
for f in data/raw/fasta/*.fna; do
  id="$(basename "$f" .fna)"
  if [ ! -s "data/interim/amrfinder/${id}.tsv" ]; then
    amrfinder -n "$f" -O Escherichia --plus -o "data/interim/amrfinder/${id}.tsv"
  fi
done
```

Build AMRFinderPlus gene/mutation features from the annotation TSVs:

```bash
conda run -n genome python scripts/build_amrfinder_features.py
```

Build an interpretable reduced AMRFinderPlus training matrix:

```bash
conda run -n genome python scripts/build_refined_amrfinder_features.py
```

This keeps the raw AMRFinderPlus evidence table unchanged and writes a separate
`data/interim/ecoli_amrfinder_refined_features.csv` matrix with `refined_`
columns. It removes constant/singleton source features, excludes AMR classes
outside the configured drug scope, collapses allele-level markers into
class/subclass/family features, and writes an audit catalogue plus manifest.
Use it with the generic trainer as follows:

```bash
conda run -n genome python scripts/train_models.py \
  --features data/interim/ecoli_amrfinder_refined_features.csv \
  --feature-prefix refined_ \
  --feature-type amrfinderplus_refined_biological_aggregates
```

The reduced CSV is useful for exploratory training and pipeline wiring. For a
final held-out evaluation, fit any prevalence or label-driven selection inside
each grouped training fold rather than once across the whole cohort.

Fetch BV-BRC genome annotations for target compatibility checks:

```bash
conda run --no-capture-output -n genome python scripts/fetch_bvbrc_genome_features.py
```

Build the target-presence table:

```bash
conda run -n genome python scripts/build_target_presence.py
```

Train the challenge-aligned model from AMRFinderPlus features:

```bash
conda run -n genome python scripts/train_models.py --features data/interim/ecoli_amrfinder_features.csv --feature-prefix amr_ --feature-type amrfinderplus_gene_mutation
```

Preferred one-command AMRFinderPlus model build:

```bash
conda run -n genome python scripts/train_amrfinder_models.py
```

Preferred full AMRFinderPlus pipeline after TSVs exist:

```bash
conda run -n genome python scripts/run_amrfinder_pipeline.py
```

Tune no-call thresholds:

```bash
conda run -n genome python scripts/run_amrfinder_pipeline.py --likely-to-work-max 0.25 --likely-to-fail-min 0.75
```

Generate AMRFinderPlus cohort predictions with drug-relevant supporting AMR hits:

```bash
conda run -n genome python scripts/predict_amrfinder_models.py
```

Generate evaluation plots and model card:

```bash
conda run -n genome python scripts/generate_evaluation_report.py
```

Generate AMRFinderPlus-specific evaluation plots and model card:

```bash
conda run -n genome python scripts/generate_evaluation_report.py --metrics reports/metrics/amrfinder_model_metrics.csv --predictions reports/metrics/amrfinder_heldout_predictions.csv --model-card reports/amrfinder_model_card.md --prefix amrfinder --feature-source "AMRFinderPlus AMR gene and mutation presence features"
```

Compare classical models on refined AMRFinderPlus aggregate features:

```bash
conda run -n genome python scripts/compare_classical_models.py
```

Run feature-group ablations with one fixed calibrated logistic model:

```bash
conda run -n genome python scripts/ablate_refined_features.py
```

Temporary plumbing baseline while AMRFinderPlus is unavailable:

```bash
conda run -n genome python scripts/build_kmer_features.py --k 4
conda run -n genome python scripts/train_models.py
conda run -n genome python scripts/predict_from_features.py
```

The k-mer baseline is not the final scientific approach. It exists to validate project wiring, splitting, metrics, and model serialization until AMRFinderPlus features are available.

## Current Local State

Current prepared cohort:

- Quality-filtered FASTA-backed genomes currently represented in metadata: `497`
- FASTA files currently present on disk: `500`
- AMRFinderPlus TSVs currently present on disk: `499`
- Cohort labels: `1,136`
- Cohort labels file: `data/processed/ecoli_cohort_labels.csv`
- Cohort metadata file: `data/processed/ecoli_genome_metadata.csv`
- Skipped FASTA log: `reports/metrics/ecoli_skipped_fasta.csv`

The cohort script resumes downloads, skips genomes with missing sequence records, filters for good-quality genomes, and keeps filling until the target number of FASTA-backed genomes is reached.

Current cohort label balance:

| Antibiotic | Resistant | Susceptible | Modeling status |
| --- | ---: | ---: | --- |
| amoxicillin | 187 | 134 | usable, weak separation |
| ampicillin | 295 | 188 | usable |
| cefotaxime | 51 | 37 | usable |
| ceftazidime | 47 | 26 | usable |
| chloramphenicol | 60 | 27 | usable |
| gentamicin | 37 | 47 | usable |

Temporary k-mer baseline outputs:

- Feature matrix: `data/interim/ecoli_kmer_features.csv`
- Metrics: `reports/metrics/model_metrics.csv`
- Held-out predictions: `reports/metrics/heldout_predictions.csv`
- Full cohort predictions: `reports/metrics/cohort_predictions.csv`
- Models: `models/genome_firewall/`

Current package-backed temporary baseline metrics:

| Antibiotic | Balanced accuracy | Resistant recall | Susceptible recall | F1 | AUROC | PR-AUC |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cefotaxime | 0.583 | 0.667 | 0.500 | 0.571 | 0.667 | 0.810 |
| chloramphenicol | 0.500 | 1.000 | 0.000 | 0.727 | 0.083 | 0.451 |
| gentamicin | 0.500 | 1.000 | 0.000 | 0.600 | 0.833 | 0.867 |

Current no-call behavior on held-out rows:

| Antibiotic | No-call rate | Called rows | Called accuracy |
| --- | ---: | ---: | ---: |
| cefotaxime | 1.000 | 0 | unavailable |
| chloramphenicol | 0.429 | 4 | 0.250 |
| gentamicin | 1.000 | 0 | unavailable |

These scores are not submission-quality. They use capped 4-mer frequencies over the first 500,000 bases per genome only to validate the data/model/reporting plumbing. The challenge-aligned model should use AMRFinderPlus-derived genes and mutations as features.

Current demo/reporting outputs:

- Streamlit app: `app.py`
- Antibiotic target-gate scaffold: `configs/antibiotics.yaml`
- Model card: `reports/model_card.md`
- Decision summary: `reports/metrics/decision_summary.csv`
- Reliability plots: `reports/figures/reliability_cefotaxime.png`, `reports/figures/reliability_chloramphenicol.png`, `reports/figures/reliability_gentamicin.png`
- Called-decision confusion plot: `reports/figures/called_confusion_chloramphenicol.png`

Current AMRFinderPlus model path status:

- AMRFinderPlus TSVs for the current local FASTA set are in `data/interim/amrfinder/`.
- AMRFinderPlus feature matrix exists at `data/interim/ecoli_amrfinder_features.csv`.
- Refined AMRFinderPlus feature matrix exists at `data/interim/ecoli_amrfinder_refined_features.csv`.
- AMRFinderPlus evidence table exists at `data/interim/ecoli_amrfinder_evidence.csv`.
- Full AMRFinderPlus hit table exists at `data/interim/ecoli_amrfinder_all_hits.csv`.
- AMRFinderPlus models are saved under `models/amrfinder/`.
- Refined AMRFinderPlus models are saved under `models/amrfinder_refined/`.
- Classical model comparison artifacts are saved under `models/classical_comparison/`.
- AMRFinderPlus metrics are saved at `reports/metrics/amrfinder_model_metrics.csv`.
- Refined AMRFinderPlus metrics are saved at `reports/metrics/amrfinder_refined_model_metrics.csv`.
- Classical model comparison is saved at `reports/metrics/classical_model_comparison.csv`.
- Best classical model per antibiotic is saved at `reports/metrics/best_classical_models.csv`.
- AMRFinderPlus held-out predictions are saved at `reports/metrics/amrfinder_heldout_predictions.csv`.
- AMRFinderPlus cohort prediction report with drug-relevant supporting hits is saved at `reports/metrics/amrfinder_cohort_predictions.csv`.
- AMRFinderPlus model card is saved at `reports/amrfinder_model_card.md`.
- BV-BRC target annotation table exists at `data/interim/bvbrc_genome_features.csv`.
- Target-presence table exists at `data/interim/target_presence.csv`.
- AMRFinderPlus prediction report now includes target-gate columns: `target_gate_status`, `target_gate_reason`, `target_present`, `matched_target_markers`, and `molecular_targets`.

Current raw AMRFinderPlus logistic-regression metrics after the 497-genome cohort sync:

| Antibiotic | Balanced accuracy | Resistant recall | Susceptible recall | F1 | AUROC | PR-AUC | No-call rate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| amoxicillin | 0.552 | 0.925 | 0.179 | 0.731 | 0.595 | 0.678 | 0.957 |
| ampicillin | 0.931 | 0.947 | 0.914 | 0.947 | 0.946 | 0.972 | 0.043 |
| cefotaxime | 0.742 | 0.667 | 0.818 | 0.571 | 0.848 | 0.667 | 0.714 |
| ceftazidime | 0.742 | 0.769 | 0.714 | 0.800 | 0.857 | 0.940 | 0.700 |
| chloramphenicol | 0.900 | 1.000 | 0.800 | 0.957 | 0.982 | 0.992 | 0.375 |
| gentamicin | 0.923 | 1.000 | 0.846 | 0.917 | 0.930 | 0.854 | 0.667 |

Best classical model per antibiotic on refined AMRFinderPlus aggregate features:

| Antibiotic | Best model | Balanced accuracy | Brier score | No-call rate |
| --- | --- | ---: | ---: | ---: |
| amoxicillin | gradient boosting | 0.558 | 0.229 | 0.880 |
| ampicillin | linear SVC | 0.959 | 0.039 | 0.033 |
| cefotaxime | extra trees | 0.909 | 0.085 | 0.071 |
| ceftazidime | extra trees | 0.819 | 0.165 | 0.400 |
| chloramphenicol | gradient boosting | 0.955 | 0.066 | 0.000 |
| gentamicin | gradient boosting | 0.962 | 0.066 | 0.042 |

These metrics use grouped splits with `cgmlst_hc100`; all current trained antibiotics report `0` overlapping train/test groups. The weaker `amoxicillin` result should be shown honestly as high-resistant-recall but poor susceptible recall.

Feature-engineering ablation with the same calibrated logistic model:

| Antibiotic | Best refined feature group | Balanced accuracy | Interpretation |
| --- | --- | ---: | --- |
| amoxicillin | subclass only | 0.562 | Very small gain; still weak |
| ampicillin | family only | 0.968 | Family-level beta-lactamase/marker aggregates helped most |
| cefotaxime | all refined / family / subclass tied | 0.909 | Mechanism aggregates fixed much of the raw sparse-marker weakness |
| ceftazidime | subclass only / all refined tied | 0.819 | Subclass-level beta-lactam signal helped |
| chloramphenicol | class/family/subclass all tied | 0.955 | Phenicol signal is strong and simple |
| gentamicin | family only / class+family tied | 0.962 | Aminoglycoside family aggregates helped most |

The ablation suggests the useful feature engineering is mainly collapsing sparse allele-level AMRFinder hits into `family` and `subclass` mechanism groups. Broad `summary_only` burden counts are weaker by themselves and should be presented as statistical context, not known biological evidence.

## Safety Notes

- Use only lab-measured AMR labels for supervised training.
- Do not train on BV-BRC computational AMR predictions as labels.
- Do not claim broad species or antibiotic coverage.
- Keep no-call as a first-class outcome.
- Do not present SHAP or feature importance as proof of biological causality.
- Always show: `Research prototype only. Confirm with standard laboratory testing.`
