# Team video presentation — 60 seconds

Use this as a simple four-slide storyboard plus an end card. Keep each slide on screen while the
matching speaker talks. The four speakers can record their sections separately
and combine them with one short end card.

## Slide 1 — Data pipeline (0:00–0:14)

**Visual:** BV-BRC cohort → FASTA files → quality-control check.

**Speaker 1:**

> “Hi, we are the four-person team behind Genome Firewall. I worked on the data
> pipeline: selecting the BV-BRC genomes, organizing the cohort, checking FASTA
> inputs, and making the raw data reproducible for the rest of the team.”

## Slide 2 — Feature engineering (0:14–0:28)

**Visual:** AMRFinderPlus evidence → engineered feature matrix.

**Speaker 2:**

> “I worked on feature engineering. I ran AMRFinderPlus evidence through our
> feature builder, encoded resistance determinants and target signals, and
> created the 115-feature matrix used consistently by every model.”

## Slide 3 — Machine learning and validation (0:28–0:44)

**Visual:** Feature matrix → grouped split → model cards.

**Speaker 3:**

> “I built and evaluated the machine-learning layer: grouped train, calibration,
> and test splits prevent genetic leakage; calibrated antibiotic models produce
> confidence scores; and we compare logistic regression, random forest, extra
> trees, and XGBoost with a conservative no-call policy.”

## Slide 4 — UI and dashboard (0:44–0:58)

**Visual:** HTML dashboard beside a JSON response and model-comparison table.

**Speaker 4:**

> “I built the product layer: the FastAPI JSON backend, HTML dashboard, model
> comparison view, downloadable results, and Streamlit fallback. The optional
> GPT explanation summarizes evidence but cannot change a prediction. Together,
> we made the pipeline easy to inspect and safe to demonstrate.”

## End card — Safety message (0:58–1:00)

**Visual:** Product logo, “Research prototype — confirm with laboratory AST”.

**Speaker 4:**

> “Together, we turned a reproducible genomic pipeline into an inspectable demo.
> Genome Firewall: evidence first, uncertainty visible.”

## Recording checklist

- Show the FASTA upload and QC result during Speaker 1.
- Show AMRFinderPlus evidence and the feature matrix during Speaker 2.
- Show the model-comparison table briefly; do not imply clinical validation.
- Keep “Research prototype — confirm with standard laboratory testing” visible.
- Aim for 125–140 spoken words per minute and leave the final two seconds for
  the end card.
