# 60-second video scripts

These scripts are written for four speakers. At a natural pace, each speaker
has roughly 12–15 seconds. Keep the product UI, architecture diagram, and one
JSON result visible while speaking.

## Technical video — stack, architecture, implementation

**Speaker 1 — 0:00–0:14**

“Genome Firewall is a defensive genomic decision-support prototype. A user
uploads one reconstructed bacterial FASTA, and the system returns a separate
antibiotic prediction with confidence, evidence, target status, and a safe
no-call option.”

**Speaker 2 — 0:14–0:28**

“The pipeline starts with FASTA quality control. AMRFinderPlus identifies
resistance genes and mutations, and our Python feature builder creates a
reproducible binary feature matrix plus engineered evidence counts. Genome
groups keep near-identical strains out of both train and test.”

**Speaker 3 — 0:28–0:43**

“For each configured species and antibiotic, we train a calibrated logistic
model and compare it with random forest, extra trees, and XGBoost. The target
router selects the correct reader, model artifact, and molecular target rules.
Low confidence or missing target compatibility produces no-call instead of an
overconfident answer.”

**Speaker 4 — 0:43–0:58**

“The backend is FastAPI with a JSON prediction endpoint and a lightweight HTML
interface; Streamlit remains as a fallback demo. An optional fast GPT-4o-mini
layer explains structured results only—it cannot change the model decision.
Every result clearly says to confirm with standard laboratory AST.”

**End card — 0:58–1:00**

“Genome Firewall: genomic evidence first, uncertainty made visible.”

## Team video — introduction and roles

**Speaker 1 — 0:00–0:14**

“Hi, we are the four-person team behind Genome Firewall. We built a defensive
AI system to help interpret genomic resistance evidence faster, while keeping
laboratory confirmation and human oversight at the center.”

**Speaker 2 — 0:14–0:28**

“I worked on the data and bioinformatics pipeline: BV-BRC cohort selection,
FASTA quality control, AMRFinderPlus annotation, evidence extraction, and the
reproducible feature artifacts.”

**Speaker 3 — 0:28–0:43**

“I worked on machine learning and validation: grouped train/test splits,
calibrated per-antibiotic models, model comparisons, confidence thresholds,
target compatibility, and the conservative no-call policy.”

**Speaker 4 — 0:43–0:58**

“I worked on product integration: the FastAPI JSON backend, HTML dashboard,
Streamlit fallback, downloadable reports, and the optional evidence-grounded
LLM explanation. Together, we turned the research pipeline into a demo people
can inspect.”

**End card — 0:58–1:00**

“Genome Firewall—defending treatment decisions with transparent genomic AI.”

## Recording tips

- Show the FASTA upload and QC cards during Speaker 2.
- Show the model comparison table and a `no-call` result during Speaker 3.
- Show the HTML dashboard and JSON response during Speaker 4.
- Keep the phrase “research prototype” visible on screen; do not describe it as
  an autonomous clinical treatment system.
