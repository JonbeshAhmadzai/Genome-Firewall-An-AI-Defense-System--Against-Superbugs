# Team video presentation — 60 seconds

Use this as a simple six-slide storyboard. Keep each slide on screen while the
matching speaker talks. The four speakers can record their sections separately
and combine them with one short end card.

## Slide 1 — The team and the problem (0:00–0:10)

**Visual:** Four team names/photos and the Genome Firewall title.

**Speaker 1:**

> “Hi, we are the four-person team behind Genome Firewall. We built a defensive
> genomic AI prototype to help interpret antibiotic-resistance evidence faster,
> while keeping laboratory confirmation and human oversight at the center.”

## Slide 2 — Data and bioinformatics (0:10–0:20)

**Visual:** FASTA upload → QC → AMRFinderPlus → evidence table.

**Speaker 2:**

> “I built the data and bioinformatics pipeline: selecting the BV-BRC cohort,
> checking FASTA quality, running AMRFinderPlus, and creating reproducible
> resistance features and evidence artifacts.”

## Slide 3 — Machine learning (0:20–0:32)

**Visual:** Feature matrix → grouped split → model cards.

**Speaker 3:**

> “I developed the machine-learning and validation layer. We use genetic-group
> separated evaluation, calibrated per-antibiotic models, and comparisons with
> random forest, extra trees, and XGBoost.”

## Slide 4 — Safety and routing (0:32–0:43)

**Visual:** Species + antibiotic → router → model artifact → prediction/no-call.

**Speaker 3 (continued):**

> “The configuration-driven router selects the compatible species, antibiotic,
> model, and molecular targets. Low confidence or incompatible targets produce
> a transparent no-call instead of an overconfident recommendation.”

## Slide 5 — Product integration (0:43–0:55)

**Visual:** HTML dashboard beside a JSON response and model-comparison table.

**Speaker 4:**

> “I integrated the FastAPI JSON backend and lightweight HTML dashboard, with
> Streamlit retained as a fallback. An optional fast GPT-4o-mini explanation
> summarizes structured evidence but cannot change the model decision.”

## Slide 6 — Close and safety message (0:55–1:00)

**Visual:** Product logo, “Research prototype — confirm with laboratory AST”.

**Speaker 4:**

> “Together, we turned a reproducible genomic pipeline into an inspectable demo.
> Genome Firewall: evidence first, uncertainty visible.”

## Recording checklist

- Show one FASTA upload, one QC result, and one prediction card.
- Show the model-comparison table briefly; do not imply clinical validation.
- Keep “Research prototype — confirm with standard laboratory testing” visible.
- Aim for 125–140 spoken words per minute and leave the final two seconds for
  the end card.
