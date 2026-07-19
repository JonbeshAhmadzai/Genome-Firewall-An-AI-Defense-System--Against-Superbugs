# Live demo walkthrough script — approximately 75 seconds

Use `data/demo/mock_genome.fna` for this walkthrough. It is synthetic input
created for demonstration and is not a clinical isolate.

## 0:00–0:10 — Open the app

> “This is Genome Firewall, a research prototype for transparent genomic
> antibiotic-resistance screening. The main page is focused on running one
> genome analysis, while the full validation view is available in a separate
> coverage window.”

## 0:10–0:20 — Show coverage

Click **Open full coverage**.

> “Here we can see the current scope: 100 genomes, three trained antibiotic
> models, 115 engineered AMR features, and the held-out model comparison. The
> validation uses genetic-group separation with 60 percent training, 15 percent
> calibration, and 25 percent testing.”

Return to the analysis page.

## 0:20–0:32 — Upload and QC

Upload `mock_genome.fna`, select **Escherichia coli**, and click **Run genome
analysis**.

> “I’ll upload a small synthetic FASTA file. The backend first checks sequence
> quality, reports contig count, assembly size, N50, and ambiguous bases, then
> runs AMRFinderPlus to extract resistance evidence.”

## 0:32–0:48 — Read the prediction cards

> “The system returns a separate result for each configured antibiotic. Each
> card shows the model decision, probability, confidence, target status, and
> the reason for any no-call. A high model probability cannot bypass the target
> compatibility gate, so uncertainty stays visible instead of becoming an
> overconfident recommendation.”

## 0:48–1:00 — Explain AMR evidence

Scroll to **AMR evidence**.

> “This input has no reportable AMRFinderPlus determinant. That does not prove
> susceptibility; it only means no matching resistance gene or mutation was
> detected in this run. The UI makes that limitation explicit.”

## 1:00–1:15 — Optional explanation and close

Enable **Generate optional GPT explanation** and run again, if time allows.

> “The optional GPT layer summarizes the structured QC, prediction, and evidence
> fields. It cannot change the model output. Every result ends with the same
> safety boundary: this is a research prototype, and findings must be confirmed
> with standard laboratory AST.”

## Fast fallback version — approximately 30 seconds

> “Genome Firewall takes a FASTA genome, checks its quality, runs AMRFinderPlus,
> and scores each configured antibiotic with calibrated models. The dashboard
> shows confidence, target compatibility, AMR evidence, and a conservative
> no-call when evidence is insufficient. The coverage window exposes our
> 100-genome validation and model comparison. An optional GPT layer explains
> structured results but cannot change them. This is a research prototype—every
> result requires laboratory confirmation.”
