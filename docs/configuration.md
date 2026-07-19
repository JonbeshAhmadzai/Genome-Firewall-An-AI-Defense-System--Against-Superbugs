# Extending coverage

`targets_config.py` is the single registry for pathogen and drug coverage.
Run the validation command before fetching or training:

```bash
python data_fetchers/validate_config.py
```

## Add a bacterial antibiotic

Add an alias, molecular target genes, annotation markers, and a target-pair
entry. Set `enabled: true` only after enough resistant/susceptible laboratory
labels are available:

```python
ANTIBIOTIC_ALIASES["drug-name"] = "drug-name"
MOLECULAR_TARGETS["drug-name"] = ("target_gene",)
TARGET_ANNOTATION_MARKERS["drug-name"] = {
    "genes": ("target_gene",),
    "product_keywords": ("target protein",),
}
TARGET_DEFINITIONS[("Escherichia coli", "drug-name")] = {"enabled": True}
```

Then select genomes, annotate them, build features/target presence, and run
`data_fetchers/train_cohort.py`. The trained artifact and report are discovered
by the router and Streamlit app from the species/drug names.

## Add a bacterial species

Add one `PATHOGEN_CONFIG` entry with its BV-BRC taxon ID and AMRFinderPlus
organism name, then add the species/drug pair to `TARGET_DEFINITIONS`. The
router selects the organism reader and model path from that configuration.

## Add a virus

Register the virus with `kind: "virus"` and a virus-specific `reader` name.
The current AMRFinderPlus reader intentionally rejects virus routes. To enable
one, implement that reader/backend and its validated feature/target rules first;
then the same router can dispatch the configured virus route without treating it
as a bacterium.
