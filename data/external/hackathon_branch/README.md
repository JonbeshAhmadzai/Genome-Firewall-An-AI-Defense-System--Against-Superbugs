# External data payload from `genome-firewall-hackathon`

Source archive: `/home/becode/Downloads/data.zip`  
SHA-256: `334804e3ca212aed126a51e5353cc99ac5fbff615d9273c7979971bebec4bca1`

This payload is kept separate from the Blueprint 1 cohort. Its processed
labels cover 50 genomes and 195 laboratory-labelled rows; its AMRFinder matrix
covers 53 genome IDs, and its BV-BRC annotation/target tables cover 50 genome
IDs. The archive also contains 500 FASTA files, most of which are not part of
the labelled cohort.

The IDs are distinct from the current `cohort100` selection (only one FASTA ID
overlaps), so labels and features must not be merged by row position. Use the
tables here as a reproducible external benchmark/reference source unless an
explicit cross-cohort ID audit is performed.
