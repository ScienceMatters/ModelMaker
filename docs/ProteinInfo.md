# ProteinInfo

## Purpose

ProteinInfo is the fundamental biological object within ModelMaker.

It stores information that is intrinsic to a protein and independent of any
genome editing strategy.

ProteinInfo should answer questions such as:

- What protein is this?
- What organism does it come from?
- What is its amino acid sequence?
- What domains and motifs does it contain?
- How similar is it to homologous proteins?

ProteinInfo should NOT contain:

- genomic coordinates
- exon structure
- CRISPR guide RNAs
- donor templates
- editing strategies
- variant engineering logic

Those belong to later layers of the framework.

---

## Initial Fields

- protein_id
- gene_symbol
- species
- sequence
- domains
- motifs
- biological_role
- alignment statistics
    - best_match
    - percent_identity
    - coverage
    - evalue
- notes

---

## Design Principle

ProteinInfo should represent biological truth about a protein.

Everything related to manipulating DNA should live elsewhere.

This separation keeps the architecture modular, testable, and extensible.