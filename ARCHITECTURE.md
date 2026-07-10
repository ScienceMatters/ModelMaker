# ModelMaker Architecture

## Philosophy

ModelMaker is not a collection of scripts.

It is a framework for translating human genetic variants into experimentally actionable animal disease models.

The architecture is designed around biological entities rather than computational tasks.

The central biological object is the **Variant**.

Every computational module enriches the Variant with additional biological knowledge.

No module should duplicate work already performed by another module.

---

# Core Design

The ModelMaker core remains intentionally small.

It contains only:

- Variant
- Plugin
- Pipeline

Everything else is implemented as modules.

---

# Pipeline

Human Variant

↓

Orthology

↓

Protein Mapping

↓

Genome Mapping

↓

Annotation

↓

Editing

↓

Validation

↓

Report

Each module accepts a Variant and returns the same Variant enriched with new information.

---

# Guiding Principles

## Biology First

Recommendations should be biologically interpretable.

Every decision should include evidence and confidence.

---

## Reproducibility

Every report records:

- ModelMaker version
- UniProt release
- Ensembl release
- Database versions
- Alignment parameters

---

## Modularity

Every biological capability should exist as an independent module.

Examples:

- Orthology
- Alignment
- Genome Mapping
- Base Editing
- Prime Editing
- Annotation
- Structure
- Reporting

---

## Extensibility

New species should require adding a species module rather than modifying existing code.

New editing technologies should require adding a new editing module rather than changing the pipeline.

---

# Long-Term Vision

ModelMaker should answer one question:

"What is the best way to model this human genetic variant in an experimental organism?"

The answer should include:

- orthology
- conservation
- structure
- genome coordinates
- editing strategy
- validation strategy
- model suitability
- publication-ready figures

