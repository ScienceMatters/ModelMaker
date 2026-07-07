# ModelMaker

## Vision

ModelMaker is an open-source platform for translating human genetic variants into experimentally actionable animal disease models.

Our goal is to reduce the time between identifying a clinically relevant human mutation and generating a genetically faithful experimental model by integrating comparative genomics, structural biology, genome editing, and biological validation into a single reproducible workflow.

Rather than functioning as a collection of alignment or CRISPR design scripts, ModelMaker serves as an intelligent decision-support system for disease model development.

Ultimately, ModelMaker should answer one question:

> "What is the best way to model this human genetic variant in my experimental organism?"

---

# Philosophy

Building an animal model involves much more than identifying an ortholog.

Researchers must determine

- whether the affected residue is evolutionarily conserved,
- whether nearby residues may provide a better functional model,
- whether the mutation is technically editable,
- which genome editing technology is most appropriate,
- how the mutation should be validated,
- and ultimately whether the resulting model is likely to reproduce human biology.

ModelMaker aims to integrate these questions into a single transparent workflow.

---

# Guiding Principles

## Biology First

Every recommendation should be biologically interpretable.

ModelMaker should explain **why** it recommends an edit rather than acting as a black box.

---

## Reproducibility

Every output should be reproducible.

Each report should document

- sequence versions
- transcript versions
- database sources
- alignment parameters
- software version

allowing results to be regenerated years later.

---

## Extensibility

Every major capability should exist as an independent module.

Examples include

- Orthology
- Protein Alignment
- Functional Conservation
- Structural Annotation
- Genome Mapping
- Base Editing
- Prime Editing
- HDR Design
- Validation
- Reporting

Future species and editing technologies should require minimal changes to the existing codebase.

---

## Community

ModelMaker is intended to be useful beyond zebrafish.

Although zebrafish is the primary development organism, the architecture should support additional vertebrate and invertebrate model systems.

---

# Long-Term Goal

ModelMaker should become the standard workflow for computational disease-model design.

Researchers should be able to begin with a human genetic variant and receive a complete experimental blueprint including

- ortholog identification
- residue conservation
- structural interpretation
- genome editing strategy
- validation plan
- model confidence score
- publication-quality figures

within a single reproducible pipeline.

---

# Current Development Roadmap

## Core

✓ Ortholog identification

✓ Protein alignment

✓ Residue mapping

✓ Mature ↔ canonical numbering

✓ Nearby functional residue identification

---

## Next

□ Protein domain annotation

□ Model suitability scoring

□ Base editing design

□ Prime editing design

□ Guide RNA ranking

□ Validation primer design

---

## Future

□ AlphaFold/PDB integration

□ Multi-species support

□ ClinVar integration

□ OMIM integration

□ ZFIN integration

□ Behavioral phenotype integration

□ Brain imaging integration

□ Drug-response integration

□ Interactive HTML reports

---

# Mission

Accelerate translational biology by making precise disease-model design accessible, reproducible, and biologically informed.

# Workflow

Human Variant
        │
        ▼
 Comparative Genomics
        │
        ▼
 Functional Conservation
        │
        ▼
 Structural Biology
        │
        ▼
 Genome Editing
        │
        ▼
 Experimental Validation
        │
        ▼
 Disease Model
        │
        ▼
 Phenotype
        │
        ▼
 Therapeutic Discovery


 ## Non-goals

ModelMaker is not intended to replace genome browsers, alignment software, or CRISPR design tools.

Instead, ModelMaker integrates existing resources into a coherent, biologically informed workflow that helps researchers make better decisions when designing disease models.
