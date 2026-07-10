# Variant Engine

Status: Draft

Author: Matt McCarroll

Version: 0.1

---

## Motivation

Every capability in ModelMaker ultimately operates on a single biological concept:

A human genetic variant.

Current bioinformatics software often passes dictionaries, JSON blobs, or disconnected result objects between analysis steps.

As ModelMaker expands to include orthology, structural biology, genome mapping, editing design, and phenotyping, this approach becomes difficult to maintain.

The Variant Engine introduces a single canonical biological object that is progressively enriched throughout the pipeline.

---

## Design Goals

The Variant should represent biology.

It should never contain algorithmic logic.

Algorithms belong to modules.

The Variant stores biological knowledge.

---

## Guiding Principle

The Variant owns the biology.

Modules own the algorithms.

The Pipeline owns execution.

Reports own presentation.

---

## Lifecycle

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

Each stage enriches the same Variant object.

Nothing is discarded.

---

## Core Objects

Variant

VariantInput

OrthologyInfo

ProteinInfo

GenomeInfo

AnnotationInfo

EditingInfo

ValidationInfo

ReportInfo

Metadata

---

## Plugins

Every module receives

Variant

↓

returns

Variant

Modules should never directly communicate.

The Variant is the only shared object.

---

## Future

Eventually Variant should become serializable.

A Variant file should completely describe the computational state of an analysis.

Example

CHRNE_P121L.mm

could be reopened years later with all annotations preserved.


## Non-responsibilities

Variant never performs alignments.

Variant never downloads data.

Variant never calls APIs.

Variant never designs guides.

Variant never writes reports.

Variant never performs calculations.

Variant only stores biological state.