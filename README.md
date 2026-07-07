# ModelMaker / varianttool

A small Python package for mapping human disease variants onto orthologous zebrafish proteins, then scoring whether an exact or nearby conserved site is a plausible modeling target.

## Example files

- `examples/targets.csv` — human/zebrafish ortholog targets
- `examples/variants.csv` — mutation list
- `examples/cms_example.csv` — single-target CMS example

## Multi-target mapping

```bash
varianttool map --targets examples/targets.csv --variants examples/variants.csv --outdir out
```

This writes a combined `variant_map.csv` plus per-target alignment files in `out/<gene>/`.

## Single-target mapping

```bash
varianttool map --human P02708 --zfish "gene:chrna1 AND organism_id:7955" --mutations examples/cms_example.csv --outdir out
```

## Nearby functional-site logic

For residues that are not perfectly conserved, ModelMaker searches a local window around the target site for the nearest conserved aligned residue and reports it as a nearby functional-site candidate.
