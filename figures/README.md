# Figure code

Figure workflows use the panel-oriented layout:

```text
figures/<figure>/<panel>/<numbered script>
```

Run numbered scripts in order after reading the panel README and configuring
the required inputs and environment-specific paths. Final multipanel assembly,
typography and labels may be completed outside the computational scripts.

Reusable plotting and preprocessing utilities are stored in
`figures/Shared_Workflows/` or a figure-level `Shared_Workflows/` directory.

See [`../FIGURE_CODE_MAP.tsv`](../FIGURE_CODE_MAP.tsv) for the figure-panel
script index and [`../docs/FIGURE_CODE_NOTES.md`](../docs/FIGURE_CODE_NOTES.md)
for relationship labels.
