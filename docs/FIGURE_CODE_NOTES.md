# Figure workflows

Figure workflows use the layout:

```text
figures/<figure>/<panel>/<numbered script>
```

Numbered scripts encode execution order. Each panel README describes the input
files, expected outputs and any environment-specific settings.

## Relationship labels

- `PANEL_CODE`: calculates or plots content used in the named panel.
- `PANEL_INPUT`: generates an intermediate used by panel code.
- `PARTIAL_PANEL`: generates a component of a multipanel output.
- `VALIDATION`: checks required inputs or generated outputs.
- `ANALYSIS_SUPPORT`: supplies an upstream or shared analysis step.

Use [`../FIGURE_CODE_MAP.tsv`](../FIGURE_CODE_MAP.tsv) to locate scripts by
figure and panel. Shared locus, enhancer and association plotting utilities are
stored under `figures/Shared_Workflows/`.
