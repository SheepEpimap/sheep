# Figure-code interpretation guide

Figure workflows use the layout:

```text
figures/<figure>/<panel>/<numbered script>
```

Numbered scripts encode execution order, not panel order. Read the panel README
before running any command.

## Relationship labels

- `PANEL_CODE`: directly calculates or plots content used in the named panel.
- `PANEL_INPUT`: generates a required intermediate used by panel code.
- `PARTIAL_PANEL`: generates only part of the final panel or additional outputs
  that must be selected during assembly.
- `VALIDATION`: checks required inputs or outputs but does not generate a panel.
- `ANALYSIS_SUPPORT`: supports the scientific analysis without reproducing the
  final panel layout.

The final Illustrator assembly, typography, annotations and panel labels are not
assumed to be reproducible unless a panel README explicitly says otherwise.

## Mapping policy

Assignments were checked against the final RGB PDF filename and visual content,
then against the later manuscript legends. The separate figure-caption document
was used where consistent. Conflicts are recorded in
[`CAPTION_RECONCILIATION.md`](CAPTION_RECONCILIATION.md).

Code was not assigned to a panel solely because it used a related biological
method. Reusable locus and regional-association plotters remain in
`figures/Shared_Workflows/` when exact regions, traits or command lines were not
supplied.

The detailed mapping is in [`../FIGURE_CODE_MAP.tsv`](../FIGURE_CODE_MAP.tsv),
and missing exact matches are recorded in
[`UNMATCHED_PANELS.tsv`](UNMATCHED_PANELS.tsv).
