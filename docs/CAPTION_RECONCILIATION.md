# Figure and caption reconciliation

The supplied materials are not perfectly synchronized. To avoid attaching code
to the wrong figure number, this repository uses the following precedence:

1. final RGB PDF filename and visual content;
2. the later manuscript file (`Epimap manuscript(2026.6.13)-LF(7.14).docx`);
3. the figure-and-caption file (`Sheep Epimap Fig (2026.7.13).docx`) where it is
   consistent with the final PDFs.

The source documents and PDFs are not included in this code repository.

## Resolved numbering conflicts

| Final PDF identity | Final subject used in this repository | Conflicting source text | Resolution |
|---|---|---|---|
| Extended Data Fig. 2 | Super-enhancers and their potential functions | The figure-caption document labels the enhancer-link evaluation as Extended Data Fig. 2 | The RGB PDF named `Extended Data Fig.2.pdf` shows super-enhancer analyses, matching the later manuscript legend |
| Extended Data Fig. 3 | ATAC-guided selection and characterization of enhancer-target gene links | The figure-caption document labels super-enhancer analyses as Extended Data Fig. 3 | The RGB PDF named `Extended Data Fig.3.pdf` shows enhancer-link evaluation, matching the later manuscript legend |
| Supplementary Fig. 8 | Conservation, accessibility and DNA methylation profiles around chromatin states | Captions for Supplementary Figs. 8-11 are out of sequence in the figure-caption document | The RGB PDF identity is used; only panel 8d has a defensible supplied code match |
| Supplementary Fig. 9 | Summary characteristics of super-enhancers | The caption order places this content elsewhere | The RGB PDF identity is used |
| Supplementary Fig. 10 | Distribution of tissue-specific regulatory elements across tissues | The caption order places this content elsewhere | The RGB PDF identity is used |
| Supplementary Fig. 11 | Functional and developmental characterization of sheep regulatory elements | The caption order contains a duplicated Supplementary Fig. 11 label | The RGB PDF identity is used |

## Excluded outdated mapping

An older `Supplementary_Figure_08/08e` workflow plotted state-level methylation
after excluding E11 and values at or above 99.5%. It does not reproduce the
final Supplementary Fig. 8e profile plot and was therefore excluded. The updated
E1-E11 state-level methylation workflow is correctly retained under
`figures/Figure_02/02c/` for main Figure 2c.

This reconciliation is intentionally conservative. A related method or shared
plotting engine is not treated as an exact panel match without the required
panel-specific inputs and invocation.
