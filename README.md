# Chemical biology reveals the direct redeployment of functional target cell surface proteins to macrophages

This repository contains the final limma analyses used for the manuscript tables and figures.

Folders
- methods: statistical methods write-up for the limma pipeline
- scripts: R scripts used to generate the final results
- results: CSV outputs for the final comparisons only

Final comparisons included
- Surface enrichment Eaters vs Non-Eaters: 4A Light, 4B Heavy
- Whole-protein Eaters vs Non-Eaters: 3A Light, 3B Heavy, 3C Light+Heavy
- Co-culture exposed vs naive macrophages: EXMac vs Mac for 5A Light and 5B Heavy
- SW620 control: UE vs SW for 5A Light

Notes
- logFC is the model-estimated mean log2 difference between conditions.
- adj.P.Val is Benjamini–Hochberg FDR on limma p-values.
- Missing values are LOD-imputed and then log2-transformed as described in methods.
