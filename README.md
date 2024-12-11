<div align="center">
  <img src="man/figures/xcell2logo.png" width="30%">
</div>

# xCell 2.0: Robust Cell Type Enrichment Analysis

## Overview
**xCell 2.0** is an advanced R package and web application for cell type enrichment analysis, building upon the widely-used xCell algorithm. It introduces a training function that allows the utilization of custom reference datasets, making it adaptable to diverse tissue types and experimental conditions.

## Installation
The package is currently available on GitHub and is being prepared for Bioconductor submission:

```R
devtools::install_github('AlmogAngel/xCell2')
library(xCell2)
```

## Important Notes

<ul>
<li>xCell 2.0 produces enrichment scores, not percentages. It's designed for comparing across samples, not cell types.</li>
<li>xCell 2.0 requires variability among samples for accurate linear transformation.</li>
<li>The method is designed for mixed samples, not for determining cell of origin or single-cell analysis.</li>
</ul>

## Vignettes

[Introduction to Cell Type Enrichment Analysis with xCell 2.0](https://aran-lab.com/xcell2-vignette)


## Contributors

xCell 2.0 is developed by the Aran lab at the Technion - Israel Institute of Technology. Contact: Almog Angel (almog.angel at campus.technion.ac.il) Dvir Aran (dvir.aran at technion.ac.il)

## Citation

If you use xCell 2.0 in your research, please cite our paper:
[Angel A, Naom L, Nabel-Levy S, Aran D. xCell 2.0: Robust Algorithm for Cell Type Proportion Estimation Predicts Response to Immune Checkpoint Blockade. bioRxiv 2024.](https://doi.org/10.1101/2024.09.06.611424)

## License

GLP 3.0
