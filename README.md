# Landslide Mechanics: Stability on Curved Basal Surfaces

This repository contains the numerical modeling tools and analytical solutions developed to investigate the stability of deep-seated landslides on curved basal surfaces, with a specific focus on the role of glacial buttressing. 

The core of this work utilizes a depth-integrated limit equilibrium framework to explore how basal curvature and ice-rock density ratios govern the evolution of landslide topography and stability thresholds.

## Key Features
* **Analytical Solutions:** Implementations of the Taylor-expanded equilibrium thickness equations for ice-free landslides.
* **Numerical Models:** Solvers for the non-linear ODEs governing landslide geometry under glacial buttressing.
* **Stability Envelopes:** Tools to map the $(f(1-\\lambda), \\alpha)$ parameter space to identify regimes of \"self-stabilization.\"
* **Glacial Coupling:** Modules to simulate the transition from buttressed to unbuttressed states during glacier retreat.

## Citation
If you use this code or the associated models in your research, please cite the repository using the following DOI:

**DOI:** [10.5281/zenodo.19674416](https://doi.org/10.5281/zenodo.19674416)

**Reference:** Mallick, R., Finnegan, N., & Fielding, E. (2026). Deep-seated landslide stability on curved basal surfaces with and without glacier buttressing. Zenodo. https://doi.org/10.5281/zenodo.19674416

### Installation
```bash
git clone [https://github.com/mallickrishg/landslides.git](https://github.com/mallickrishg/landslides.git)
cd landslides
