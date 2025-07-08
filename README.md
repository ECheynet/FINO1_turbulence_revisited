# FINO1 Turbulence Revisited

This repository revisits wind turbulence statistics from the FINO1 offshore research platform, based on sonic anemometer data collected in 2007 and 2008. The original post-processed time series are openly available on Zenodo:

- [Zenodo Record 15826899](https://zenodo.org/records/15826899)
- [Zenodo Record 15826678](https://zenodo.org/records/15826678)

The focus is on reanalysing turbulence characteristics using updated methods, comparing mean flow and friction velocity estimates from ERA5 with in situ measurements, and fitting wind spectra.

---

## Contents

The repository includes:

- **`dataAnalysed/`**  
  Contains all processed data:
  - `data1.nc`, `data2.nc`: ERA5 reanalysis data for 2007 and/or 2008.
  - 24 `.mat` files: Monthly flow statistics derived from the FINO1 time series.

- **`Documentation1.mlx`**  
  A MATLAB Live Script that demonstrates how to reproduce comparisons of mean wind speed and friction velocity between ERA5 and sonic anemometer data.

- **`Documentation2.mlx`**  
  A MATLAB Live Script focused on fitting turbulence spectra and extracting key parameters for monthly analysis.

- **`fit_and_plot_spectra.m`**  
  A supporting MATLAB function used in `Documentation2.mlx` for fitting and visualising turbulence spectra.

---

## Relevance and Applications

This dataset and analysis are relevant for:

- Validating turbulence models such as Kaimal or Mann for wind load calculations.
- Analysing single-point auto-spectral and cross-spectral densities of wind turbulence.
- Evaluating wind turbine design standards, such as IEC 61400-1 and IEC 61400-3.

---

## References

The data have been used in the following peer-reviewed studies:

- Cheynet, E., Jakobsen, J. B., & Obhrai, C. (2017). *Spectral characteristics of surface-layer turbulence in the North Sea*. Energy Procedia, 137, 414–427. Elsevier.
- **Cheynet, E., Jakobsen, J. B., & Reuder, J. (2018). *Velocity spectra and coherence estimates in the marine atmospheric boundary layer*. Boundary-Layer Meteorology, 169(3), 429–460. Springer Netherlands.**  
  → This is the main paper in which the original data analysis was conducted.
- Cheynet, E. (2019). *Influence of the Measurement Height on the Vertical Coherence of Natural Wind*. Lecture Notes in Civil Engineering, 27, 207–221. Springer.

This repository provides an updated look at the FINO1 dataset and complements the work published in the 2018 study. This is the first version of this repository; Some bugs may still be present


---

