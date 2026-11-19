# Two-dimensional reactive solute transport

Code, numerical results, and figures accompanying the manuscript **Dispersion and Inlet History Effects on Two-Dimensional Reactive Solute Transport under Transient Darcy Flow**.

This repository contains the complete verified numerical supplement: the corrected solver, all seven completed production/convergence runs, independent verification results, Figures 2–16, and the additional normalized-inventory analysis for Table 4. The 171 files from `Manuscript_Numerical_Supplement_Updated.zip` are preserved byte for byte. This deposit adds navigation and an independent file-integrity check; it does not change the model or rerun the simulations.

## Start here

- [Verification report](Verification_report.txt): numerical methods, convergence, conservation, and limitations.
- [Study protocol](Study_protocol.txt) and [methods](Methods_for_manuscript.txt): model assumptions and verification criteria.
- [Convergence table](Final_convergence_table.csv) and [acceptance assessment](assessment.json).
- [Figures](figures/) in PNG and vector PDF, with [captions](figures/Figure_captions.txt). Download and open `figures/Gallery.html` locally to browse the gallery.
- [Concentration normalization](manuscript_normalization/NORMALIZATION_README.txt) and [Table 4 percentages](manuscript_normalization/Table_4_percentages.csv).
- [Original technical README](README.txt): solver components and output formats.

## Model and scenario design

Transient pressure determines longitudinal and transverse Darcy velocities. Conservative transport couples dispersion with adsorption in active and passive zones and first-order degradation. The inlet opening is fixed at 0.04 m; the downstream pressure is ambient pressure. The selected domain is a numerical truncation of the semi-infinite model.

Each completed run contains these four scenarios:

| Output prefix | Dispersion | Inlet concentration |
| --- | --- | --- |
| `constant` | Constant anisotropic | Constant |
| `variable` | Spatially variable scalar | Constant |
| `exp_0.001` | Constant anisotropic | Exponentially increasing, rate 0.001 s^-1 |
| `exp_0.0005` | Constant anisotropic | Exponentially increasing, rate 0.0005 s^-1 |

The dispersion comparison changes anisotropy and coefficient magnitudes together with spatial variation. It does not isolate heterogeneity alone. Adsorption and degradation parameters are fixed across these comparisons.

All concentration variables use **one common reference concentration**. They are dimensionless; coordinates are in metres and time is in seconds. Area-integrated normalized inventory, input, outlet loss, and degradation have units of m². Their ratios are dimensionless. The normalization note supersedes earlier statements that the concentration basis was unspecified.

## Saved runs

All runs reach 15 h and include flow fields, concentration arrays, histories, and a parameter/source-hash summary.

| Directory under `runs/` | Purpose |
| --- | --- |
| `mesh0_dt2` | Coarsest spatial level |
| `mesh1_dt2` | Intermediate spatial level |
| **`mesh2_dt2_t8`** | **Selected finest run for Figures 2–16 and Table 4** |
| `mesh2_dt4_t8` | Finest-mesh temporal comparison |
| `mesh0_dt4` | Coarse-mesh temporal comparison |
| `mesh0_dt1` | Coarse-mesh temporal comparison |
| `domain128_dt2` | Doubled downstream domain |

`dt` describes the transport step after the refined startup schedule, not the first step. Pressure uses its own adaptive implicit integration. `_t8` records eight Numba threads; the equations are unchanged. Array time keys are described in [README.txt](README.txt); use the saved physical times rather than interpreting the key as an actual startup-step count.

## Installation and reproduction

Use a Python environment compatible with the exact versions in [requirements.txt](requirements.txt). The archived results were produced on Windows with working Numba/LLVM and writable cache directories. Recheck the independent tests if changing dependencies; the pressure implementation uses SciPy BDF internal interfaces.

```bash
git clone https://github.com/fayzievbm/2D_Solute_Transport.git
cd 2D_Solute_Transport
python -m venv .venv
```

Activate the environment (`.venv\Scripts\Activate.ps1` in PowerShell, or `source .venv/bin/activate` on Unix-like systems), then:

```bash
python -m pip install -r requirements.txt
python tools/check_integrity.py
python manuscript_normalization/normalized_inventory_analysis.py . local_reproductions/table4
```

The last command reproduces Table 4 from existing arrays without running a simulation. The output directory must be new. Compare its CSV and JSON with the archived files in `manuscript_normalization/`.

To rerun all four scenarios on the selected mesh:

```bash
python run.py --name repeat_finest --level 2 --L 64 --dt 2 --pressure-rtol 1e-8
```

The solver refuses an existing run name. It defaults to two Numba threads; the archived finest runs used `NUMBA_NUM_THREADS=8` and `OPENBLAS_NUM_THREADS=2`. A full finest-mesh simulation is computationally substantial.

Independent pressure and transport checks are available as `verify_pressure.py` and `verify_transport.py`. These write into `verification/`; run them in a separate working copy to retain the archived verification files unchanged. Likewise, analysis/report scripts write output files, and `make_verified_figures.py` requires an accepted assessment and refuses an existing `figures/` directory.

### Historical helper scripts

The solver and normalized-inventory analysis use this repository's data. Some original workflow helpers also preserve earlier local projects: `audit_results.py` calls both `audit()` and `preservation()`, and the latter requires the author's earlier workspace. For a numerical-data audit in a separate working copy, use `python -c "import audit_results; audit_results.audit()"`. The archived preservation reports document the original checks; `tools/check_integrity.py` verifies the portable repository files without those external paths. The older README and normalization note contain historical statements about a repository deposit not yet existing; this README records the current GitHub deposit.

## Interpretation and status

The selected results pass the stated numerical acceptance checks. Both dispersion cases have approximately 67% stored parent solute and 33% cumulative degradation at the final time. Their small difference in degraded fraction changes sign with mesh refinement and is not evidence of improved degradation performance. Increased retained inventory must be interpreted alongside the increased cumulative input.

These are idealized numerical results, not independent experimental validation or a demonstrated remediation benefit. Application to a specific aquifer requires a dimensional concentration scale, justified parameters, calibration, and observational validation. The accompanying manuscript remains in preparation.

For traceability, [deposit_manifest.json](deposit_manifest.json) records the source archive hash and each deposited scientific file's SHA-256. Cite this repository with the commit identifier used in your work. No manuscript DOI or repository DOI has been assigned here.
