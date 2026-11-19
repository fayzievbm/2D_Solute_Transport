Startup pressure correction, convergence checks and replacement figures

All files in this directory are NEW. Original scripts, data, manuscripts,
the previous corrected package and the previous convergence review are
preserved. See preservation_audit.json for SHA-256 verification.

Main entry point: run.py. For example, from this directory:
  python run.py --name repeat_finest --level 2 --L 64 --dt 2 --pressure-rtol 1e-8
The requested output name must be new. Existing runs are never overwritten.
All four scenarios are run by default. requirements.txt pins the runtime.
OPENBLAS_NUM_THREADS=2 was used. NUMBA_NUM_THREADS=8 was used for the two
finest production runs (suffix _t8), and 2 for the other runs. A short
thread-count comparison produced identical concentration arrays; the
larger count reduced runtime. See verification/thread_benchmark.json.
The Windows
runtime requires working Numba/LLVM and writable cache directories. run.py
creates a separate short cache under ../scn for each process. Caches and
installed dependencies are not part of the scientific output package.

Solver files:
- corrected_solver.py: physical clock, startup schedule, coupling, outputs.
- pressure_implicit.py: complete pressure operator and adaptive implicit
  solve. The pressure equation is no longer advanced using ADI.
- blocked_solver.py: retained mesh/state/parameter/field utilities. Its
  legacy step() and simulate() functions are NOT used by run.py.
- blocked_transport.py: retained conservative ADI transport and ledger.
  Its legacy advance() pressure update is NOT used by run.py.
- blocked_base.py, reference_solver.py: retained flux/reaction utilities
  and physical parameters. Their legacy entry points are not production.
- floating_point.py, worker_floating_point.py: documented Windows handling
  of subnormal double values only. No fastmath or meaningful-value cutoff.

Production and verification runs:
  mesh0_dt2, mesh1_dt2, mesh2_dt2_t8: three nested meshes; 2 s subsequent step.
  mesh2_dt4_t8: finest mesh, 4 s subsequent step, pressure rtol 1e-7.
  mesh0_dt4, mesh0_dt1: additional temporal levels, pressure rtol 1e-8.
  domain128_dt2: doubled domain on level 0, all four scenarios.
All other runs use pressure rtol 1e-8. All pressure atol values are 1e-7 Pa.
The pressure solver chooses its own adaptive steps. The transport step
listed in a run name applies after the smooth first 800 s, not from t=0.

Run output format:
- summary.json: all parameters, mesh sizes, actual step count, startup
  schedule, pressure and solute balance statistics, runtime, source hashes.
- flow_snapshots.npz: coordinates, volumes, inlet fractions, pressure in Pa
  and Darcy face velocities. Prefix keys are physical time / subsequent
  transport dt, NOT the actual startup step count. *_t stores time in s.
- <case>_snapshots.npz: complete c, ca and cp arrays at saved physical times.
  Full arrays include algebraic boundary nodes; storage uses interiors.
- <case>_history.csv: time, inventories, boundary integrals, degradation,
  residual and physical-point values. Headers identify all columns.
No concentration mass unit has been invented. Read Methods_for_manuscript.txt.

Independent tests:
  python verify_pressure.py
  python verify_transport.py
  python pressure_benchmark.py --level 2 --rtol 1e-7
  python pressure_benchmark.py --level 2 --rtol 1e-8
  python compare_pressure.py
Their results are in verification/. Pressure-only NPZ q fields are gauge
pressure and contain INTERIOR nodes; this differs from flow_snapshots.npz.

Analysis sequence after all runs complete:
  python analyze.py
  python convergence_order.py
  python assess_histories.py
  python audit_results.py
  python assess.py
Only after assessment.json reports accepted=true:
  python make_verified_figures.py
  python make_report.py
The figure generator refuses a failed assessment and an existing figures
directory. Scientific figures are supplied as 300 dpi PNG and vector PDF,
with separate captions. The Word document is not modified by these scripts.

Review Study_protocol.txt, Verification_report.txt, Final_convergence_table.csv,
Methods_for_manuscript.txt and figures/Figure_captions.txt first. Full detailed
metrics and machine-readable audit results accompany the concise report.
Open figures/Gallery.html to browse all numbered figures and their captions.
The figure-only ZIP includes this gallery and PNG/PDF exports; the full ZIP
also includes all completed simulation arrays, histories, source and audits.

Known reproducibility dependency: the structured pressure solve replaces
SciPy BDF linear-solver callbacks and uses its dense-output polynomial
representation. These interfaces were checked in the pinned SciPy version.
Repeat the independent tests when changing SciPy versions.
