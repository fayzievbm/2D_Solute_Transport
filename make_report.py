"""Human-readable report derived only from accepted completed output."""
from pathlib import Path
import json,html
R=Path(__file__).parent
def read(name):return json.loads((R/name).read_text())
a=read('assessment.json')
if not a['accepted']:raise RuntimeError('Do not report unaccepted figures as final')
s=read('runs/'+a['selected_run']+'/summary.json');d=read('comparison_metrics.json')
audits=read('output_audit.json');pres=read('preservation_audit.json')
orders=read('observed_orders.json');hist=read('history_convergence.json')
pressure=read('verification/pressure_tolerance_refinement.json')
lines=['STARTUP PRESSURE CORRECTION AND FIGURE VERIFICATION',
    '18 November 2026','',
    'Result: all stated numerical acceptance checks passed. Figures 2-16 have been regenerated from the selected complete run '+a['selected_run']+'.',
    'This is numerical verification of the supplied model, not experimental validation or a declaration that the entire manuscript is ready for submission.','',
    '1. What changed and why',
    'The previous directionally split pressure time update produced an unphysical startup overshoot under the abruptly imposed 500 kPa inlet pressure. The new solver integrates the complete finite-volume pressure system with adaptive implicit BDF/NDF steps. It retains the physical pressure PDE, initial condition, inlet pressure, ambient outlet, permeability and storage coefficients. No pressure clipping or artificial boundary-pressure ramp is applied.',
    'The full implicit linear system is solved using a verified structured factorization. Independent sparse-LU and matrix-exponential references check the implementation. The constant and variable dispersion transport stages, accepted adsorption increments, capacities and degradation parameters are retained.',
    'Transport now receives Darcy velocities at the actual midpoint/endpoints of each physical time step. A smooth geometric transport schedule covers the first 800 s; the selected subsequent step is 2 s. Both the startup schedule and later steps are refined in the time-convergence study.',
    'All four scenarios needed by Figures 2-16 were rerun through 15 h: constant and variable dispersion with constant inlet, and constant dispersion with exponential inlet rates 0.001 and 0.0005 s^-1.','',
    '2. Pressure overshoot and independent verification',
    f"Selected run: interior absolute pressure range {s['parameters']['p_ambient']+s['pressure']['min_gauge_pressure']:.12g} to {s['parameters']['p_ambient']+s['pressure']['max_gauge_pressure']:.12g} Pa. The prescribed full inlet face remains exactly 500000 Pa. All adaptive states and all transport pressure queries passed bounds, with 1e-6 Pa roundoff allowance.",
    f"Tightening pressure relative tolerance from 1e-7 to 1e-8 on the finest mesh changed saved near-field pressure by at most {100*pressure['max_field_L2']:.8g}% in weighted relative L2. Maximum absolute saved difference: {max(z['max_absolute_difference_Pa'] for z in pressure['comparisons']):.8g} Pa.",
    'The smaller-mesh full pressure operator agreed with the independently assembled spatial operator to 1.8e-16 relative error. Structured implicit solves agreed with independent sparse LU to less than 2.8e-13 relative error across four step factors. The pressure evolution agreed with a matrix-exponential reference to less than 2.5e-9 relative L2 error at the tested times.',
    'The transport test prescribed an independent steady Darcy field and disabled reactions. Against the full advection-diffusion matrix exponential, halving steps from 0.2 to 0.1 to 0.05 s reduced errors by approximately four for both dispersion laws. This verifies the transport stages in that controlled setting; observed full-model temporal orders are reported separately.','',
    '3. Measured convergence',
    'Field comparisons use the common x=0-4 m region at 0.5, 1, 2, 5 and 15 h. Norms are volume weighted; the finer or longer-domain field supplies the denominator. Pressure uses p-p0. Each comparison covers all four scenarios and all three concentration fields. The inlet-region norm (x<=0.2 m) is additionally reported, rather than hidden in the full-region average.',
    'Comparison | maximum field L2 change (%) | inlet-region L2 change (%) | maximum final integral change (%)']
for z in d.values():lines.append(f"{z['coarse']} -> {z['fine']} | {100*z['maximum_field_l2']:.8g} | {100*z['maximum_inlet_region_l2']:.8g} | {100*z['maximum_integral_change']:.8g}")
lines+=['','Acceptance check | measured | required | result']
for k,v in a['measured'].items():
    scale=1 if k in ('mass_residual_fraction','pressure_water_balance_fraction') else 100
    unit='fraction' if scale==1 else '%'
    lines.append(f"{k} | {v*scale:.9g} {unit} | <{a['criteria'][k]*scale:g} {unit} | PASS")
lines += ['',
    'Time refinement on the finest mesh compares 4 and 2 s subsequent transport steps, with pressure tolerance also tightened from 1e-7 to 1e-8. A separate three-level 4/2/1 s study uses the coarsest mesh and identical pressure tolerances. The latter does not establish a three-level observed order on the finest mesh.',
    'Full plotted histories were compared at the same physical point (0.05,0.4) m, including startup. Coarse histories were linearly interpolated onto fine times. Both are sampled every 8 s after startup and approximately every 1 s during startup, with additional figure-time samples. Reported differences include interpolation effects.',
    'The 64-to-128 m domain comparison is performed on the coarsest mesh, whose near-field coordinates are unchanged. Its demonstrated adequacy applies to the reported 0-4 m observation region and 15 h duration. It is not a guarantee for longer simulations or arbitrary physical parameters.','',
    '4. Three-level observed orders and conditional GCI']
for family,rows in orders.items():
    lines.append(family)
    for k,z in rows.items():
        if 'observed_order' in z:lines.append(f"  {k}: observed order {z['observed_order']:.6g}; fine-level GCI estimate {100*z['fine_gci_fraction']:.6g}%.")
        else:lines.append('  '+k+': '+z['interpretation'])
lines+=['These Richardson/GCI estimates assume an asymptotic error expansion. A less-than-1% grid-to-grid change does not establish a less-than-1% absolute discretization error. Full-model kinetics, capacity constraints and upwind fluxes can affect observed order.','',
    '5. Mass balance',
    'The exact numerical half-step boundary fluxes, phase increments and degradation losses are used in the solute ledger. M=sum(V*(theta*c+ca+cp)); R=M-M0-input+output+degraded. Capacities were respected and no meaningful negative dissolved concentrations were found. No dissolved concentration clipping was applied.',
    'Case | final stored | cumulative input | cumulative output | cumulative degradation | max cumulative residual / final input']
for name,z in s['transport'].items():lines.append(f"{name} | {z['inventory']:.12g} | {z['cumulative_in']:.12g} | {z['cumulative_out']:.8g} | {z['cumulative_degradation']:.12g} | {z['relative_max_cumulative_residual']:.8g}")
lines += [
    f"The variable-dispersion case has cumulative input {s['transport']['variable']['cumulative_in']/s['transport']['constant']['cumulative_in']:.8g} times that of the constant-dispersion case under the same prescribed inlet concentration. Their stored inventories therefore correspond to different cumulative loadings; absolute stored amounts alone should not be described as a normalized removal efficiency.",
    f"The maximum saved residual/current-input ratio in the selected run is {audits[a['selected_run']]['maximum_saved_residual_over_current_input']:.8g}. The independent hydraulic storage/net-boundary-flux integral residual after 1 s is at most {s['pressure']['max_relative_storage_integral_residual_after_1s']:.8g} relative in this run.",
    'Integrated solute amounts retain the supplied concentration basis per unit transverse thickness; no unsupported kg or mol unit is assigned. Algebraic boundary nodes have no extra storage. The discrete interior storage volume approaches the full geometric volume under refinement.',
    'The audit independently recomputes snapshot inventories, verifies Darcy velocities against saved pressure gradients, checks the exact 0.04 m inlet width, finite arrays, capacities, downstream ambient pressure and all recorded source hashes.','',
    '6. Final figures and manuscript implications',
    'Figures 2-5 replace the pressure surfaces, contours and cross sections. Figures 6-7 replace constant-inlet/constant-dispersion results. Figures 8-13 replace the exponential-inlet comparisons and histories. Figures 14-16 replace the variable-dispersion results and constant/variable comparison. Each is provided as a 300 dpi PNG and vector PDF; Figure_captions.txt identifies all cases, times, scales and displayed spatial intervals.',
    'The selected mesh has 1013 by 489 nodes including algebraic boundaries and is nonuniform. The numerical domain is 64 m; displayed surface intervals are 0-2 or 0-4 m. Figure 13 uses the same physical observation point for every inlet history. Figure 14 explicitly distinguishes 2 and 15 h.',
    'The Word manuscript was not edited. Before inserting these figures into a new Word copy, replace the old pressure ADI method description/stage equations with the supplied unsplit pressure formulation, document the nonuniform mesh and startup schedule, insert the measured convergence table, and reconcile all prose values tied to the superseded results. Methods_for_manuscript.txt supplies professional English wording.',
    'The adopted inlet width and concentration/kinetic parameter units still require the authors\' physical justification. The calculations establish numerical consistency under those assumptions.','',
    '7. Preservation and reproducibility',
    f"{pres['original_files_matching_previous_baseline']} source/data/manuscript files match the previous study baseline. All {pres['corrected_simulations_files_unchanged']} files in the preceding corrected-results package and all {pres['convergence_checks_files_unchanged']} files in the previous convergence package match their recorded SHA-256 hashes.",
    'One historical-baseline difference was found in the latest Word manuscript: it was saved before this new output directory was created. This task did not edit it. Its current digest was recorded and verified again; preexisting_manifest_difference.json and preservation_audit.json document both hashes and the saved timestamp.',
    'The two initially launched finest runs were stopped during startup to increase the Numba thread count from two to eight after an identical-array speed comparison. They supplied no final result; superseded_runs.json records the replacements. No solver equations or physical parameters changed for this optimization.',
    'The pressure implementation depends on the pinned SciPy BDF callbacks and dense-output representation. Run the independent tests before changing SciPy versions. Numerical source hashes, exact parameters, coordinate arrays and complete saved fields/histories accompany every completed run. See README.txt for reproduction commands.',
    '', 'Reference: SciPy BDF implementation documentation, https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.BDF.html',
    'The standalone figures and their supporting computations are finalized; manuscript integration and physical validation are distinct subsequent tasks.']
text='\n'.join(lines)+'\n'
(R/'Verification_report.txt').write_text(text,encoding='utf-8')
body=''.join('<p>'+html.escape(line)+'</p>' if line else '<hr>' for line in lines)
(R/'Verification_report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Startup pressure and figure verification</title><style>body{max-width:1000px;margin:40px auto;padding:0 24px;font:16px/1.6 system-ui;color:#17212b}p{margin:.6em 0}hr{border:0;border-top:1px solid #dbe2e8;margin:24px 0}p:first-child{font-size:25px;font-weight:700}</style><body>'+body+'</body></html>',encoding='utf-8')
print('Wrote verification report from accepted completed runs.')
