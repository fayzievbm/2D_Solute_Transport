Concentration normalization and additional manuscript results
18 November 2026

The authors confirmed dimensionless concentration variables with one common
reference concentration. This clarification applies to all saved scenarios:
  c  = dimensional dissolved concentration per water volume / C_ref
  ca = dimensional active adsorbed amount per bulk volume / C_ref
  cp = dimensional passive adsorbed amount per bulk volume / C_ref
The same positive C_ref scales inlet concentration, capacities and stage
thresholds. It is not a separate capacity normalization for each phase.
No numerical value or physical solute identity is assigned to C_ref.
Writing m^3/m^3 would specifically impose a volume-fraction interpretation;
the manuscript uses the more general common-reference concentration ratios.

Only concentrations are normalized. Coordinates are in metres, time is in
seconds and the remaining physical parameters retain their stated units.
M = sum(area * (porosity*c + ca + cp)) therefore has units of m^2.
Cumulative input I, outlet O and degradation G have the same units.
C_ref*M is dimensional parent-solute amount per unit transverse thickness.
All inventory/input and degradation/input ratios are dimensionless.

No solver equations, parameter values, figures or saved simulation arrays were
changed by this clarification. Earlier explanatory files in startup_correction
describe the amount basis as unspecified. This note supplies the concentration
normalization and integrated-quantity units for the revised manuscript.

Table 4 is new post-processing of existing data, not an additional simulation.
The calculation sums only interior control volumes; algebraic boundary nodes
have no storage. It checks the phase sum against the summary and CSV ledger.
It also compares degraded/input fractions on the three existing meshes and
the finest time-step pair. The variable-minus-constant degraded fraction
changes sign between the intermediate and finest meshes. A small difference
in degradation performance is consequently not claimed.

To reproduce, use the pinned numerical environment supplied in the original
startup_correction/requirements.txt and run:
  python normalized_inventory_analysis.py PATH_TO_STARTUP_CORRECTION NEW_OUTPUT_DIR
The output directory must not exist. Inputs are read-only. The script produces
normalized_inventory_results.json and Table_4_percentages.csv.

The updated supplementary ZIP retains all original archive members unchanged
and adds this normalization note, post-processing script and results under
manuscript_normalization/. The main manuscript is supplied separately.
No repository deposit or DOI has been created; repository publication remains
an author action before submission.
