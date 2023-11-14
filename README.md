# MPM_MATLAB

This MATLAB code demonstrates the following variants of the material point method (MPM) for analyzing large deformation problems:
1. The B-spline material point method (BSMPM),
2. The Convected particle domain interpolation (CPDI) method, and
3. The B-spline convected particle domain interpolation (BSCPDI) method.

This code is associated with examples 5.1 & 5.2 of the following paper: 
A. Sadeghirad, "B-spline convected particle domain interpolation method", Engineering Analysis with Boundary Elements, 2023.
Example 5.1: Axis-aligned displacement in a square
Example 5.2: Expanding ring

SOURCE CODE:
This code includes one main source file, 'CPDI.m', and six auxiliary source files:
1. CPDI.m: main file to run the simulation, including the loop over time steps
2. Preprocess.m: prerocess tasks, including background grid generation, particle generation, define bounadry and initial conditions, ...
3. CalcGradSF.m: effective basis function calculation for the CPDI, BSMPM, and BSCPDI methods
4. Bspline.m: B-spline basis function calculation, which is called from CalcGradSF.m
5. Material.m: constitutive model (the neo-Hookean material model)
6. AnalyticalSolution.m: analytical solution of the problem
7. Snapshot.m: to plot simulation snapshots

HOW TO RUN:
Simply run CPDI.m in MATLAB.

OUTPUTS OF THE SIMULATION:
1. Four snapshots of the simulation will be plotted after each quarter of the total simulation time,
2. Displacement & energy error norms will be given upon the completion of the simulation
