# Mathemtical Modeling of NFkB Signaling Dynamics in Polarized Macrophages
Code and instructions to perform the parameter fitting for NFkB trajectories in 6 polarization states under Pam3CSK stimulation and subsequent downstream analyses as presented in "Stimulus-Response signaling dynamics characterize macrophage polarization states" (Singh et al, 2023).
## Things to know
1. Polarization state shorthands and order used:
     1. P3K = M0 (naive) cell with Pam3CSK ligand
     2. P3Kib = M:IFNβ cell with Pam3CSK ligand
     3. P3Kig = M:IFNγ cell with Pam3CSK ligand
     4. P3Ki0 = M:IL10 cell with Pam3CSK ligand
     5. P3Ki3 = M:IL13 cell with Pam3CSK ligand
     6. P3Ki4 = M:IL4 cell with Pam3CSK ligand
2. Biochemical parameters fit and order used:
     1. IκBα mRNA NFkB induced synthesis Km
     2. IκBα protein degradation rate
     3. TAK1 inactivation rate
     4. TLR2 synthesis rate
     5. ligand-receptor complex degradation rate
     6. IκBα transcriptional delay
     7. NFκB initial abundance
4. RMSD = root-mean-square deviation
5. ODE = ordinary differential equation
## Set-up
1. Install MATLAB on your computer (version R2020b or later)
2. Download the contents of this repository

## Parameter fitting (~1 week to run)
1. In MATLAB, navigate to the `param_fitting` folder that you downloaded from this repository
2. Run `parallel_runs_full` in the command line
     - NOTE: code set up to leverage parallelization with 50 cores. Change if needed.
   ### Output:
   6 folders, one for each polarization state, including:
      - line graphs with the trajectories of the top 10 parameter-fit simulations vs. experimental data (1 graph per cell)
      - `runtime.mat` with the total runtime per polarization state
      - `result_cell.mat` with 7 matrices:
        1. Row numbers of the cells used
        2. Randomly sampled initial parameters
        3. All sets of parameter fits per cell
        4. RMSDs for all parameter fits per cell
        5. End objective values for all parameter fits per cell (weighted RMSD + penalty)
        6. Sorted top 10 RMSDs per cell
        7. Sorted top 10 parameter fits per cell (corresponding to #6)
## ALTERNATE: Parameter fitting unit test (~7 min to run)
1. In MATLAB, navigate to the `param_fitting` folder that you downloaded from this repository
2. Run `unit_test` in the command line
     - NOTE: code set up to leverage parallelization with 4 cores. Change if needed.
   ### Output:
   Same as Part 1, except only 1 folder containing results from 2 initial conditions and 2 experimental cells

## Distance matrix calculation
- Code for calculating distance matrix from parameter distributions can be found in the `param_distance_calculation.m` file

## Data for downstream analyses found in `data` folder
- `params_3d.mat`: optimized parameter distributions in 3D format (1800 cells by 7 params by 10 best fits)
       - every 300 rows is a new polarization state. The order is as listed in [Things to know](#Things-to-know)
- `params_2d.csv`: the same optimized parameter distributions, but in 2D format (18000 rows by 7 params) (every 10 rows are the best fits for 1 cell)
       - every 3000 rows is a new polarization state. The order is as listed in [Things to know](#Things-to-know)
- `param_dist_matrix.csv`: distance matrix calculated based on average Jenson-Shannon distances for optimized param distributions between cells (1800 cells by 1800 cells)
       - every 300 rows is a new polarization state. the order is as listed in [Things to know](#Things-to-know)
- `feat_dist_matrix.csv`: distance matrix calculated based on euclidean distances for 10 experimental trajectory features between cells  (1800 cells by 1800 cells)
       - every 300 rows is a new polarization state. the order is as listed in [Things to know](#Things-to-know)

