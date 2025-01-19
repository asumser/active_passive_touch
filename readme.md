Code to run analyses and generate figures for Sumser et al. "Differential representation of active and passive touch in mouse somatosensory thalamus". 
Tested on MATLAB 2021b and requires MATLAB's Parallel Computing and Signal Processing Toolboxes, as well as  "CircStat for Matlab" to run (https://de.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics) [1]

(1) run "whisker_decomposition.m" with the respective path to downloaded data. The script extracts amplitude, phase, velocity parameters from whisker angle data and saves to file
(2) run AnalysisFigs.m which extracts the relevant info and plots figures with the same style as in the paper


References:
[1] P. Berens, CircStat: A Matlab Toolbox for Circular Statistics, Journal of Statistical Software, Volume 31, Issue 10, 2009