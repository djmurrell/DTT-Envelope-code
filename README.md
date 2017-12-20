# DTT-Envelope-code

This code is used to produce the figures in the revised manuscript "A global envelope test to detect non-random bursts of trait evolution" (an earlier version of which can be found here doi: https://doi.org/10.1101/175968).

The file 
  figures.true-false.positives.R 
 is R script to produce figures 2 and 3 that investigate the false and true positive rates of the candidate null model tests.

The file 
  empirical.tests.R 
is used to produce figure 4 and is the recommended pipeline for using datasets that would be analysed via the function dtt in the R library geiger.

The file
  multi-trait.rank.test.power.R
is used to produce the data for the figure 4 that investigates the power to detect non-random trait evolution when there are two traits using the rank envelope test on the multivariate disparity value.
