# LinearCorrelator

Usage: linear-correlator [-h] <scriptname>

A basic linear correlator. Reads the following commands in <scriptname>.
Each command must be specified exactly once.
infile          Input data file
outfile         Output results file
dt              Time between each row of the input data file
max_time        Max time gap for correlations
columns         A string of characters, either 'x' (exclude the column),
                'c' (correlate the column normally), or 'f' (correlate the
                column after subtracting the mean) (e.g., xffffccc).
correlations    One pair of zero-indexed column indices for each correlation
                to calculate, notated as a,b for cross-correlating
                column a with later values of column b, e.g.,
                1,1 2,2 3,3 4,4 5,5 6,6 7,7 1,2 1,3 1,4 1,5 1,6 1,7 etc.

The data in the input file is assumed to have the following properties:
 - Whitespace-delimited
 - At least two columns; the first is time and the rest are the instantaneous values to be correlated
 - The time column is assumed to be evenly-spaced, such that a correlation between lines `i` and `i+k` is the same
as a correlation between lines `i+n` and `i+n+k`.

## Requirements

Requires an OpenMP >= 4.0.

## Examples
Given a data file `big_test.txt` with the following contents, where the columns correspond to
simulation time, hydrostatic pressure, and pressure tensor components 
$P_{xx}$, $P_{yy}$, $P_{zz}$, $P_{xy}$, $P_{xz}$, $P_{yz}$, respectively, 
```
# Fix print output for fix 3
0.000 1.5881 1.7441 1.4002 1.6201 0.025027 -0.230762 0.125430
0.005 1.5897 1.7491 1.3870 1.6330 0.020649 -0.243693 0.134716
0.010 1.5917 1.7495 1.3726 1.6529 0.015236 -0.255164 0.144528
...
49.990 1.3705 1.6128 1.1443 1.3544 0.133676 0.158138 -0.131165
49.995 1.3450 1.5989 1.0932 1.3428 0.121541 0.155157 -0.123050
50.000 1.3205 1.5862 1.0482 1.3271 0.106898 0.149418 -0.119044
```
we can compute the following autocorrelation functions using the respective scripts:

 1. The hydrostatic pressure fluctuations, up to and including a time gap of 2.0:
```
infile big_test.txt
outfile corr_hydrostatic.txt
dt 0.005
max_time 2.0
columns xfxxxxxx
correlations 1,1
```
 2. The shear pressure components, with 1001 data points in the result:
```
infile big_test.txt
outfile corr_shear.txt
dt 0.005
max_time 5.0
columns xxxxxccc
correlations 5,5 6,6 7,7
```
 3. All columns, up to a time gap of 10 (columns 1-4 should be
flagged with `-m` to remove equilibrium pressure contribution):
```
infile big_test.txt
outfile corr_all.txt
dt 0.005
max_time 10.0
columns xffffccc
correlations 1,1 2,2 3,3 4,4 5,5 6,6 7,7
```

## Benchmarking
In the `test/` directory, there is a large file named `big_test.txt` which contains data 
similar to that described above. Running the script `corr_script.txt`
should give the same result in `corr_test_check.txt` as is already written in `corr_test.txt`,
with or without the OpenMP acceleration, and independent of the optimization level.
Timing results of this benchmark will be examined in the near future and included in
the `test/` directory.

## To-Do
 -[x] Cross-correlators
 -[x] OMP multi-threading
 -[ ] Integrate long correlator and allow as optional command

