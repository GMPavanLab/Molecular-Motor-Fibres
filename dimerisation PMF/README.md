
# Dimerisation PMF

Compute and plot dimerisation PMF for all combinations of motor cores, both for AA and CG.

There is one main runfile: 

```run_all.sh``` Runs a 500 metadynamics simulation on a selected pair of motors. Also reweights the bias.

```analysis.ipynb``` is used to extract and plot the potentials against each other.

NOTE that there are two types of CG motors, the standard and v2. To compute PMF on v2, one needs to change the imports of itps in structure_files/CG_*.top 