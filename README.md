# SphericalHarmonics
Using recurrence relations to generate scalar and vector Spherical Harmonics functions

#######################
# Jih-Wang Aaron Wang #
# 2023/5/29           #
#######################

  This program was created by Jih-Wang Aaron Wang. It is to generate scalar
and vector Spherical Harmonics Functions. Only two (Y and PSI) of the three
(Y, PHI, and PSI) functions are generated and enough to represent the system.
The recurrence relation is the key to generate the sequential functions, given
the functions at two previous wavenumbers. Note that there are several
recurrence relationship equations, but the ones that are presented here have
been tested extensively to be the most reliable and stable ones.

  The latitudes and longitudes of all the grid points, as well as the area of
each grid box, are needed beforehand. The *irregular* grid that is used here
is the Spectral Element grid of resolution ne120 (nearly 1/4 degree), which
is the standard grid for CESM, commonly-used NCAR GCM.

  The code was written in IDL, and was designed to run on single CPU on a single
node. The work started around the middle of 2016 and ended on 2021/8/14.
Please cite my name if you use the program.

###############

Important files:

SpherHarm_recur2.pro

generate_SpherHarmFunc_ne120.pro
generate_SpherHarmFunc_ne120.run
generate_SpherHarmFunc_ne120.csh

##################################

To run the IDL program on a computing node, simply do

sbatch generate_SpherHarmFunc_ne120.csh
