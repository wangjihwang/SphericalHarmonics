#!/bin/csh

#SBATCH -p orion
#SBATCH --qos=batch
#SBATCH -A ome
#SBATCH -J generate_SpherHarmFunc_ne120
#SBATCH -o generate_SpherHarmFunc_ne120.out
#SBATCH -e generate_SpherHarmFunc_ne120.err
##SBATCH --mail-type=FAIL
##SBATCH --mail-user=Aaron.Wang@noaa.gov
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
##SBATCH --mem=32G

module load idl

idl << EOF 
.compile SpherHarm_recur2.pro
.run generate_SpherHarmFunc_ne120.pro
EOF

exit
