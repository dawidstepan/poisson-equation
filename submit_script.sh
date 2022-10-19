#!/bin/bash
# Number of cores
#SBATCH -c 20
# Runtime of this jobs is less then 20 minutes
#            (hh:mm:ss)
#SBATCH --time=00:20:00
# Clear the environment
module purge > /dev/null 2>&1

rm results.txt

for s in 1 2 4 8 10 20
do
	for res in 256 512 1024 2048
	do
		export OMP_NUM_THREADS=$s
		export OMP_SCHEDULE=static
		./stenciljacobi $res $s static

		export OMP_NUM_THREADS=$s
		export OMP_SCHEDULE=static,1
		./stenciljacobi $res $s static1
		
		export OMP_NUM_THREADS=$s
		export OMP_SCHEDULE=dynamic
		./stenciljacobi $res $s dynamic
	done
done
