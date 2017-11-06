#!/bin/csh
### Set the job name
#PBS -N changs

### Request email when job begins and ends
#PBS -m bea

### Specify email address to use for notification.
#PBS -M changs@email.arizona.edu

### Specify the PI group found with va command
#PBS -W group_list=ece677

### Set the queue to submit this job.
#PBS -q windfall
### Set the number of cpus up to a maximum of 128 #################
#PBS -l select=1:ncpus=4:mem=4gb

### Specify up to a maximum of 1600 hours total cpu time for the job
#PBS -l cput=864:0:0

### Specify up to a maximum of 240 hours walltime for the job
#PBS -l walltime=24:0:0

cd ~changs/ECE677/Project/svd/output

### Include this only if you need to convert data file from big to little endian
setenv F_UFMTENDIAN big  

source /usr/share/Modules/init/csh
module load intel-mpi

date
/usr/bin/time -o mpi_prog.timing mpirun -np 64 ~changs/ECE677/Project/svd/svd_mpi ~changs/ECE677/Project/svd/data/inputMatrices/inputMatrix2048.txt outLeftMatrix.txt outRightMatrix.txt outSingularValue.txt 100 reducedRightMatrix.txt simMatrix.txt markSV.txt
date
