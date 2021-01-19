#!/bin/bash
#SBATCH -J S16SCD
#SBATCH -n 1
#SBATCH -t 240:00:00
#SBATCH --mem=500m

./dcd_traj_shape ../dcd/cM0.0000.dcd 350.0 1037 cM0.0000.shape
./dcd_traj_shape ../dcd/cM0.0001.dcd 350.0 1037 cM0.0001.shape
./dcd_traj_shape ../dcd/cM0.0002.dcd 350.0 1037 cM0.0002.shape
./dcd_traj_shape ../dcd/cM0.0004.dcd 350.0 1037 cM0.0004.shape
./dcd_traj_shape ../dcd/cM0.0006.dcd 350.0 1037 cM0.0006.shape
./dcd_traj_shape ../dcd/cM0.0008.dcd 350.0 1037 cM0.0008.shape
./dcd_traj_shape ../dcd/cM0.0010.dcd 350.0 1037 cM0.0010.shape
./dcd_traj_shape ../dcd/cM0.0012.dcd 350.0 1037 cM0.0012.shape
./dcd_traj_shape ../dcd/cM0.0015.dcd 350.0 1037 cM0.0015.shape
./dcd_traj_shape ../dcd/cM0.0020.dcd 350.0 1037 cM0.0020.shape
./dcd_traj_shape ../dcd/cM0.0025.dcd 350.0 1037 cM0.0025.shape
./dcd_traj_shape ../dcd/cM0.0030.dcd 350.0 1037 cM0.0030.shape
./dcd_traj_shape ../dcd/cM0.0040.dcd 350.0 1037 cM0.0040.shape
./dcd_traj_shape ../dcd/cM0.0050.dcd 350.0 1037 cM0.0050.shape
./dcd_traj_shape ../dcd/cM0.0100.dcd 350.0 1037 cM0.0100.shape
./dcd_traj_shape ../dcd/cM0.0200.dcd 350.0 1037 cM0.0200.shape
./dcd_traj_shape ../dcd/cM0.0300.dcd 350.0 1037 cM0.0300.shape
