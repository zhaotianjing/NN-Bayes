#! /bin/bash

#this file to run simulation analysis
nRepy=20
for repy in $( eval echo {1..$nRepy} )
do
	cd y"$repy"
	sbatch /home/tianjing/paper_simu/simulation.sbatch $repy
	cd ..
done