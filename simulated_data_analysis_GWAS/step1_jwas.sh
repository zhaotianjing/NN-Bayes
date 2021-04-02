#! /bin/bash

#this file is to run JWAS
nRepy=20
methodAll=("hmc_bayescpi_2" "bayescpi")

for repy in $( eval echo {1..$nRepy} )
do  
	cd y$repy
	for method in "${methodAll[@]}"
	do
		cd run_"$method"
		sbatch /group/qtlchenggrp/tianjing/GWAS/jwas.sbatch $repy $method
		cd ..
	done
	cd ..
done




