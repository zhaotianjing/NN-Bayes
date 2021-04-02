#!/bin/bash

#this file is to calculate prediction accuracy 
repyAll=20 #y1,..,y20

methodAll=("hmc_rrblup_2" "hmc_rrblup_3" "hmc_rrblup_5" "hmc_rrblup_10" "rrblup"
	       "hmc_bayesa_2" "hmc_bayesa_3" "hmc_bayesa_5" "hmc_bayesa_10" "bayesa"
	       "hmc_bayesb_2" "hmc_bayesb_3" "hmc_bayesb_5" "hmc_bayesb_10" "bayesb"
	       "hmc_bayescpi_2" "hmc_bayescpi_3" "hmc_bayescpi_5" "hmc_bayescpi_10" "bayescpi"
	       "hmc_bayesl_2" "hmc_bayesl_3" "hmc_bayesl_5" "hmc_bayesl_10" "bayesl")

for method in "${methodAll[@]}"
do
	for repy in $( eval echo {1..$repyAll} )
	do
		cd y"$repy"/accuracy_"$method"
		sbatch /home/tianjing/paper_simu/accuracy.sbatch $repy $method
		cd ..
		cd ..
	done
done
