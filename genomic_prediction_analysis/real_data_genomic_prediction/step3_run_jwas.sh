#! /bin/bash

#this file is to run JWAS for one species
species="maize"
nRep=20 #rep1,...,rep20

methodAll=("hmc_rrblup_2" "hmc_rrblup_3" "hmc_rrblup_5" "hmc_rrblup_10" "rrblup"
		   "hmc_bayesa_2" "hmc_bayesa_3" "hmc_bayesa_5" "hmc_bayesa_10" "bayesa" 
	       "hmc_bayesb_2" "hmc_bayesb_3" "hmc_bayesb_5" "hmc_bayesb_10" "bayesb" 
	       "hmc_bayescpi_2" "hmc_bayescpi_3" "hmc_bayescpi_5" "hmc_bayescpi_10" "bayescpi"
	       "hmc_bayesl_2" "hmc_bayesl_3" "hmc_bayesl_5" "hmc_bayesl_10" "bayesl")  #e.g., hmc_bayescpi_2 is NN-BayesCPi-2

cd "$species"

for method in "${methodAll[@]}"
do
	cd run_"$method"

	for rep in $( eval echo {1..$nRep} )  #rep1,...,rep20
	do
	    sbatch /group/qtlchenggrp/tianjing/paper_real/jwas.sbatch $species $method $rep
	done

	cd .. 

done



