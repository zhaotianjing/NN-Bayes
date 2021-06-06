#! /bin/bash

#this file is to calculate prediction accuracy for one species
species="maize"
methodAll=("hmc_rrblup_2" "hmc_rrblup_3" "hmc_rrblup_5" "hmc_rrblup_10" "rrblup"
		   "hmc_bayesa_2" "hmc_bayesa_3" "hmc_bayesa_5" "hmc_bayesa_10" "bayesa" 
	       "hmc_bayesb_2" "hmc_bayesb_3" "hmc_bayesb_5" "hmc_bayesb_10" "bayesb" 
	       "hmc_bayescpi_2" "hmc_bayescpi_3" "hmc_bayescpi_5" "hmc_bayescpi_10" "bayescpi"
	       "hmc_bayesl_2" "hmc_bayesl_3" "hmc_bayesl_5" "hmc_bayesl_10" "bayesl")

nRep=20  #rep1,...,rep20

cd "$species"
for method in "${methodAll[@]}"
do
	cd accuracy_"$method"
	sbatch /group/qtlchenggrp/tianjing/paper_real/accuracy.sbatch $species $method
	cd ..  

done
