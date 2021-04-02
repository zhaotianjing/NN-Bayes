#!/bin/bash

#this file is to make folder for following analysis
nRep=20
species=("maize" "pig" "soy" "spruce" "switchgrass")
methodAll=("hmc_rrblup_2" "hmc_rrblup_3" "hmc_rrblup_5" "hmc_rrblup_10" "rrblup"
	       "hmc_bayesa_2" "hmc_bayesa_3" "hmc_bayesa_5" "hmc_bayesa_10" "bayesa"
	       "hmc_bayesb_2" "hmc_bayesb_3" "hmc_bayesb_5" "hmc_bayesb_10" "bayesb"
	       "hmc_bayescpi_2" "hmc_bayescpi_3" "hmc_bayescpi_5" "hmc_bayescpi_10" "bayescpi"
	       "hmc_bayesl_2" "hmc_bayesl_3" "hmc_bayesl_5" "hmc_bayesl_10" "bayesl")

for name in "${species[@]}"
do
    mkdir -p "$name"
    cd "$name"
    mkdir -p train_test_ID

    for method in "${methodAll[@]}"
    do
	    mkdir -p "$method"
		cd "$method"
	    for j in $( eval echo {21..$nRep} ) 
	    do
	      mkdir -p rep"$j"
	    done
		cd ..

		mkdir -p run_"$method"
		mkdir -p accuracy_"$method"
	done
	cd ..
done