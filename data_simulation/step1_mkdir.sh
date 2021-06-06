#!/bin/bash

#this file is to create new folders for analysis

nRep=50  #rep1,..,rep50
repy=20  #y1,..,y20
methodAll=("hmc_rrblup_2" "hmc_rrblup_3" "hmc_rrblup_5" "hmc_rrblup_10" "rrblup"
	       "hmc_bayesa_2" "hmc_bayesa_3" "hmc_bayesa_5" "hmc_bayesa_10" "bayesa"
	       "hmc_bayesb_2" "hmc_bayesb_3" "hmc_bayesb_5" "hmc_bayesb_10" "bayesb"
	       "hmc_bayescpi_2" "hmc_bayescpi_3" "hmc_bayescpi_5" "hmc_bayescpi_10" "bayescpi"
	       "hmc_bayesl_2" "hmc_bayesl_3" "hmc_bayesl_5" "hmc_bayesl_10" "bayesl")  #e.g., hmc_bayescpi_2 is for NN-BayesCPi-2


for repy in $( eval echo {1..$repy} )   #y1..y20
do
    mkdir -p y"$repy"
    cd y"$repy"
    for method in "${methodAll[@]}"     #e.g., "NN-BayesCpi-2"
    do
	    mkdir -p "$method"
		cd "$method"
	    for j in $( eval echo {1..$nRep} ) 
	    do
	      mkdir -p rep"$j"               #rep1,...,rep50
	    done
		cd ..

		mkdir -p run_"$method"           #folder to run analysis
		mkdir -p accuracy_"$method"      #folder to calculate accuracy
	done
	cd ..
done