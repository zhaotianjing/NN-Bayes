#! /bin/bash

#this file is to calculate AUC for linear and non-linear models
nRepy=20  #y1,...,y20
methodAll=("hmc_bayescpi_2" "bayescpi")

for repy in $( eval echo {1..$nRepy} ) #y1,...,y20
do
    cd y"$repy"
    for method in "${methodAll[@]}"
    do  
    	cd "$method"
		sbatch /group/qtlchenggrp/tianjing/GWAS/roc.sbatch $repy $method
		cd ..
    done
	cd ..
done




