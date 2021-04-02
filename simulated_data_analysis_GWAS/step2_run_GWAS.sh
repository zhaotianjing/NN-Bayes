#! /bin/bash

nRepy=20

methodAll=("hmc_bayescpi_2" "bayescpi")

for repy in $( eval echo {1..$nRepy} )
do
	for method in "${methodAll[@]}"
	do
	    cd y"$repy"

		sbatch /group/qtlchenggrp/tianjing/GWAS/gwas.sbatch $repy $method

		cd .. 
	done
done




