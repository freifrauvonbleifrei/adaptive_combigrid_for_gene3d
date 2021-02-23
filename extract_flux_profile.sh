#!/usr/bin/env bash
#module unload python-site/3.6
#module load tools/conda/default #cray-python
module load python/3.6_intel

PROB_PREPATH="/hppfs/scratch/02/di39qun2/gene3d-flw-simulations"
DIAG_PATH="/hppfs/work/pn34mi/di39qun2/gene-diag"
DIAGNOSTICS=DiagFluxesgene3d
for PROBNAME in $PROB_PREPATH/prob_*; do
        PROBNAME=$(basename $PROBNAME)
	#PROBNAME='prob_5_5_5_5_3'
	#PROBNAME='kjm_prob_5_7_8_5_3'
	if [ -d $PROB_PREPATH/${PROBNAME} ]; then
		echo $PROBNAME
		#totallasttime=0.0
		#for f in `ls $PROBNAME/nrg*`; do
		#    lasttime=$(tail -n 2 $f | head -n 1)
		#    echo "lasttime of file $f is $lasttime"
		#    if (( $(echo "$lasttime > $totallasttime" |bc -l) )); then
		#	endtime=$(tail -n 4 $f |head -n 1)
		#    fi
		#done
		OUTDIR='./flux_diags/'
		#python3 ./diag-python/GENE_gui_python/gene_cl.py -i ./$PROBNAME/out -o $OUTDIR --starttime=800.0 --endtime=$endtime -e '_1' '_2' '_3' '_4' '_5' '_6' '_7' '_8' '_9' '_10' '_11' '_12' '_13' '_14' '_15' '_16' '_17' '_18' '_19' '_20' '_21' '_22' '_23' '.dat'
		#python3 ./diag-python/GENE_gui_python/gene_cl.py -i $PROB_PREPATH/$PROBNAME/out -o $OUTDIR --starttime=150.0 -d${DIAGNOSTICS} -e '_1' '_2' '_3' '_4' '_5' '_6' '_7' '_8' '_9' '_10' '_11' '_12' '_13' '_14' '_15' '_16' '_17' '_18' '_19' '_20' '_21' '_22' '_23' #'.dat'
		python3 ${DIAG_PATH}/GENE_gui_python/gene_cl.py -i $PROB_PREPATH/$PROBNAME/out -o $OUTDIR --starttime=150.0 -d${DIAGNOSTICS} -e '.dat'
		
        mv $OUTDIR/flux_profile_ions_.h5 $OUTDIR/flux_profile_ions_$PROBNAME.h5
        mv $OUTDIR/flux_profile_electrons_.h5 $OUTDIR/flux_profile_electrons_$PROBNAME.h5
		RETURNCODE=$?
		if [ $RETURNCODE -ne 0 ]; then
			echo $PROBNAME >> $OUTDIR/probnames_that_did_not_work
		fi
	else 
		echo $PROBNAME not there!
	fi
done
