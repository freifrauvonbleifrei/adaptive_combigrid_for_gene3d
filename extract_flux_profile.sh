#!/usr/bin/env bash

. ./adaptation_parameters.sh

PROB_PREPATH=${ADAPTATION_PROB_PATH}
DIAG_PATH="./gene-diag"
DIAGNOSTICS=DiagFluxesgene3d
OUTDIR='./flux_diags/'
mkdir -p $OUTDIR
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
        #       endtime=$(tail -n 4 $f |head -n 1)
        #    fi
        #done
        # only re-run the diagnostic if the nrg.dat file is newer than the existing flux profile output
        if [ ! -e $OUTDIR/flux_profile_ions_$PROBNAME.h5 ] || [ $PROB_PREPATH/$PROBNAME/out/nrg.dat -nt $OUTDIR/flux_profile_ions_$PROBNAME.h5 ] ; then
            #python3 ./diag-python/GENE_gui_python/gene_cl.py -i ./$PROBNAME/out -o $OUTDIR --starttime=800.0 --endtime=$endtime -e '_1' '_2' '_3' '_4' '_5' '_6' '_7' '_8' '_9' '_10' '_11' '_12' '_13' '_14' '_15' '_16' '_17' '_18' '_19' '_20' '_21' '_22' '_23' '.dat'
            NRGFILESLIST=""
            # https://unix.stackexchange.com/questions/239772/bash-iterate-file-list-except-when-empty
            shopt -s nullglob
            for NRGFILENUM in $PROB_PREPATH/$PROBNAME/out/nrg_*; do
                echo $NRGFILENUM
                NRGFILENUM=$(basename $NRGFILENUM)
                NRGSUBSTRING=$(echo $NRGFILENUM | cut -d'_' -f 2)
                NRGFILESLIST+=" _${NRGSUBSTRING}"
            done
            NRGFILESLIST+=" .dat"
            #python3 ./diag-python/GENE_gui_python/gene_cl.py -i $PROB_PREPATH/$PROBNAME/out -o $OUTDIR --starttime=150.0 -d${DIAGNOSTICS} -e '_1' '_2' '_3' '_4' '_5' '_6' '_7' '_8' '_9' '_10' '_11' '_12' '_13' '_14' '_15' '_16' '_17' '_18' '_19' '_20' '_21' '_22' '_23' #'.dat'
            python3 ${DIAG_PATH}/GENE_gui_python/gene_cl.py -i $PROB_PREPATH/$PROBNAME/out -o $OUTDIR --starttime=150.0 -d${DIAGNOSTICS} -e ${NRGFILESLIST}

            mv $OUTDIR/flux_profile_ions_.h5 $OUTDIR/flux_profile_ions_$PROBNAME.h5
            if [ ${ADAPTATION_NUMBER_OF_SPECIES} = "2" ] ; then
                mv $OUTDIR/flux_profile_electrons_.h5 $OUTDIR/flux_profile_electrons_$PROBNAME.h5
            fi
        fi
        RETURNCODE=$?
        if [ $RETURNCODE -ne 0 ]; then
            echo $PROBNAME >>$OUTDIR/flux_probnames_that_did_not_work
        fi
    else
        echo $PROBNAME not there!
    fi
done
