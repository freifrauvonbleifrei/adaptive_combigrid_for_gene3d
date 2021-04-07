#!/usr/bin/env bash

. ./adaptation_parameters.sh

PROB_PREPATH=${ADAPTATION_PROB_PATH}
DIAG_PATH="./gene-diag"
export PYTHONPATH=$DIAG_PATH/GENE_gui_python:$PYTHONPATH

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
            
            NRGFILESLIST=""
            # https://unix.stackexchange.com/questions/239772/bash-iterate-file-list-except-when-empty
            shopt -s nullglob
            for NRGFILENUM in $PROB_PREPATH/$PROBNAME/out/nrg_*; do
                echo $NRGFILENUM
                NRGFILENUM=$(basename $NRGFILENUM)
                NRGSUBSTRING=$(echo $NRGFILENUM | cut -d'_' -f 2)
                NRGFILESLIST+=" _${NRGSUBSTRING}"
                #echo $NRGFILESLIST
            done
            if [ -e $PROB_PREPATH/$PROBNAME/out/nrg.dat ] ; then
                NRGFILESLIST+=" .dat"
            fi
            
            #python3 ${DIAG_PATH}/GENE_gui_python/gene_cl.py -i $PROB_PREPATH/$PROBNAME/out -o $OUTDIR --starttime=${ADAPTATION_PARAM_CROP_TIME} -d${DIAGNOSTICS} -e ${NRGFILESLIST}
            python3 ./additional-diagnostics/additional_diagnostics.py -i $PROB_PREPATH/$PROBNAME/out -o $OUTDIR --starttime=${ADAPTATION_PARAM_CROP_TIME} -d${DIAGNOSTICS} -e ${NRGFILESLIST}

            #mv $OUTDIR/flux_profile_ions_.h5 $OUTDIR/flux_profile_ions__$PROBNAME.h5
            mv $OUTDIR/flux_profile_ions.h5 $OUTDIR/flux_profile_ions_$PROBNAME.h5
            mv $OUTDIR/flux_spectra_Gem_ions.h5 $OUTDIR/flux_spectra_Gem_ions_$PROBNAME.h5
            mv $OUTDIR/flux_spectra_Ges_ions.h5 $OUTDIR/flux_spectra_Ges_ions_$PROBNAME.h5
            mv $OUTDIR/flux_spectra_Qem_ions.h5 $OUTDIR/flux_spectra_Qem_ions_$PROBNAME.h5
            mv $OUTDIR/flux_spectra_Qes_ions.h5 $OUTDIR/flux_spectra_Qes_ions_$PROBNAME.h5
            mv $OUTDIR/profile_ions.h5 $OUTDIR/profile_ions_$PROBNAME.h5
            if [ ${ADAPTATION_NUMBER_OF_SPECIES} = "2" ] ; then
                #mv $OUTDIR/flux_profile_electrons_.h5 $OUTDIR/flux_profile_electrons__$PROBNAME.h5
                mv $OUTDIR/flux_profile_electrons.h5 $OUTDIR/flux_profile_electrons_$PROBNAME.h5
                mv $OUTDIR/flux_spectra_Gem_electrons.h5 $OUTDIR/flux_spectra_Gem_electrons_$PROBNAME.h5
                mv $OUTDIR/flux_spectra_Ges_electrons.h5 $OUTDIR/flux_spectra_Ges_electrons_$PROBNAME.h5
                mv $OUTDIR/flux_spectra_Qem_electrons.h5 $OUTDIR/flux_spectra_Qem_electrons_$PROBNAME.h5
                mv $OUTDIR/flux_spectra_Qes_electrons.h5 $OUTDIR/flux_spectra_Qes_electrons_$PROBNAME.h5
                mv $OUTDIR/profile_electrons.h5 $OUTDIR/profile_electrons_$PROBNAME.h5

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
