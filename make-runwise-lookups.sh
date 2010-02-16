#!/bin/sh

TOOLSDIR=~/Source/Working/PyFITSTools

for ii in ~/Data/FITS/HESS_Simulations/*.fits.gz; do
    python ${TOOLSDIR}/generate-lookup-tables.py $ii
done

