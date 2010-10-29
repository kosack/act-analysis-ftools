#!/bin/sh

for ii in ~/Data/FITS/HESS/Simulations/Phase1b/Gamma/0.7deg/*.fits.gz; do
    python ${TOOLSDIR}/generate-lookup-tables.py $ii
done

