#!/bin/sh

# just a hack right now, should do this with a proper makefile:
# generates the final lookup tables from the sums (see sum-lookups.sh
# and make-runwise-lookups.sh)

TOOLSDIR=~/Source/Working/PyFITSTools



#for tel in CT001 CT002 CT003 CT004; do
for tel in TYPE01_00; do
    for val in HIL_TEL_WIDTH HIL_TEL_LENGTH MC_ENERGY; do

	echo LOOKUPS FOR CT $tel $val  ...
	
	# mean and stddev:

	ftpixcalc ${tel}-${val}-lookup-mean.fits 'A/(B+1e-10)' \
	    a=${tel}-${val}-lookup-sum.fits \
	    b=${tel}-${val}-lookup-count.fits clobber=true

	ftpixcalc ${tel}-${val}-lookup-stddev.fits 'sqrt(C/(B+1e-10)-A**2)' \
	    a=${tel}-${val}-lookup-mean.fits \
	    b=${tel}-${val}-lookup-count.fits \
	    c=${tel}-${val}-lookup-sum2.fits  clobber=true


    done
done
