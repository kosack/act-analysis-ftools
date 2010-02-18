#!/bin/sh

# just a hack right now, should do this with a proper makefile:
# generates the final lookup tables from the sums (see sum-lookups.sh
# and make-runwise-lookups.sh)

TOOLSDIR=~/Source/Working/PyFITSTools



for tel in 001 002 003 004; do
    for val in HIL_TEL_WIDTH HIL_TEL_LENGTH MC_ENERGY; do

	echo LOOKUPS FOR CT$tel $val  ...
	
	# mean and stddev:

	ftpixcalc CT${tel}-${val}-lookup-mean.fits 'A/(B+1e-10)' \
	    a=CT${tel}-${val}-lookup-sum.fits \
	    b=CT${tel}-${val}-lookup-count.fits clobber=true

	ftpixcalc CT${tel}-${val}-lookup-stddev.fits 'sqrt(C/(B+1e-10)-A**2)' \
	    a=CT${tel}-${val}-lookup-mean.fits \
	    b=CT${tel}-${val}-lookup-count.fits \
	    c=CT${tel}-${val}-lookup-sum2.fits  clobber=true

	# gaussian smoothing:

	gaussargs="2.0 nullval=1e-5 ratio=1.0 theta=0.0  nsigma=4.0"
	gaussargs="${gaussargs} boundary=nearest constant=0.0 clobber=yes"

	fgauss CT${tel}-${val}-lookup-mean.fits tmp1.fits $gaussargs
	mv tmp1.fits CT${tel}-${val}-lookup-mean-gauss.fits

	fgauss CT${tel}-${val}-lookup-stddev.fits  tmp2.fits $gaussargs
	mv tmp2.fits CT${tel}-${val}-lookup-stddev-gauss.fits

        fgauss CT${tel}-${val}-lookup-count.fits tmp3.fits $gaussargs
	mv  tmp3.fits CT${tel}-${val}-lookup-count-gauss.fits

    done
done
