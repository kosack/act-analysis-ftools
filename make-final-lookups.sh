#!/bin/sh

# just a hack right now, should do this with a proper makefile:
# generates the final lookup tables from the sums (see sum-lookups.sh
# and make-runwise-lookups.sh)

TOOLSDIR=~/Source/Working/PyFITSTools



for tel in 001 002 003 004; do
    for val in HIL_TEL_WIDTH HIL_TEL_LENGTH; do

	echo LOOKUPS FOR CT$tel $val  ...
	
	ftpixcalc CT${tel}-${val}-lookup-mean.fits 'A/(B+1e-10)' \
	    a=CT${tel}-${val}-lookup-sum.fits \
	    b=CT${tel}-${val}-lookup-count.fits clobber=true

	fgauss CT${tel}-${val}-lookup-mean.fits \
	    CT${tel}-${val}-lookup-mean-gauss.fits sigma=2.0 \
	    nullval=1e-10 ratio=1.0 theta=0.0 nsigma=4.0 boundary=nearest \
	    constant=0.0 clobber=true

	ftpixcalc CT${tel}-${val}-lookup-stddev.fits 'sqrt(C/(B+1e-10)-A**2)' \
	    a=CT${tel}-${val}-lookup-mean.fits \
	    b=CT${tel}-${val}-lookup-count.fits \
	    c=CT${tel}-${val}-lookup-sum2.fits  clobber=true


	fgauss CT${tel}-${val}-lookup-stddev.fits \
	    CT${tel}-${val}-lookup-stddev-gauss.fits sigma=2.0 \
	    nullval=1e-10 ratio=1.0 theta=0.0 nsigma=4.0 boundary=nearest \
	    constant=0.0 clobber=true

	fgauss CT${tel}-${val}-lookup-count.fits \
	    CT${tel}-${val}-lookup-count-gauss.fits sigma=2.0 \
	    nullval=1e-10 ratio=1.0 theta=0.0 nsigma=4.0 boundary=nearest \
	    constant=0.0 clobber=true


    done
done
