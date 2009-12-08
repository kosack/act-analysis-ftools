full-analysis:
	process_runs.pl


testflat:
	rm -fv tmp*.fits
	python2.4 make-flat-eventlist.py \
		~/Analysis/FITSEventLists/Analysis/cmap_sum.fits tmpflat.fits
	ftselect tmpflat.fits tmpexcl.fits 'regfilter("excluded.reg",RA,DEC)'
	python2.4 makefits.py --output tmpflatmap.fits --geom 256,256 --fov 6.0,6.0 tmpflat.fits
	python2.4 makefits.py --output tmpexclmap.fits --geom 256,256 --fov 6.0,6.0 tmpexcl.fits

