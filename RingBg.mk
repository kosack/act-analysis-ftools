#-------------------------------------------------------------------------
# ring background model (off events are sum of on events in annulus about
# each bin)
#-------------------------------------------------------------------------

RINGAREAFACTOR ?= 7.0
RINGGAP ?= 0.2                # gap between ring and ON radius in degrees

TARGETS += ringbg_excess_gauss.fits 
TARGETS += ringbg_significance.fits 
TARGETS += ringbg_diagnostic_significance.ps

ringbg: $(TARGETS)


ring.fits: exclmap.fits
	@echo "GENERATING RING: $@"
	@$(MAKERING) --output $@ --onradius=$(strip $(ONRADIUS)) \
		--gap $(RINGGAP) \
		--areafactor $(RINGAREAFACTOR) $^ $(REDIRECT)


exclmap_ring.fits: exclmap.fits ring.fits
	@echo "CONVOLVE RING EXCLMAP: $@"
	@$(CONVOLVE) --output $@ $^ $(REDIRECT)

flatmap_ring.fits: flatmap.fits ring.fits
	@echo "CONVOLVE RING FLATMAP: $@"
	@$(CONVOLVE) --output $@ $^ $(REDIRECT)

# the RingBg off map is the sum of events in a ring about each
# position, taking into account exclusions.  Therefore it's the
# convolution of the excluded countmap with a ring, with the
# ring-convolved exclusion map used as a correction.
%_offmap_ring.fits: %_cmap_excluded.fits ring.fits exclmap_ring.fits
	@echo "RING OFF MAP: $*"
	@$(CONVOLVE) --output tmp_$@ $^ $(REDIRECT)
	@ftimgcalc $@ 'A*sum(C)/B' \
		a=tmp_$@ \
		b=exclmap_ring.fits \
		c=ring.fits \
		clobber=true $(REDIRECT)
	@$(RM) tmp_$@


%_accmap_ring.fits: %_accmap.fits ring.fits flatmap_ring.fits
	@echo "CONVOLVE RING ACCMAP: $*"
	@$(CONVOLVE) --output $@ $^ $(REDIRECT)
#	@$(CONVOLVE) --output tmp_$@ $^ $(REDIRECT)
#	@ftpixcalc $@ 'A/B' a=tmp_$@ b=flatmap_ring.fits clobber=true $(REDIRECT)
#	@$(RM) tmp_$@


ringbg_alpha.fits: accmap_sum.fits accmap_ring_sum.fits
	@echo "RING ALPHA MAP"
	@ftpixcalc $@ 'A/B' a=accmap_sum.fits b=accmap_ring_sum.fits \
		clobber=true $(REDIRECT)


ringbg_excess.fits: cmap_sum.fits ringbg_alpha.fits offmap_ring_sum.fits
	@echo "RING EXCESS: $@"
	@ftpixcalc $@ 'A-B*C' \
		a=cmap_sum.fits \
		b=offmap_ring_sum.fits \
		c=ringbg_alpha.fits \
		clobber=yes $(REDIRECT)


ringbg_alpha_tophat.fits: accmap_sum_tophat.fits accmap_ring_sum_tophat.fits
	@echo "RING ALPHA MAP TOPHAT: $@"
	@ftpixcalc $@ 'A/B' a=accmap_sum_tophat.fits b=accmap_ring_sum_tophat.fits \
		clobber=true $(REDIRECT)
ringbg_significance.fits: cmap_sum_tophat.fits ringbg_alpha_tophat.fits offmap_ring_sum_tophat.fits
	@echo "RING SIGNIFICANCE: $@"
	@echo "DEBUG: $(LIMA)"
	@ftpixcalc alt_$@ '$(LIMA)' \
		a="NON=cmap_sum_tophat.fits" \
		b="NOFF=offmap_ring_sum_tophat.fits" \
		c="ALPHA=ringbg_alpha_tophat.fits" \
		clobber=yes $(REDIRECT)
	@ftpixcalc $@ '$(SIGNIF)' \
		a="NON=cmap_sum_tophat.fits" \
		b="NOFF=offmap_ring_sum_tophat.fits" \
		c="ALPHA=ringbg_alpha_tophat.fits" \
		clobber=yes $(REDIRECT)

# TODO: need to divide by exclmap_tophat to get the right values in the exclusion region!
EXCLCORR=(EXCL/max(EXCL))

ringbg_significance_exmasked.fits: cmap_sum_exmasked_tophat.fits ringbg_alpha_exmasked_tophat.fits offmap_ring_sum_exmasked_tophat.fits exclmap_tophat.fits
	@echo "RING MASKED SIGNIFICANCE: $@"
	@echo "LIMA_EXCL: $(LIMA_EXCL)"
	@ftpixcalc alt_$@ '$(LIMA_EXCL)' \
		a="NON=cmap_sum_exmasked_tophat.fits" \
		b="NOFF=offmap_ring_sum_exmasked_tophat.fits" \
		c="ALPHA=ringbg_alpha_exmasked_tophat.fits" \
		d="EXCL=exclmap_tophat.fits" \
		clobber=yes $(REDIRECT)
	@ftimgcalc $@ '$(SIGNIF_EXCL)'\
		a="NON=cmap_sum_exmasked_tophat.fits" \
		b="NOFF=offmap_ring_sum_tophat.fits" \
		c="ALPHA=ringbg_alpha_tophat.fits" \
		d="EXCL=exclmap_tophat.fits" \
		clobber=yes $(REDIRECT)


#-------------------------------------------------------------------------
# Diagnostic plots (some use gnuplot)
#-------------------------------------------------------------------------

ringbg_diagnostic_significance.ps: $(TOOLSDIR)/ringbg_diagnostic_significance.gpl ringbg_significance_imhist.fits ringbg_significance_exmasked_imhist.fits
	@echo "DIAGNOSTIC: $@"
	gnuplot $< > $@


%_imhist.fits: %.fits
	@echo "HISTOGRAM: $@"
	@fimhisto $< $@ range=-10,100 binsize=1.0 clobber=yes  $(REDIRECT)
