# The template method just selects events from an OFF region in
# parameter space (Defined using the TEMPLATECUTS variable).  New
# acceptances must be calculated for OFF events, independently of ON
# events.


TEMPLATECUTS ?= '(HIL_MSW>2.0&&HIL_MSW<10.0)&&(HIL_MSL>-2.0&&HIL_MSL<2.0)&&MULTIP>=2'

TARGETS+=templatebg_excess_gauss.fits 
TARGETS+=templatebg_significance.fits

templatebg: $(TARGETS)



%_template_event_selected.fits: $(SOURCEDIR)/%_eventlist.fits.gz %_event_verify.txt
	@echo "TEMPLATE EVENT SELECTION: $*"
	@ftselect $< $@ $(TEMPLATECUTS) clobber=true $(REDIRECT)

RUNS_TEMPLATE_CMAP=$(addsuffix _template_cmap.fits,$(BASERUNS))	
RUNS_ACCMAP_TEMPLATE=$(addsuffix _template_accmap.fits,$(BASERUNS))

templatebg_cmap_sum.fits: $(RUNS_TEMPLATE_CMAP)
accmap_template_sum.fits: $(RUNS_ACCMAP_TEMPLATE)

templatebg_alpha.fits: accmap_sum.fits accmap_template_sum.fits
	@echo "TEMPLATE ALPHA MAP: $@"
	@ftpixcalc $@ 'A/B' a=accmap_sum.fits b=accmap_template_sum.fits \
		clobber=true $(REDIRECT)

templatebg_excess.fits: cmap_sum.fits templatebg_alpha.fits templatebg_cmap_sum.fits
	@echo "TEMPLATE EXCESS: $@"
	@ftpixcalc $@ 'A-B*C' \
		a=cmap_sum.fits \
		b=templatebg_cmap_sum.fits \
		c=templatebg_alpha.fits \
		clobber=yes $(REDIRECT)

templatebg_significance.fits: cmap_sum_tophat.fits templatebg_alpha_tophat.fits templatebg_cmap_sum_tophat.fits
	@echo "TEMPLATE SIGNIFICANCE: $@"
	@echo "DEBUG: $(LIMA)"
	@ftpixcalc alt_$@ '$(LIMA)' \
		a="NON=cmap_sum_tophat.fits" \
		b="NOFF=templatebg_cmap_sum_tophat.fits" \
		c="ALPHA=templatebg_alpha_tophat.fits" \
		clobber=yes $(REDIRECT)
	@ftpixcalc $@ '(NON-NOFF*ALPHA)/sqrt(ALPHA*(NON+NOFF))'\
		a="NON=cmap_sum_tophat.fits" \
		b="NOFF=templatebg_cmap_sum_tophat.fits" \
		c="ALPHA=templatebg_alpha_tophat.fits" \
		clobber=yes $(REDIRECT)




CLEANUP_RUNS+=$(RUNS_TEMPLATE_CMAP) $(RUNS_ACCMAP_TEMPLATE)
CLEANUP_EVENTS+=$(addsuffix _template_event_selected.fits,$(BASERUNS))
CLEANUP_BG+=templatebg_alpha.fits templatebg_excess.fits templatebg_cmap_sum.fits accmap_template_sum.fits