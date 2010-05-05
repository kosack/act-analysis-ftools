# The template method just selects events from an OFF region in
# parameter space (Defined using the TEMPLATECUTS variable).  New
# acceptances must be calculated for OFF events, independently of ON
# events.

TEMPLATECUTS ?= '(HIL_MSW>2.0&&HIL_MSW<10.0)&&(HIL_MSL>-2.0&&HIL_MSL<2.0)&&MULTIP>=2'
%_template_event_selected.fits: $(SOURCEDIR)/%_eventlist.fits.gz %_event_verify.txt
	@echo "TEMPLATE EVENT SELECTION: $*"
	@ftselect $< $@ $(TEMPLATECUTS) clobber=true $(REDIRECT)

RUNS_TEMPLATE_CMAP=$(addsuffix _template_cmap.fits,$(BASERUNS))	
RUNS_ACCMAP_TEMPLATE=$(addsuffix _template_accmap.fits,$(BASERUNS))

template_cmap_sum.fits: $(RUNS_TEMPLATE_CMAP)
accmap_template_sum.fits: $(RUNS_ACCMAP_TEMPLATE)

template_alpha.fits: accmap_sum.fits accmap_template_sum.fits
	@echo "TEMPLATE ALPHA MAP: $@"
	@ftpixcalc $@ 'A/B' a=accmap_sum.fits b=accmap_template_sum.fits \
		clobber=true $(REDIRECT)

template_excess.fits: cmap_sum.fits template_alpha.fits template_cmap_sum.fits
	@echo "TEMPLATE EXCESS: $@"
	@ftpixcalc $@ 'A-B*C' \
		a=cmap_sum.fits \
		b=template_cmap_sum.fits \
		c=template_alpha.fits \
		clobber=yes $(REDIRECT)

TARGETS+=template_excess_gauss.fits
CLEANUP_RUNS+=$(RUNS_TEMPLATE_CMAP) $(RUNS_ACCMAP_TEMPLATE)
CLEANUP_EVENTS+=$(addsuffix _template_event_selected.fits,$(BASERUNS))
CLEANUP_BG+=template_alpha.fits template_excess.fits template_cmap_sum.fits accmap_template_sum.fits