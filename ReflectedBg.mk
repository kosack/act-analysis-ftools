# rules for reflected-region background

MAKEREGIONS=$(PYTHON) $(TOOLSDIR)/regionbg.py $(ONTH2)

#-------------------------------------------------------------------------
# Region Background
#-------------------------------------------------------------------------
%_event_region_on.fits: %_event_selected.fits %_ON.reg
	@echo ON REGION EVENTS: $*
	@ftselect $< $@ 'regfilter("$*_ON.reg", DETX,DETY)' \
		clobber=true $(REDIRECT)

%_event_region_off.fits: %_event_selected.fits %_OFF.reg
	@echo OFF REGION EVENTS: $*
	@ftselect $< $@ 'regfilter("$*_OFF.reg", DETX,DETY)' \
		clobber=true $(REDIRECT)

#%_event_region_off.fits: %_event_selected.fits %_OFF.reg

# TODO: check for exclusion regions in regionbg.py! perhaps leave it
# the way it is, but add to the region filter the -regions from
# excluded.reg. Then calculate the correct area for ON and OFF using
# the flatlist! Not so good, since that would give some partial regions...

# TODO: naming option in regionbg.py, write Area on and Area off
# somewhere (then use 'fthedit run_023309_event_region_off.fits
# keyword=REGAREA value=13 operation=add' to add it to the header of
# selected events 
%_ON.reg: %_event_selected.fits
	@echo "REFLECTED REGIONS: $*"
	@$(MAKEREGIONS) $< $(REDIRECT)

%_OFF.reg: %_event_selected.fits
	@echo "REFLECTED REGIONS: $*"
	@$(MAKEREGIONS) $< $(REDIRECT)

# TODO: implement a tool that loops over ON/OFF event datasets and
# makes a FITS table with runwise statistics (N_on, N_off, Exp_on,
# Exp_off, Alpha, excess, signif), or alternately add these to the
# header of each region_on/off file (probably more flexible)
