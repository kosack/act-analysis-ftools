# =========================================================================
# ANALYSIS PARAMETERS:
# =========================================================================
CUTS='(HIL_MSW>-2.0&&HIL_MSW<0.7)&&(HIL_MSL>-2.0&&HIL_MSL<2.0)'
SOURCEDIR=$(HOME)/Analysis/FITSEventLists/HESS_Crab
EXCLMASK='regfilter("excluded.reg",RA,DEC)'
FOVX=7.0
FOVY=7.0
GEOMX=301
GEOMY=301
CENTERRA=83.633333
CENTERDEC=22.014444

# =========================================================================
# runlist
# =========================================================================
LISTS=$(wildcard $(SOURCEDIR)/*.fits.gz)
BASERUNS=$(patsubst %_eventlist.fits,%,$(basename $(notdir $(LISTS))))
RUNS=$(addsuffix _accmap.fits,$(BASERUNS))
RUNS_CMAP=$(addsuffix _cmap.fits,$(BASERUNS))

# =========================================================================
# utility parameters
# =========================================================================
MAPARGS=--fov $(FOVX),$(FOVY) --geom $(GEOMX),$(GEOMY) \
	--center $(CENTERRA),$(CENTERDEC)

MAKEFITS=$(HOME)/Source/PyFITSTools/makefits.py
ACCEPTANCE=$(HOME)/Source/PyFITSTools/acceptance.py


# =========================================================================
# Rules to generate various outputs:
# =========================================================================


.SECONDARY: # clear secondary rule, so intermediate files aren't deleted
	echo "Secondary $@"

all: $(RUNS)
	@echo "Processing runs"

# Gamma-hadron separated eventlist
%_event_selected.fits: $(SOURCEDIR)/%_eventlist.fits.gz
	@echo *** EVENT SELECTION $*
	ftselect $< $@ $(CUTS)

# masked eventlist
%_event_excluded.fits: %_event_selected.fits
	@echo ===========================================
	@echo EVENT EXCLUSION $*
	@echo ===========================================
	ftselect $< $@ $(EXCLMASK)

# countmap 
%_cmap.fits: %_event_selected.fits
	@echo ===========================================
	@echo COUNT MAP $*
	@echo ===========================================
	python $(MAKEFITS) $(MAPARGS) --output $@ $<

#excluded count map
%_cmap_excluded.fits: %_event_excluded.fits
	@echo ===========================================
	@echo  EXCLUDED COUNT MAP $*
	@echo ===========================================
	python $(MAKEFITS) $(MAPARGS) --output $@ $<

# acceptance map
%_accmap.fits: %_event_excluded.fits %_cmap.fits 
	@echo ===========================================
	@echo ACCEPTANCE MAP $*
	@echo ===========================================
	python $(ACCEPTANCE) --output $@ $^

sum_cmap.fits: $(RUNS_CMAP)
	@echo "RUNS: '$(RUNS_CMAP)'"
	for ii in $(RUNS_CMAP); do echo $$ii ;done

# do 
#   echo $$ii
#   ftpixcalc sum_cmap.tmp.fits 'A+B' a=$@ b=$$ii 
#   rm -f $@
#   mv sum_cmap.tmp.fits $@
# done

clean:
	 rm -fv run_*_*.fits