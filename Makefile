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


REDIRECT= >> output.log 2>&1 # set to blank to get all output

# =========================================================================
# runlist
# =========================================================================
LISTS=$(wildcard $(SOURCEDIR)/*.fits.gz)
BASERUNS=$(patsubst %_eventlist.fits,%,$(basename $(notdir $(LISTS))))

# =========================================================================
# utility parameters
# =========================================================================
MAPARGS=--fov $(FOVX),$(FOVY) --geom $(GEOMX),$(GEOMY) \
	--center $(CENTERRA),$(CENTERDEC)

PYTHON=python
TOOLSDIR=$(HOME)/Source/PyFITSTools
MAKEFITS=$(PYTHON) $(TOOLSDIR)/makefits.py
ACCEPTANCE=$(PYTHON) $(TOOLSDIR)/acceptance.py
SUMMER=$(TOOLSDIR)/sum_maps.pl
FLATLIST=$(PYTHON) $(TOOLSDIR)/make-flat-eventlist.py -s 1 

# =========================================================================
# Rules to generate various outputs:
# =========================================================================


.SECONDARY: # clear secondary rule, so intermediate files aren't deleted
	echo "Secondary $@"

all: fov_excess.fits
	@echo "Done processing runs"

# Gamma-hadron separated eventlist
%_event_selected.fits: $(SOURCEDIR)/%_eventlist.fits.gz
	@echo EVENT SELECTION $*
	@ftselect $< $@ $(CUTS) $(REDIRECT)

# masked eventlist
%_event_excluded.fits: %_event_selected.fits
	@echo EVENT EXCLUSION $*
	@ftselect $< $@ $(EXCLMASK) $(REDIRECT)

# countmap 
%_cmap.fits: %_event_selected.fits
	@echo COUNT MAP $*
	@$(MAKEFITS) $(MAPARGS) --output $@ $< $(REDIRECT)

#excluded count map
%_cmap_excluded.fits: %_event_excluded.fits
	@echo  EXCLUDED COUNT MAP $*
	@$(MAKEFITS) $(MAPARGS) --output $@ $< $(REDIRECT)

# acceptance map
%_accmap.fits: %_event_excluded.fits %_cmap.fits 
	@echo ACCEPTANCE MAP $*
	@$(ACCEPTANCE) --output $@ $^ $(REDIRECT)

# ==============================================================
# Summed maps
# ==============================================================

RUNS_CMAP=$(addsuffix _cmap.fits,$(BASERUNS))
RUNS_CMAP_EXCL=$(addsuffix _cmap_excluded.fits,$(BASERUNS))
RUNS_ACCMAP=$(addsuffix _accmap.fits,$(BASERUNS))

sum_cmap.fits: $(RUNS_CMAP)
	@$(RM) $@
	@$(SUMMER) -v -o $@ $^ 


sum_cmap_excluded.fits: $(RUNS_CMAP_EXCL)
	@$(RM) $@
	@$(SUMMER) -v -o $@ $^


sum_accmap.fits: $(RUNS_ACCMAP)
	@$(RM) $@
	@$(SUMMER) -v -o $@ $^


# ==============================================================
# Final maps
# ==============================================================

# exclusion map
FIRSTCMAP=$(firstword $(RUNS_CMAP))

flatlist.fits: $(FIRSTCMAP)
	@echo "FLAT EVENTLIST: $@ using $^"
	@$(FLATLIST) $(FIRSTCMAP) $@ $(REDIRECT)

flatlist_excluded.fits: flatlist.fits
	@echo "FLAT EVENTLIST EXCLUDED: $@"
	@ftselect $< $@ $(EXCLMASK) $(REDIRECT)	


exclmap.fits: flatlist_excluded.fits
	@echo "EXCLUSION MAP: $@"
	@$(MAKEFITS) $(MAPARGS) --output $@ $< $(REDIRECT)

flatmap.fits: flatlist.fits
	@echo "FLAT MAP: $@"
	@$(MAKEFITS) $(MAPARGS) --output $@ $< $(REDIRECT)


# Field-of-view background model map (background is assumed to be the
# scaled acceptance map)
fov_excess.fits: sum_accmap.fits sum_cmap.fits
	@echo FOV EXCESS MAP: $@
	@ftpixcalc $@ 'A-B' a=sum_cmap.fits b=sum_accmap.fits \
		clobber=true $(REDIRECT)


# ==============================================================
# cleanup
# ==============================================================

clean: clean-runs clean-sums clean-fov clean-excl
	$(RM) output.log

clean-runs:
	$(RM) run_*_*.fits

clean-sums:
	$(RM) sum_*.fits

clean-fov:
	$(RM) fov_excess.fits	

clean-excl:
	$(RM) exclmap.fits flatmap.fits