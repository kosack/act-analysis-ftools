# =========================================================================
# Set defaults:
# =========================================================================

TOOLSDIR ?= $(HOME)/Source/PyFITSTools
CUTS ?= '(HIL_MSW>-2.0&&HIL_MSW<0.7)&&(HIL_MSL>-2.0&&HIL_MSL<2.0)'
SOURCEDIR ?=$(HOME)/Analysis/FITSEventLists/HESS_Crab
FOVX ?= 7.0                   # Field of view of output map (X degrees)
FOVY ?= 7.0                   # Field of view of output map (Y degrees)
GEOMX ?= 301                  # number of X bins in map (integer)
GEOMY ?= 301		      # number of Y bins in map (integer)
CENTERRA ?= 83.633333	      # center of output map in RA
CENTERDEC ?= 22.014444	      # center of output map in Dec
EXCLUSIONFILE ?= excluded.reg # exclusion region file in ascii region format

# =========================================================================
# runlist
# =========================================================================
BASERUNS=$(patsubst %_eventlist.fits,%,$(basename $(notdir $(EVENTLISTS))))

# =========================================================================
# utility parameters
# =========================================================================
REDIRECT= >> output.log 2>&1 # set to blank to get all output

EXCLMASK='regfilter("$(strip $(EXCLUSIONFILE))",RA,DEC)' # spatial exclusion filter

# parameters passed to MAKEMAP to generate the images
MAPARGS=--fov $(strip $(FOVX)),$(strip $(FOVY)) \
	--geom $(strip $(GEOMX)),$(strip $(GEOMY)) \
	--center $(strip $(CENTERRA)),$(strip $(CENTERDEC))

PYTHON=python
TOOLSDIR=$(HOME)/Source/PyFITSTools
MAKEMAP=$(PYTHON) $(TOOLSDIR)/make-fits-image.py $(MAPARGS)
ACCEPTANCE=$(PYTHON) $(TOOLSDIR)/acceptance.py
SUMMER=$(TOOLSDIR)/sum_maps.pl
FLATLIST=$(PYTHON) $(TOOLSDIR)/make-flat-eventlist.py --oversample 1 


.SECONDARY: # clear secondary rule, so intermediate files aren't deleted
	echo "Secondary"

all:  fov_excess.fits 
	@echo "Done processing runs"

# =========================================================================
# Rules to analyze each run
# =========================================================================

# Gamma-hadron separated eventlist
%_event_selected.fits: $(SOURCEDIR)/%_eventlist.fits.gz
	@echo EVENT SELECTION $*
	@ftselect $< $@ $(CUTS) clobber=true $(REDIRECT)

# masked eventlist
%_event_excluded.fits: %_event_selected.fits excluded.reg
	@echo EVENT EXCLUSION $*
	@ftselect $< $@ $(EXCLMASK) clobber=true $(REDIRECT)

# countmap 
%_cmap.fits: %_event_selected.fits
	@echo COUNT MAP $*
	@$(MAKEMAP) --output $@ $< $(REDIRECT)

# excluded count map
%_cmap_excluded.fits: %_event_excluded.fits
	@echo  EXCLUDED COUNT MAP $*
	@$(MAKEMAP) --output $@ $< $(REDIRECT)

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

sum_%: 
	@$(RM) $@
	@$(SUMMER) -v -o $@ $^ 

sum_cmap.fits: $(RUNS_CMAP)
sum_cmap_excluded.fits: $(RUNS_CMAP_EXCL)
sum_accmap.fits: $(RUNS_ACCMAP)

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
	@ftselect $< $@ $(EXCLMASK) clobber=true $(REDIRECT)	


exclmap.fits: flatlist_excluded.fits
	@echo "EXCLUSION MAP: $@"
	@$(MAKEMAP) --output $@ $< $(REDIRECT)

flatmap.fits: flatlist.fits
	@echo "FLAT MAP: $@"
	@$(MAKEMAP) --output $@ $< $(REDIRECT)


# Field-of-view background model map (background is assumed to be the
# scaled acceptance map)
fov_excess.fits: sum_accmap.fits sum_cmap.fits
	@echo "FOV EXCESS MAP: $@"
	@ftpixcalc $@ 'A-B' a=sum_cmap.fits b=sum_accmap.fits \
		clobber=true $(REDIRECT)

fov_excess_excluded.fits: sum_accmap.fits sum_cmap.fits exclmap.fits
	@echo "FOV EXCESS MAP EXCLUDED: $@"
	@ftpixcalc $@ '(A*C)-(B*C)' a=sum_cmap.fits b=sum_accmap.fits \
		c=exclmap.fits clobber=true $(REDIRECT)

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