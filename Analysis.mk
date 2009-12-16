# ------------------------------------------------
# Cherenkov Analysis Makefile script
# Author: Karl Kosack <karl.kosack@cea.fr>
# ------------------------------------------------
#
# This Analysis.mk file should be included from a Makefile that sets
# some or all of the variables in the "Set Defaults" section as well
# as a EVENTLISTS variable that lists what eventlists to process. A
# set of rules for tranforming inputs to outputs are defined here (using
# FTOOLS and various specific tools). These are used to generate
# high-level analysis outputs (like excess or significance maps).
#
# Please set the TOOLSDIR variable to where the cherenkov analysis
# tools are located
#
# Note that you should run 'make -j2' or higher if you have a
# multi-core/processor machine, as it will speed up the analysis
# greatly since multiple tools can run simultaneously (the
# dependencies are automatically taken into account)

# some high-level targets are for example:
#   make ring_excess.fits
#   make fov_excess.fits

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
ONRADIUS ?= 0.1                  # on-region theta^2

TOOLSDIR ?= $(HOME)/Source/PyFITSTools
PYTHON ?= python

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



MAKEMAP=$(PYTHON) $(TOOLSDIR)/make-fits-image.py $(MAPARGS)
ACCEPTANCE=$(PYTHON) $(TOOLSDIR)/acceptance.py
SUMMER=$(TOOLSDIR)/sum_maps.pl
FLATLIST=$(PYTHON) $(TOOLSDIR)/make-flat-eventlist.py --oversample 1 
MAKERING=$(PYTHON) $(TOOLSDIR)/make-ring.py
CONVOLVE=$(PYTHON) $(TOOLSDIR)/convolve-images.py

.SECONDARY: # clear secondary rule, so intermediate files aren't deleted
	echo "Secondary"

.PHONY: setup clean help all verify clean clean-runs clean-sums clean-fov clean-excl


all:  setup ring_significance.fits 
	@echo "Done processing runs"


setup:
	@echo FTOOLS
	@fversion
	@echo PYTHON '$(PYTHON)'
	@$(PYTHON) -V

help:
	@echo "=============================================================="
	@echo "COMMON TARGETS:            FUNCTION:"
	@echo "    make verify            - check that all event-lists are ok"
	@echo "    make fov_excess.fits   - generate FOV model excess map"
	@echo "    make ring_excess.fits  - generate Ring model excess map"
# =========================================================================
# Rules to analyze each run
# =========================================================================


# Gamma-hadron separated eventlist
%_event_selected.fits: $(SOURCEDIR)/%_eventlist.fits.gz %_event_verify.txt
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
RUNS_OFFMAP_RING=$(addsuffix _offmap_ring.fits,$(BASERUNS))
RUNS_CMAP_EXCL=$(addsuffix _cmap_excluded.fits,$(BASERUNS))
RUNS_ACCMAP=$(addsuffix _accmap.fits,$(BASERUNS))
RUNS_ACCMAP_RING=$(addsuffix _accmap_ring.fits,$(BASERUNS))

%_sum.fits: 
	@$(RM) $@
	@$(SUMMER) -v -o $@ $^ 

cmap_sum.fits: $(RUNS_CMAP)
cmap_excluded_sum.fits: $(RUNS_CMAP_EXCL)
accmap_sum.fits: $(RUNS_ACCMAP)
offmap_ring_sum.fits: $(RUNS_OFFMAP_RING)
accmap_ring_sum.fits: $(RUNS_ACCMAP_RING)

# ==============================================================
# Background models and final maps
# ==============================================================

#-------------------------------------------------------------------------
# maps common to all Background methods (exclusion map, etc)
#-------------------------------------------------------------------------

FIRSTCMAP=$(firstword $(RUNS_CMAP))

flatlist.fits: $(FIRSTCMAP)
	@echo "FLAT EVENTLIST: $@ using $^"
	@$(FLATLIST) $(FIRSTCMAP) $@ $(REDIRECT)

flatlist_excluded.fits: flatlist.fits $(EXCLUSIONFILE)
	@echo "FLAT EVENTLIST EXCLUDED: $@"
	@ftselect $< $@ $(EXCLMASK) clobber=true $(REDIRECT)	


exclmap.fits: flatlist_excluded.fits
	@echo "EXCLUSION MAP: $@"
	@$(MAKEMAP) --output $@ $< $(REDIRECT)

flatmap.fits: flatlist.fits
	@echo "FLAT MAP: $@"
	@$(MAKEMAP) --output $@ $< $(REDIRECT)

# tophat correlation of a map
%_corr.fits: %.fits tophat.fits
	@echo "TOPHAT CONVOLUTION: $*"
	@$(CONVOLVE) --output $@ $< tophat.fits $(REDIRECT)

tophat.fits: exclmap.fits
	$(MAKERING) --output $@ --onradius=$(strip $(ONRADIUS)) \
		    --make-on $^ $(REDIRECT)	

#-------------------------------------------------------------------------
# Field-of-view background model maps (background is assumed to be the
# scaled acceptance map)
#-------------------------------------------------------------------------
fov_excess.fits: accmap_sum.fits cmap_sum.fits
	@echo "FOV EXCESS MAP: $@"
	@ftpixcalc $@ 'A-B' a=cmap_sum.fits b=accmap_sum.fits \
		clobber=true $(REDIRECT)

fov_excess_excluded.fits: accmap_sum.fits cmap_sum.fits exclmap.fits
	@echo "FOV EXCESS MAP EXCLUDED: $@"
	@ftpixcalc $@ '(A*C)-(B*C)' a=cmap_sum.fits b=accmap_sum.fits \
		c=exclmap.fits clobber=true $(REDIRECT)

#-------------------------------------------------------------------------
# ring background model (off events are sum of on events in annulus about
# each bin)
#-------------------------------------------------------------------------

ring.fits: exclmap.fits
	@echo "GENERATING RING: $@"
	@$(MAKERING) --output $@ --onradius=$(strip $(ONRADIUS)) \
		--areafactor 7.0 $^ $(REDIRECT)


exclmap_ring.fits: exclmap.fits ring.fits
	@echo "CONVOLVE RING EXCLMAP: $*"
	@$(CONVOLVE) --output $@ $^ $(REDIRECT)

%_offmap_ring.fits: %_cmap_excluded.fits ring.fits exclmap_ring.fits
	@echo "RING OFF MAP: $*"
	@$(CONVOLVE) --output tmp_$@ $^ $(REDIRECT)
	@ftpixcalc $@ 'A/B' a=tmp_$@ b=exclmap_ring.fits clobber=true $(REDIRECT)
	@$(RM) tmp_$@

%_accmap_ring.fits: %_accmap.fits ring.fits exclmap_ring.fits
	@echo "CONVOLVE RING ACCMAP: $*"
	@$(CONVOLVE) --output tmp_$@ $^ $(REDIRECT)
	@ftpixcalc $@ 'A/B' a=tmp_$@ b=exclmap_ring.fits clobber=true $(REDIRECT)
	@$(RM) tmp_$@

# TODO: this doesn't take in acceptance right, I think! (it works
# somewhat beacuse the ring is scaled by the area factor, so it's
# basically offmap is really alpha*off)
ring_excess.fits: cmap_sum.fits offmap_ring_sum.fits
	@echo "RING EXCESS"
	@ftpixcalc $@ 'A-B' a=cmap_sum.fits \
		b=offmap_ring_sum.fits clobber=yes $(REDIRECT)


# TODO: need to properly calcualte alpha and use Li/Ma formula! 
ring_significance.fits: cmap_sum_corr.fits offmap_ring_sum_corr.fits
	@echo "RING SIGNIFICANCE: $@"
	@ftpixcalc $@ '(A-B*1.0)/sqrt(1.0*(A+B))' \
		a=cmap_sum_corr.fits \
		b=offmap_ring_sum_corr.fits clobber=yes $(REDIRECT)

# TODO: make a rule like:
#  %_significance.fits: cmap_sum_corr.fits %_offmap_sum_corr.fits
# (will need to rename some things like ring_offmap_sum_corr.fits

# ==============================================================
# Utilities
# ==============================================================

%_event_verify.txt: $(SOURCEDIR)/%_eventlist.fits.gz
	@echo "VERIFY $*"
	@ftverify $< > $@

verify: $(addsuffix _event_verify.txt, $(BASERUNS))
	@echo "VERIFY: Your runlist passes verification"
clean::
	$(RM) run*_event_verify.txt

# ==============================================================
# cleanup
# ==============================================================



clean:: clean-runs clean-sums clean-fov clean-excl
	$(RM) output.log

clean-runs:
	$(RM) run_*_*.fits

clean-sums:
	$(RM) *_sum.fits

clean-fov:
	$(RM) fov_excess.fits	

clean-excl:
	$(RM) exclmap.fits flatmap.fits