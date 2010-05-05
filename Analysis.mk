# ------------------------------------------------
# Cherenkov Analysis Makefile script
# Author: Karl Kosack <karl.kosack@cea.fr>
# ------------------------------------------------
#
# This Analysis.mk file should be included from a Makefile that sets
# some or all of the variables in the "Set Defaults" section as well
# as a EVENTLISTS variable that lists what eventlists to process. 
#
# A set of rules for tranforming inputs to outputs are defined here
# (using FTOOLS and various specific tools). These are used to
# generate high-level analysis outputs (like excess or significance
# maps).
#
# Please set the TOOLSDIR variable to where the cherenkov analysis
# tools are located
#
# Note that you should run 'make -j2' or higher if you have a
# multi-core/processor machine, as it will speed up the analysis
# greatly since multiple tools can run simultaneously (the
# dependencies are automatically taken into account)
#
# Also note that you can stop the analysis at any time with CTRL-C,
# and restart it again later and it will start from where it left
# off. There is no need to start from scratch unless you change
# parameters in the Makefile.
#
# For detailed help and instructions, see the "instructions.org" file
# (which is in emacs org-mode format, and can be exported to latex,
# html, ASCII, etc using the 'C-c C-e' command in emacs)
#

# =========================================================================
# ANALYSIS PARAMETERS AND THEIR DEFAULTS: the following parameters are
# user-definable and may be set in the users's analysis Makefile that
# includes Analysis.mk. If they are not set, the default values listed
# below will be used.
# =========================================================================

TOOLSDIR ?= $(HOME)/Source/PyFITSTools
CUTS ?= '(HIL_MSW>-2.0&&HIL_MSW<0.7)&&(HIL_MSL>-2.0&&HIL_MSL<2.0)&&MULTIP>=2'
SOURCEDIR ?=$(HOME)/Analysis/FITSEventLists/HESS_Crab
FOVX ?= 7.0                   # Field of view of output map (X degrees)
FOVY ?= 7.0                   # Field of view of output map (Y degrees)
GEOMX ?= 301                  # number of X bins in map (integer)
GEOMY ?= 301		      # number of Y bins in map (integer)
CENTERRA ?= 83.633333	      # center of output map in RA
CENTERDEC ?= 22.014444	      # center of output map in Dec
EXCLUSIONFILE ?= excluded.reg # exclusion region file in ascii region format
ONRADIUS ?= 0.1               # on-region theta^2
SMOOTHRAD ?= 0.1              # smoothing radius in degrees (for gauss smoothing)
RINGAREAFACTOR ?= 7.0
RINGGAP ?= 0.2                # gap between ring and ON radius in degrees
PROJECTION ?= CAR             # WCS projection type for the map 
MAXEVENTRADIUS ?= 3.0         # maximum event radius in degrees (psi cut)

TOOLSDIR ?= $(HOME)/Source/PyFITSTools
PYTHON ?= python
PYTHONPATH+=:$TOOLSDIR

TARGETS = diagnostic_significance.ps 
TARGETS += fov_excess_gauss.fits fov_significance.fits
TARGETS += ring_excess_gauss.fits

# =========================================================================
# CALCULATED PARAMETERS (not user-definable)
# =========================================================================
GAUSSSIG=$(shell perl -e 'print ($(SMOOTHRAD)*$(GEOMX)/($(FOVX)+1e-12));')
ONTH2=$(shell perl -e 'print $(ONRADIUS)**2;')

# =========================================================================
# format the runlist properly from the user-specified EVENTLISTS list
# =========================================================================
BASERUNS=$(sort $(patsubst %_eventlist.fits,%,$(basename $(notdir $(EVENTLISTS)))))

# =========================================================================
# Sanity checks:
# =========================================================================

ifeq ($(words $(BASERUNS)),0)
$(error RUNLIST IS EMPTY! Check your EVENTLISTS variable and make sure you have not mistyped a directory or something)
endif

# =========================================================================
# utility parameters
# =========================================================================
REDIRECT= >> output.log 2>&1 # set to blank to send output to stdout

# spatial exclusion filter
EXCLMASK='regfilter("$(strip $(EXCLUSIONFILE))",RA,DEC)' 

# parameters passed to MAKEMAP to generate the images
MAPARGS=--fov $(strip $(FOVX)),$(strip $(FOVY)) \
	--geom $(strip $(GEOMX)),$(strip $(GEOMY)) \
	--center $(strip $(CENTERRA)),$(strip $(CENTERDEC)) \
	--proj $(strip $(PROJECTION))

MAKEMAP:=$(PYTHON) $(TOOLSDIR)/make-fits-image.py $(MAPARGS)
ACCEPTANCE:=$(PYTHON) $(TOOLSDIR)/acceptance.py --rmax=$(strip $(MAXEVENTRADIUS))
SUMMER:=$(PYTHON) $(TOOLSDIR)/sum-maps.py
FLATLIST:=$(PYTHON) $(TOOLSDIR)/make-flat-eventlist.py --oversample 1 
MAKERING:=$(PYTHON) $(TOOLSDIR)/make-ring.py
CONVOLVE:=$(PYTHON) $(TOOLSDIR)/convolve-images.py
FOVMASK:=$(PYTHON) $(TOOLSDIR)/make-radial-cutmask.py 
ifndef NOVERIFY
VERIFY:=$(PYTHON) $(TOOLSDIR)/verify-eventlist.py
else
VERIFY:=echo
endif 

.SECONDARY: # clear secondary rule, so intermediate files aren't deleted

.PHONY: setup clean help all verify clean clean-runs clean-sums clean-bg clean-excl clean-events clean-some clean-maps show

all::  setup show $(TARGETS)
	@echo "Done processing: $(TARGETS)"

deps.ps: Makefile
	$(TOOLSDIR)/dependencies.pl | dot -Tps -o $@

deps.pdf:
	$(TOOLSDIR)/dependencies.pl | dot -Tpdf -o $@

setup:
	@echo FTOOLS
	@fversion
	@echo PYTHON '$(PYTHON)'
	@$(PYTHON) -V

help:
	@echo "=============================================================="
	@echo "COMMON TARGETS:            FUNCTION:"
	@echo "    make show              - show the options for this analysis"
	@echo "    make verify            - check that all event-lists are ok"
	@echo "    make fov_excess.fits   - generate FOV model excess map"
	@echo "    make ring_excess.fits  - generate Ring model excess map"
	@echo "    make diagnostic_significance.ps - significance curves"
	@echo ""
	@echo "Note: if you have multiple processors, use 'make -jN' where N "
	@echo "      is the number of simultaneous processes to start."

show:
	@echo "=============================================================="
	@echo  ANALYSIS PARAMETERS:
	@echo "--------------------------------------------------------------"
	@echo "    TOOLS: $(TOOLSDIR)"
	@echo "   SOURCE: $(SOURCEDIR)"
	@echo "     RUNS: $(words $(BASERUNS))"
	@echo "      FOV: $(strip $(FOVX)) x $(strip $(FOVY)) deg"
	@echo "     GEOM: $(strip $(GEOMX)) x $(strip $(GEOMY)) pix"
	@echo "     PROJ: $(strip $(PROJECTION))"
	@echo "       RA: $(strip $(CENTERRA)) deg"
	@echo "      DEC: $(strip $(CENTERDEC)) deg"
	@echo "    R_max: $(strip $(MAXEVENTRADIUS)) deg"
	@echo "    ONRAD: $(strip $(ONRADIUS)) deg ($(strip $(ONTH2)) sq deg)"
	@echo "   SMOOTH: $(strip $(SMOOTHRAD)) deg ($(GAUSSSIG) pix)"
	@echo "=============================================================="

# =========================================================================
# exclusions
# =========================================================================

excluded.reg: $(HESSROOT)/hdanalysis/lists/ExcludedRegions_v11.dat
	$(PYTHON) $(TOOLSDIR)/hess2regfilter.py --scale 1.1 \
		--type fitsex $< > 'regtmp.reg'
	mv regtmp.reg $@

# =========================================================================
# Rules to analyze each run
# =========================================================================


# Gamma-hadron separated eventlist
%_event_selected.fits: $(SOURCEDIR)/%_eventlist.fits.gz %_event_verify.txt
	@echo EVENT SELECTION $*
	@ftselect $< $@ $(CUTS) clobber=true $(REDIRECT)


# masked eventlist
%_event_excluded.fits: %_event_selected.fits $(EXCLUSIONFILE)
	@echo EVENT EXCLUSION $*
	@ftselect $< $@ $(EXCLMASK) clobber=true $(REDIRECT)

# countmap 
%_cmap.fits: %_event_selected.fits 
	@echo COUNT MAP $*
	@$(MAKEMAP) --rmax $(strip $(MAXEVENTRADIUS)) --output $@ $< $(REDIRECT)

# excluded count map
%_cmap_excluded.fits: %_event_excluded.fits
	@echo  EXCLUDED COUNT MAP $*
	@$(MAKEMAP) --rmax $(strip $(MAXEVENTRADIUS)) --output $@ $< $(REDIRECT)

# acceptance map
%_accmap.fits: %_event_excluded.fits %_cmap.fits flatlist_excluded.fits
	@echo ACCEPTANCE MAP $*
	@$(ACCEPTANCE) --exflat flatlist_excluded.fits --output $@ \
		$*_event_excluded.fits $*_cmap.fits $(REDIRECT)

# ==============================================================
# Summed maps
# ==============================================================

RUNS_CMAP=$(addsuffix _cmap.fits,$(BASERUNS))
RUNS_OFFMAP_RING=$(addsuffix _offmap_ring.fits,$(BASERUNS))
RUNS_CMAP_EXCL=$(addsuffix _cmap_excluded.fits,$(BASERUNS))
RUNS_ACCMAP=$(addsuffix _accmap.fits,$(BASERUNS))
RUNS_ACCMAP_RING=$(addsuffix _accmap_ring.fits,$(BASERUNS))

%_sum.fits: 
	@echo "CREATING SUM: $@"
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

# a dummy eventlist that is flat across the image (each pixel is
# sampled N times in a grid, where N is set in the FLATLIST
# macro. This is useful for generating a flat map or for calculating
# the effect of exclusion regions (anywhere where you would want to do
# a loop over pixel positions)
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


# mask large radii from field of view (the psi^2 cut)
%_fovmask.fits: %_cmap.fits
	@$(FOVMASK) $(strip $(MAXEVENTRADIUS)) $< $@


# tophat correlation of a map
%_tophat.fits: %.fits tophat.fits
	@echo "TOPHAT CONVOLUTION: $*"
	@$(CONVOLVE) --output $@ $< tophat.fits $(REDIRECT)

# gaussian correlation of a map:
%_gauss.fits: %.fits 
	@echo "GAUSS SMOOTH: $* :  $(GAUSSSIG) pix"
	@fgauss $< $@ $(GAUSSSIG) ratio=1.0 theta=0.0 nsigma=4.0 \
		boundary=nearest $(REDIRECT)


%_excess_acorr.fits: %_excess.fits accmap_sum.fits
	@echo "ACCEPTANCE CORRECTED EXCESS: $@"
	@ftimgcalc $@ 'A*max(B)/B' a=$< b=accmap_sum.fits


tophat.fits: exclmap.fits
	@$(MAKERING) --output $@ --onradius=$(strip $(ONRADIUS)) \
		    --make-on $^ $(REDIRECT)	

%_exmasked.fits: %.fits exclmap.fits
	@echo "EXCLUSION MASKING: $<"
	@ftpixcalc $@ 'A*B' a=$< b=exclmap.fits clobber=yes $(REDIRECT)

%_fovmasked.fits: %.fits %_fovmask.fits
	@echo "FOV MASKING: $<"	
	@ftpixcalc $@ 'A*B' a=$< b=$*_fovmask.fits clobber=yes $(REDIRECT)

#-------------------------------------------------------------------------
# Field-of-view background model maps (background is assumed to be the
# scaled acceptance map)
#-------------------------------------------------------------------------
fov_excess.fits: accmap_sum.fits cmap_sum.fits
	@echo "FOV EXCESS MAP: $@"
	@ftpixcalc $@ 'A-B' a=cmap_sum.fits b=accmap_sum.fits \
		clobber=true $(REDIRECT)

fov_significance.fits: cmap_sum_tophat.fits accmap_sum_tophat.fits 
	@echo "FOV SIGNIFICANCE: $@"
	@ftpixcalc $@ '(A-B)/sqrt(A)'\
		a=cmap_sum_tophat.fits \
		b=accmap_sum_tophat.fits \
		clobber=yes $(REDIRECT)

#-------------------------------------------------------------------------
# ring background model (off events are sum of on events in annulus about
# each bin)
#-------------------------------------------------------------------------

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

%_offmap_ring.fits: %_cmap_excluded.fits ring.fits exclmap_ring.fits
	@echo "RING OFF MAP: $*"
	@$(CONVOLVE) --output tmp_$@ $^ $(REDIRECT)
	@ftimgcalc $@ 'A*sum(C)/B' a=tmp_$@ b=exclmap_ring.fits c=ring.fits clobber=true $(REDIRECT)
	@$(RM) tmp_$@


%_accmap_ring.fits: %_accmap.fits ring.fits flatmap_ring.fits
	@echo "CONVOLVE RING ACCMAP: $*"
	@$(CONVOLVE) --output $@ $^ $(REDIRECT)
#	@$(CONVOLVE) --output tmp_$@ $^ $(REDIRECT)
#	@ftpixcalc $@ 'A/B' a=tmp_$@ b=flatmap_ring.fits clobber=true $(REDIRECT)
#	@$(RM) tmp_$@


ring_alpha.fits: accmap_sum.fits accmap_ring_sum.fits
	@echo "RING ALPHA MAP"
	@ftpixcalc $@ 'A/B' a=accmap_sum.fits b=accmap_ring_sum.fits \
		clobber=true $(REDIRECT)


ring_excess.fits: cmap_sum.fits ring_alpha.fits offmap_ring_sum.fits
	@echo "RING EXCESS: $@"
	@ftpixcalc $@ 'A-B*C' \
		a=cmap_sum.fits \
		b=offmap_ring_sum.fits \
		c=ring_alpha.fits \
		clobber=yes $(REDIRECT)


ring_alpha_tophat.fits: accmap_sum_tophat.fits accmap_ring_sum_tophat.fits
	@echo "RING ALPHA MAP TOPHAT: $@"
	@ftpixcalc $@ 'A/B' a=accmap_sum_tophat.fits b=accmap_ring_sum_tophat.fits \
		clobber=true $(REDIRECT)


LIMA='sqrt(2)*sqrt( A*log( (1+C)/C * (A/(A+B))) + B*log((1+C)*B/(A+B)) )'
ring_significance.fits: cmap_sum_tophat.fits ring_alpha_tophat.fits offmap_ring_sum_tophat.fits
	@echo "RING SIGNIFICANCE: $@"
	@ftpixcalc alt_$@ $(LIMA)\
		a=cmap_sum_tophat.fits \
		b=offmap_ring_sum_tophat.fits \
		c=ring_alpha_tophat.fits \
		clobber=yes $(REDIRECT)
	@ftpixcalc $@ '(A-B*C)/sqrt(C*(A+B))'\
		a=cmap_sum_tophat.fits \
		b=offmap_ring_sum_tophat.fits \
		c=ring_alpha_tophat.fits \
		clobber=yes $(REDIRECT)

# TODO: need to divide by exclmap_tophat to get the right values in the exclusion region!
ring_significance_exmasked.fits: cmap_sum_tophat_exmasked.fits ring_alpha_exmasked_tophat.fits offmap_ring_sum_exmasked_tophat.fits
	@echo "RING SIGNIFICANCE: $@"
	@ftpixcalc alt_$@ $(LIMA)\
		a=cmap_sum_tophat_exmasked.fits \
		b=offmap_ring_sum_exmasked_tophat.fits \
		c=ring_alpha_exmasked_tophat.fits \
		clobber=yes $(REDIRECT)
	@ftpixcalc $@ '(A-B*C)/sqrt(C*(A+B))'\
		a=cmap_sum_tophat_exmasked.fits \
		b=offmap_ring_sum_tophat.fits \
		c=ring_alpha_tophat.fits \
		clobber=yes $(REDIRECT)


#-------------------------------------------------------------------------
# Diagnostic plots (some use gnuplot)
#-------------------------------------------------------------------------

diagnostic_significance.ps: $(TOOLSDIR)/diagnostic_significance.gpl ring_significance_imhist.fits ring_significance_exmasked_imhist.fits
	@echo "DIAGNOSTIC: $@"
	gnuplot $< > $@


%_imhist.fits: %.fits
	@echo "HISTOGRAM: $@"
	@fimhisto $< $@ range=-10,100 binsize=0.5 clobber=yes $(REDIRECT)

# TODO: make a rule like:
#  %_significance.fits: cmap_sum_tophat.fits %_offmap_sum_tophat.fits
# (will need to rename some things like ring_offmap_sum_tophat.fits



# ==============================================================
# Utilities
# ==============================================================

%_event_verify.txt: $(SOURCEDIR)/%_eventlist.fits.gz
	@echo "VERIFY $*"
	@$(VERIFY) $< $@

verify: $(addsuffix _event_verify.txt, $(BASERUNS))
	@echo "VERIFY: Your runlist passes verification"


# ==============================================================
# cleanup
# ==============================================================


clean: clean-some

distclean: clean-events clean-runs clean-sums clean-bg clean-excl 
	$(RM) output.log
	$(RM) *_tophat.fits
	$(RM) deps.ps
	$(RM) run*_event_verify.txt
	$(RM) diagnostic_*.ps

clean-some: clean-runs clean-sums clean-excl


clean-runs:
	$(RM) $(RUNS_CMAP) $(CLEANUP_RUNS)
	$(RM) $(RUNS_CMAP_EXCL)
	$(RM) $(RUNS_ACCMAP)
	$(RM) $(RUNS_ACCMAP_RING)
	$(RM) $(RUNS_OFFMAP_RING)

clean-events:
	$(RM) *_event_excluded.fits $(CLEANUP_EVENTS)
	$(RM) *_event_selected.fits

clean-sums:
	$(RM) *_sum.fits *_sum_*.fits

clean-bg:
	$(RM) fov_*.fits $(CLEANUP_BG)
	$(RM) ring*.fits *_ring*.fits
	$(RM) tophat.fits

clean-excl:
	$(RM) exclmap.fits flatmap.fits flatlist.fits flatlist_excluded.fits

# useful if you change map parameters (projection, center, etc) and
# want to redo everything
clean-maps: clean-bg clean-some
