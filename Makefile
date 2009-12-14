# =========================================================================
# ANALYSIS PARAMETERS:
# =========================================================================
CUTS='(HIL_MSW>-2.0&&HIL_MSW<0.7)&&(HIL_MSL>-2.0&&HIL_MSL<2.0)'
SOURCEDIR=$(HOME)/Analysis/FITSEventLists/HESS_Crab
FOVX=7.0		   # Field of view of output map (X degrees)
FOVY=7.0		   # Field of view of output map (Y degrees)
GEOMX=301		   # number of X bins in map (integer)
GEOMY=301		   # number of Y bins in map (integer)
CENTERRA=83.633333	   # center of output map in RA
CENTERDEC=22.014444	   # center of output map in Dec
EXCLUSIONFILE=excluded.reg # exclusion region file in ascii region format

# =========================================================================
# specify the runlist here: (or use wildcard)
#
# can use for example:
#    EVENTLISTS += $(SOURCEDIR)/run_012345_eventlist.fits
#    EVENTLISTS += $(SOURCEDIR)/run_012346_eventlist.fits
#    EVENTLISTS += $(SOURCEDIR)/run_012347_eventlist.fits
# or 
#    EVENTLISTS=$(wildcard $(SOURCEDIR)/*.fits.gz)
# =========================================================================

EVENTLISTS=$(wildcard $(SOURCEDIR)/*.fits.gz)

include Analysis.mk