

RUNS_SIMSELECTED=$(addsuffix _event_reco_selected.fits,$(BASERUNS))
RUNS_RAW=$(addsuffix _raweventlist.fits,$(BASERUNS))


%_event_reco_selected.fits: %_eventlist_reco.fits %_event_verify.txt
	@echo EVENT SELECTION $*
	@ftselect $< $@ $(CUTS) clobber=true $(REDIRECT)

# TODO: make event_selected depend on _eventlist (in the work
# directory). Also make a rule that generates _eventlist.fits.gz from
# $SOURCEDIR/%_eventlist.fits.gz by either: symlinking the file (if no
# columns need to be reconstructed), or by running col-from-lookups if
# they do.  Also, should probably change the name of runlists that
# don't have recon info in them

%_eventlist_reco.fits: $(SOURCEDIR)/%_eventlist.fits.gz
	@echo "ENERGY RECONSTRUCTION $*"
	@$(PYTHON) $(TOOLSDIR)/col-from-lookups.py --type energy $<

spec_arf.fits: $(RUNS_SIMSELECTED)
	@echo "Generating Effective Area file from runs $*..."
	@$(PYTHON) $(TOOLSDIR)/generate-response-matrix.py $(RUNS_SIMSELECTED)
