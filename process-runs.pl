#!/usr/bin/perl

use File::Basename;
use strict 'vars';

my $Dryrun=0;
my $CenterRA = 83.633333;
my $CenterDec= 22.014444;
my $GeomX = 301;
my $GeomY = 301;
my $FOVX = 7.0;
my $FOVY = 7.0;

my $indir= "$ENV{HOME}/Analysis/FITSEventLists/HESS_Crab";
my $Outdir = "$ENV{HOME}/Analysis/FITSEventLists/Analysis";
my @infilelist = <$indir/*.fits.gz>; # list of input eventlist files




my $PYTHON="python2.6";
if (system("python2.6 -V")){
    $PYTHON="python2.4";
}

my %Sumlist;
my $iteration = 0;
my $flatlistname ="$Outdir/flatlist.fits";
my $excllistname ="$Outdir/excllist.fits";
my $flatmapname = "$Outdir/flatmap.fits";
my $exclmapname = "$Outdir/exclmap.fits";


MAIN: {

    foreach my $evlist (@infilelist) {

	$iteration++;

	print "="x70,"\n";
	print " PROCESSING $iteration of ",scalar(@infilelist)," ...\n";
	print " $evlist \n";
	print "="x70,"\n";
	
	processEventList( $evlist, $iteration );

    } 


    print "/\\"x40,"\n";

    # Now do all the sums

    my %sum;
    foreach my $sumname (keys %Sumlist) {
	$sum{$sumname} = "$Outdir/sum_$sumname.fits";

	if (!-e $sum{$sumname}) {
	    print "SUMMING: $sumname ...\n";
	    foreach my $map (@{$Sumlist{$sumname}}) {
		print "\t$map\n";
		
		updateSum( $sum{$sumname}, $map );
	    }
	}
    }


    # make fov excess map from sums:
    my $fovexcessname = "$Outdir/excess_fov.fits";
    print "FOV EXCESS MAP: $fovexcessname \n";
    runtool("ftpixcalc $fovexcessname 'A-B' a=$sum{cmap} b=$sum{acc} clobber=true");


    # make a significance map from sums: (for FOV model, the error is
    # the error on N_off, and alpha=1, so the formula is very simple:)
    # TODO: do the smoothing with tophat not boxcar
    my $fovsigname = "$Outdir/significance_fov.fits";
    my $fovsignamemask = "$Outdir/significance_excl_fov.fits";
    my $fovexcessname_cor = "$Outdir/excess_fov_cor.fits";
    my $acc_cor = "$Outdir/sum_acc_cor.fits";
    runtool("fboxcar $fovexcessname $fovexcessname_cor 10 10 clobber=yes");
    runtool("fboxcar $sum{acc} $acc_cor 10 10 clobber=yes");
    print "FOV SIGNIF MAP: $fovsigname \n";
    runtool( "ftpixcalc $fovsigname 'A/(sqrt(B))' a=$fovexcessname_cor b=$acc_cor clobber=true" );

}

sub makeMap($$) {
    my $mapname=shift;
    my $eventlist=shift;
    my $args="--fov $FOVX,$FOVY --geom $GeomX,$GeomY ".
	"--center $CenterRA,$CenterDec";
     
    if (! -e $mapname ) {
	print "CMAP: --> $mapname\n";
	runtool( "$PYTHON makefits.py $args --output $mapname $eventlist" );
    }

}

sub updateSum($$){

    my $sumname = shift();
    my $input = shift();
    my $tmp = "$Outdir/temp_sum.fits";

    #print "SUM: $sumname <- $input\n";

    if (! -e $sumname) {
	system( "cp $input $sumname");
	unlink $tmp if (-e $tmp);
	print "\tCREATED $sumname on first iteration\n";
    }
    else {
	runtool("ftpixcalc $tmp 'A+B' a=$sumname b=$input clobber=true");
	runtool("rm -f $sumname && mv $tmp $sumname");
    }


}


sub runtool($) {
    my $command = shift();
    
    if ($Dryrun) {
	print "\% $command\n";
    }
    else {
	system($command);
	if ($? == -1) {
	    die "failed to execute: $!\n";
	}
	elsif ($? & 127) {
	    die "!"x70,"\nKilled at command: \n $command\n";
	}
	else {
	    sprintf("child exited with value %d\n", $? >> 8);
	}

    }

}

#
# Do all necessary processing on a single eventlist
#
#
sub processEventList($$) {

    my $evlist = shift;
    my $iteration = shift;

    my $run=0;
    $_ = $evlist; 
    $run = $1+0 if (/.*run_([0-9]+).*/);

    my $subdir = sprintf("run%06d", $run );

    my $bname=basename($evlist,"_eventlist.fits.gz");
    $bname="$Outdir/$subdir/$bname";
    mkdir("$Outdir/$subdir") if (! -e "$Outdir/$subdir");

    # apply gamma-hadron cuts
    my $selname = "${bname}_event_selected.fits.gz";
    my $inname="${evlist}[EVENTS]";
    my $cuts = '(HIL_MSW > -2.0 && HIL_MSW < 0.7) ' 
	. '&& (HIL_MSL > -2.0 && HIL_MSL < 2.0)';

    if (! -e $selname ) {
	print "SELECT: --> $selname\n";
	runtool( "ftselect $inname $selname '$cuts'");
    }

    # apply spatial cuts, masking out known sources, and generate
    # a "masked" eventlist (useful for acceptance calculation)

    my $maskname = "${bname}_event_masked.fits.gz";
    my $exclmask = "regfilter(\"excluded.reg\",RA,DEC)";
    if (! -e $maskname) {
	print "EXCLUDE: --> $maskname \n";
	runtool( "ftselect $selname $maskname '$exclmask'")
    }

    # generate a 1-D reflected-region background selector:
    # my $refregionname = "${bname}_OFF.reg";
    # if (!-e $refregionname) {
    # 	my $oldcwd;
    # 	chomp($oldcwd = `pwd`);
    # 	chdir("$Outdir/$subdir");
    # 	runtool("python $oldcwd/regionbg.py 0.01 $evlist");
    # 	chdir($oldcwd);
    # }

    # make a countmap

    my $cmapname="${bname}_cmap.fits";
    my $maskcmapname="${bname}_cmap_masked.fits";

    makeMap( $cmapname, $selname );
    makeMap( $maskcmapname, $maskname );

    # generate exclusion map, which is just a binned version of
    # the exclusion mask. This can be used for example to
    # calculate the area covered by masked regions

    if ($iteration==1) {
	unlink $flatlistname;
	unlink $excllistname;

	print "FLAT EVENTLISTS: $flatlistname \n";
	# make flat eventlist
	runtool("$PYTHON make-flat-eventlist.py -s 1 $cmapname $flatlistname");
	# make excluded flat eventlist
	runtool("ftselect $flatlistname $excllistname '$exclmask'");
	# make exclusion map and flat map
	print "EXCLUSION MAP: $exclmapname \n";
	makeMap( $flatmapname, $flatlistname );
	makeMap( $exclmapname, $excllistname );
    }

    # generate acceptance map from the data:

    my $accname = "${bname}_acc.fits";
    if (! -e $accname ){ 
	print "ACCEPTANCE MAP: $accname \n";
	runtool( "$PYTHON acceptance.py --output $accname ".
		 "$maskname $cmapname" );
    }

    # generate runwise FOV-model excess map (uses the acceptance
    # map as the background map)
    my $excessname = "${bname}_fov_excess.fits";
    my $excessnamemasked = "${bname}_fov_excess_masked.fits";
    if (! -e $excessname  ) {
	print "FOV EXCESS: $excessname\n";
	runtool( "ftpixcalc $excessname 'A-B' a=$cmapname b=$accname" );
	runtool( "ftpixcalc $excessnamemasked 'A-B' ".
		 "a=$maskcmapname b=$accname" );
    }

    # update the maps that should be summed over runs:

    push @{$Sumlist{cmap}}, $cmapname;
    push @{$Sumlist{cmap_masked}}, $maskcmapname;
    push @{$Sumlist{acc}}, $accname;
    push @{$Sumlist{fov_excess}}, $excessname;

#	updateSum( $sum{cmap}, $cmapname );
#	updateSum( $sum{cmap_masked}, $maskcmapname );
#	updateSum( $sum{acc}, $accname );
#	updateSum( $sum{fov_excess}, $excessname );

    # debug:
#	print "DEBUG: FINISHING ON FIRST RUN !!!!!! \n";
#	last; 


}
