#!/usr/bin/perl

use File::Basename;


$Dryrun=0;
$CenterRA = 83.633333;
$CenterDec= 22.014444;
$GeomX = 256;
$GeomY = 256;
$FOVX = 6.0;
$FOVY = 6.0;

$outdir = "$ENV{HOME}/Analysis/FITSEventLists/Analysis";

$sum{cmap} = "$outdir/sum_cmap.fits";
$sum{cmap_masked} = "$outdir/sum_cmap_masked.fits";
$sum{acc} = "$outdir/sum_acc.fits";
$sum{fov_excess} = "$outdir/sum_fov_excess.fits";

foreach $map (keys %sum) {
    print "Clearing sum: $map\n";
    unlink $sum{$map};
}

$PYTHON="python2.6";
if (system("python2.6 -V")){
    $PYTHON="python2.4";
}

$is_first_iter = 1;


MAIN: {

    foreach $evlist (<$ENV{HOME}/Analysis/FITSEventLists/HESS_Crab/*.fits.gz>) {

	$bname=basename($evlist,"_eventlist.fits.gz");
	$bname="$outdir/$bname";

	print "="x70,"\n";
	print " $evlist \n";
	print "="x70,"\n";

	# apply gamma-hadron cuts
	$selname = "${bname}_event_selected.fits.gz";
	$inname="${evlist}[EVENTS]";
	$cuts = 'HIL_MSW > -2.0 && HIL_MSW < 0.7 && HIL_MSL > -2.0 && HIL_MSL < 2.0';

	if (! -e $selname ) {
	    print "SELECT: --> $selname\n";
	    runtool( "ftselect $inname $selname '$cuts'");
	}

	# apply sourcemask cuts (not sure this is the best way to do this,
	# but it works)

	$maskname = "${bname}_event_masked.fits.gz";
	$exclmask = "regfilter(\"excluded.reg\",RA,DEC)";
	if (! -e $maskname) {
	    print "EXCLUDE: --> $maskname \n";
	    runtool( "ftselect $selname $maskname '$exclmask'")
	}

	# make a countmap

	$cmapname="${bname}_cmap.fits";
	$maskcmapname="${bname}_cmap_masked.fits";

	makeMap( $cmapname, $selname );
	makeMap( $maskcmapname, $maskname );

	# generate exclusion map
	$flatlistname ="$outdir/flatlist.fits";
	$excllistname ="$outdir/excllist.fits";
	$flatmapname = "$outdir/flatmap.fits";
	$exclmapname = "$outdir/exclmap.fits";
	if ($is_first_iter) {
	    unlink $flatlistname;
	    unlink $excllistname;

	    # make flat eventlist
	    runtool("$PYTHON make-flat-eventlist.py -s 1 $cmapname $flatlistname");
	    # make excluded flat eventlist
	    runtool("ftselect $flatlistname $excllistname '$exclmask'");
	    # make exclusion map and flat map
	    makeMap( $flatmapname, $flatlistname );
	    makeMap( $exclmapname, $excllistname );
	}

	# generate acceptance map

	$accname = "${bname}_acc.fits";
	if (! -e $accname ){ 
	    runtool( "$PYTHON acceptance.py --output $accname $maskname $cmapname" );
	}

	# generate runwise FOV-model excess map
	$excessname = "${bname}_fov_excess.fits";
	$excessnamemasked = "${bname}_fov_excess_masked.fits";
	if (! -e $excessname  ) {
	    print "FOV EXCESS: $excessname\n";
	    runtool( "ftpixcalc $excessname 'A-B' a=$cmapname b=$accname" );
	    runtool( "ftpixcalc $excessnamemasked 'A-B' a=$maskcmapname b=$accname" );
	}

	# update the summed countmap

	updateSum( $sum{cmap}, $cmapname );
	updateSum( $sum{cmap_masked}, $maskcmapname );
	updateSum( $sum{acc}, $accname );
	updateSum( $sum{fov_excess}, $excessname );

	$is_first_iter=0

    }

    # make fov excess map from sums:
    $fovexcessname = "$outdir/excess_fov.fits";
    print "FOV EXCESS MAP: ";
    runtool("ftpixcalc $fovexcessname 'A-B' a=$sum{cmap} b=$sum{acc} clobber=true");

}

sub makeMap($$) {
    my $mapname=shift;
    my $eventlist=shift;
    my $args="--fov $FOVX,$FOVY --geom $GeomX,$GeomY --center $CenterRA,$CenterDec";
     
    if (! -e $mapname ) {
	print "CMAP: --> $mapname\n";
	runtool( "$PYTHON makefits.py $args --output $mapname $eventlist" );
    }

}

sub updateSum($$){

    my $sumname = shift();
    my $input = shift();
    my $tmp = "$outdir/temp_sum.fits";

    print "SUM: $sumname <- $input\n";

    if (! -e $sumname) {
	system( "cp $input $sumname");
	unlink $tmp if (-e $tmp);
	print "CREATED $sumname on first iteration\n";
    }
    else {
	runtool("ftpixcalc $tmp 'A+B' a=$sumname b=$input clobber=true");
	runtool("rm -f $sumname && mv $tmp $sumname");
    }


}


sub runtool($) {
    my $command = shift();
    
    if ($Dryrun) {
	print "$command\n";
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
