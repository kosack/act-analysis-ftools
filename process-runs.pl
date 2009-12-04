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

unlink "$outdir/cmap_sum.fits";
unlink "$outdir/cmap_sum_masked.fits";

$PYTHON="python2.6";
if (system("python2.6 -V")){
    $PYTHON="python2.4";
}

$is_first_iter = 1;

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

    # update the summed countmap
    
    updateSum( "$outdir/cmap_sum.fits", $cmapname );
    updateSum( "$outdir/cmap_sum_masked.fits", $maskcmapname );
    updateSum( "$outdir/acc_sum.fits", $accname );

    $is_first_iter=0

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

    print "SUM: $sumname\n";

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
	system($command) ;#or die "Command failed: '$command'\n";
    }

}
