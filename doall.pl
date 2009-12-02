#!/usr/bin/perl

use File::Basename;


$Dryrun=0;
$CenterRA = 83.633333;
$CenterDec= 22.014444;
$GeomX = 300;
$GeomY = 300;
$FOVX = 10.0;
$FOVY = 10.0;

$outdir = "/home/kkosack/Analysis/FITSEventLists/Analysis";

$isfirst=1;

foreach $evlist (<~/Analysis/FITSEventLists/HESS_Crab/*.fits.gz>) {
    
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
    $mask = "regfilter(\"excluded.reg\",RA,DEC)";
    if (! -e $maskname) {
	print "EXCLUDE: --> $maskname \n";
	runtool( "ftselect $selname $maskname '$mask'")
    }
   
    # make a countmap

    $cmapname="${bname}_cmap.fits";
    $maskcmapname="${bname}_cmap_masked.fits";

    makeMap( $cmapname, $selname );
    makeMap( $maskcmapname, $maskname );

    # update the summed countmap
    
    updateSum( "$outdir/cmap_sum.fits", $cmapname );
    updateSum( "$outdir/cmap_sum_masked.fits", $maskcmapname );

    $isfirst=0;

}

sub makeMap($$) {
    my $mapname=shift;
    my $eventlist=shift;
    my $args="--fov $FOVX,$FOVY --geom $GeomX,$GeomY --center $CenterRA,$CenterDec";
     
    if (! -e $mapname ) {
	print "CMAP: --> $mapname\n";
	runtool( "makefits.py $args --output $mapname $eventlist" );
    }

}

sub updateSum($$){

    my $sumname = shift();
    my $input = shift();
    my $tmp = "$outdir/temp_sum.fits";

    print "SUM: $sumname\n";

    if ($isfirst) {
	system( "cp $input $sumname");
	unlink $tmp if (-e $tmp);
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
