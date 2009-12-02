#!/usr/bin/perl

use File::Basename;


$Dryrun=0;
$outdir = "/home/kkosack/Analysis/FITSEventLists/Analysis";

$isfirst=1;

foreach $evlist (<~/Analysis/FITSEventLists/HESS_Crab/*.fits.gz>) {
    
    $bname=basename($evlist,"_eventlist.fits.gz");
    $bname="$outdir/$bname";

    print " =================================================\n";
    print " $evlist \n";
    print " $bname \n";
    print " =================================================\n";

    # apply gamma-hadron cuts
    $selname = "${bname}_event_selected.fits";
    $inname="${evlist}[EVENTS]";
    $cuts = 'HIL_MSW > -2.0 && HIL_MSW < 0.7 && HIL_MSL > -2.0 && HIL_MSL < 2.0';
    print "SELECT: --> $selname\n";
    
    if (! -e $selname ) {
	runtool( "ftselect $inname $selname '$cuts'");
    }
    
    # apply sourcemask cuts (not sure this is the best way to do this,
    # but it works)

    $maskname = "${bname}_event_masked.fits";
    $mask = "regfilter(\"excluded.reg\",RA,DEC)";
    if (! -e $maskname) {
	runtool( "ftselect $selname $maskname '$mask'")
    }
   
    # make a countmap

    $cmapname="${bname}_cmap.fits";
    $maskcmapname="${bname}_cmap_masked.fits";
    $args = "--fov 10,10 --geom 300,300 --center 83.633333,22.014444";
    print "CMAP: --> $cmapname\n";

    if (! -e $cmapname ) {
	runtool( "makefits.py $args --output $cmapname $selname" );
    }
    if (! -e $maskcmapname ) {
	runtool( "makefits.py $args --output $maskcmapname $maskname" );
    }

    die "stop here for now\n";



    # update the summed countmap
    
    $sumname = "$outdir/cmap_sum.fits";
    $tmp = "$outdir/temp_sum.fits";

    if ($isfirst) {
	system( "cp $cmapname $sumname");
	unlink $tmp if (-e $tmp);
    }
    else {
	runtool("ftpixcalc $tmp 'A+B' a=$sumname b=$cmapname clobber=true");
	runtool("rm -f $sumname && mv $tmp $sumname");
    }

    
    $isfirst=0;

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
