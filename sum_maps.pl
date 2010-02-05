#!/usr/bin/perl
# sum maps given on commandline
# summaps.pl -o <sumfile> file1 [file2 ...]

use Getopt::Std;
getopts('o:v');

$output = "sum.fits";
$output = $opt_o if $opt_o;

print "SUMMING: $output \n" if $opt_v;

$count=0;
$tot= scalar(@ARGV);
foreach $image (@ARGV) {
    if (-e $output) {
	$count++;
	$pct = $count/$tot*100;
	printf("[%3d%%] ADDING $image to $output ...\n",$pct) if $opt_v;
	$cmd = "ftpixcalc temp_$output 'A+B' a=$output b=$image clobber=true";
	print `$cmd`; #	    or die "Failed to run: $cmd\n";
	rename( "temp_$output", $output);
    }
    else {
	system("cp -f $image $output");
	print "CREATED $output from $image\n" if $opt_v;
    }
}

unlink "temp_$output"; 
