#!/usr/bin/perl
# sum maps given on commandline
# summaps.pl -o <sumfile> file1 [file2 ...]

use Getopt::Std;
$opt{o} = "sum.fits";
getopt('o:', \%opt);

print "Output: $opt{o}\n";

$output = $opt{o};

foreach $image (@ARGV) {
    if (-e $output) {
	print "ADDING $image to $output...\n";
	system( "ftpixcalc temp_$output 'A+B' a=$output b=$image clobber=true");
	rename( "temp_$output", $output);
    }
    else {
	rename $image, $output;
	print "CREATED $output from $image\n";
    }
}

