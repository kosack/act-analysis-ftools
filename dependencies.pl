#!/usr/bin/perl 

# converts make dependency tree to dot (kind of a hack, but useful for
# debugging the makefile)

$usenames = 1;


$command = "make -pr ";
$command .= join(" ", @ARGV);


my %id;
my $maxid=0;
print "digraph Dependencies {\n";
print "ratio=0.7;\n";

my %rules;

sub goodname($) {
    my $name = shift;
    $name =~ s/run_[0-9]+//g;
    $name =~ s/^.*\///;
    $name =~ s/^\s+//;
    $name =~ s/[\W]//g;    
    $name =~ s/fits//g;
    $name =~ s/^_/%_/;
    return $name;
}

foreach  (`$command`) {

    next if /^[#\t]/;
    next if (!/:/);
    next if (/[A-Z ]+[:=?]/);
    next if /Analysis\.mk/;
    chomp;
    s/=//g;


    ($source,$alltargets) = split(/:/,$_) or next;
    @targets = split(" ", $alltargets);

    $source = goodname($source);

    if (!defined($id{$source})) {
	$id{$source} = $maxid;
	$maxid++;
    }



    foreach $target (sort @targets) {

	$target = goodname($target);

	if (!defined($id{$target})) {
	    $id{$target} = $maxid;
	    $maxid++;
	}

	#$rules{"$id{$source} -> $id{$target}"}++;
	$rules{"$id{$target} -> $id{$source}"}++;

    }

}

if ($usenames) {
    foreach $name (sort keys %id) {
	$num = $id{$name};

	print "\t $num [label=\"$name\"";
	print " shape=folder" if ($name =~ /^\%/);
	print "];\n";
    }
}


foreach $rule (keys %rules) {

    $count = $rules{$rule};
    $width = $count /10 + 1;
    print "\t$rule";
    if ($count > 1) {
	print " [style=\"setlinewidth($width)\"]"
    }
    print ";\n";
}


print "}\n";
