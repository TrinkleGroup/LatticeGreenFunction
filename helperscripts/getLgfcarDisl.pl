#!/usr/bin/env perl

use strict;
use warnings;

use POSIX qw(floor);

sub vecdiff($$) {
	my $atom1 = $_[0];
	my $atom2 = $_[1];

	my $atomdiff = [($$atom1[0] - $$atom2[0]),($$atom1[1] - $$atom2[1]),($$atom1[2] - $$atom2[2])];
	return($atomdiff);
}

sub loadcell($) {
	my $cellfilename = $_[0];
	my %cell;
	open(CELLFILE, "<$cellfilename") or die "Couldn't open Cell file: $cellfilename.";

	my $a0 = <CELLFILE>;
	chomp($a0);
	$a0 =~ s/\s*\#.*//;
	$cell{'a0'} = $a0;
	
	my $avec = [[],[],[]];

	for(my $idx=0; $idx<3; $idx++) {
		$_ = <CELLFILE>;
		chomp;
		s/\s*\#.*//;
		@{$$avec[$idx]} = split(/\s+/);
	}

	$cell{'avec'} = $avec;

	my $cc = <CELLFILE>;
	chomp($cc);
	$cc =~ s/\s*\#.*//;

	$cell{'cc'} = $cc;

	my $natoms = <CELLFILE>;
	chomp($natoms);
	$natoms =~ s/\s*\#.*//;

	$cell{'natoms'} = $natoms;

	my $apos = [];
	my $aposc = [];

	for(my $idx = 0; $idx < $natoms; $idx++) {
		my $atline = <CELLFILE>;
		chomp($atline);
		$atline =~ s/\s*\#.*//;
		@{$$apos[$idx]} = split(/\s+/,$atline);
		${$$aposc[$idx]}[0] = (${$$apos[$idx]}[0]*${$$avec[0]}[0]
		                      + ${$$apos[$idx]}[1]*${$$avec[1]}[0]
		                      + ${$$apos[$idx]}[2]*${$$avec[2]}[0])*$a0;
		${$$aposc[$idx]}[1] = (${$$apos[$idx]}[0]*${$$avec[0]}[1]
		                      + ${$$apos[$idx]}[1]*${$$avec[1]}[1]
		                      + ${$$apos[$idx]}[2]*${$$avec[2]}[1])*$a0;
		${$$aposc[$idx]}[2] = (${$$apos[$idx]}[0]*${$$avec[0]}[2]
		                      + ${$$apos[$idx]}[1]*${$$avec[1]}[2]
		                      + ${$$apos[$idx]}[2]*${$$avec[2]}[2])*$a0;
	}
	$cell{'apos'} = $apos;
	$cell{'aposc'} = $aposc;
	close(CELLFILE);
	return %cell;
}

sub lgfsort {
	return rveccomp($a->{'rvec'},$b->{'rvec'});
}

sub rveccomp($$) {
	my $arvec = $_[0];
	my $brvec = $_[1];
	return (abs($$arvec[0] - $$brvec[0]) < 1e-7 ?
	               ((abs($$arvec[1] - $$brvec[1]) < 1e-7) ?
	                    ((abs($$arvec[2] - $$brvec[2]) < 1e-7) ? 0 : $$arvec[2] <=> $$brvec[2])
			  : $$arvec[1] <=> $$brvec[1])
	               : $$arvec[0] <=> $$brvec[0]);
}

sub loadlgf($$) {
	my $lgffilename = $_[0];
	my $natoms = $_[1];

	open(LGFFILE, "<$lgffilename") or die "Couldn't open LGF file: $lgffilename.";

	my $lgf = [];

	while(my $line = <LGFFILE>) {
		chomp($line);
		my $rvec = [];
		@{$rvec} = split(/\s+/,$line);
		my $lgfmat = [];

		for(my $jdx = 0; $jdx < 3*$natoms; $jdx++) {
			my $lgfline = <LGFFILE>;
			chomp($lgfline);
			my @lsplt = split(/\s+/,$lgfline);
			$$lgfmat[$jdx] = [@lsplt];
		}
		my $lgfl;
		$lgfl->{'rvec'} = $rvec;
		$lgfl->{'lgfmat'} = $lgfmat;
		push(@{$lgf}, $lgfl);
	}
	close(LGFFILE);
	@{$lgf} = sort lgfsort @{$lgf};
	return $lgf;
}

sub searchlgf($$$$$) {
	my $lgf = shift;
	my $rvec = shift;
	my $atom1 = shift;
	my $atom2 = shift;
	my $rdiffs = shift;
	$atom1--;
	$atom2--;

	my $beg = 0;
	my $end = $#$lgf;

	my $lgfrec;

	while($beg <= $end) {
		my $idx = int (($end+$beg)/2);
		my $searchvec = $$lgf[$idx]->{'rvec'};
		my $compval = rveccomp($rvec,$searchvec);
		if($compval < 0) {
			$end = $idx-1;
		} elsif($compval > 0) {
			$beg = $idx+1;
		} else {
			$lgfrec = $$lgf[$idx];
			my $lgf3x3 = [];
			for(my $idx = 0; $idx < 3; $idx++) {
				for(my $jdx = 0; $jdx < 3; $jdx++) {
					$$lgf3x3[3*$idx+$jdx] = $lgfrec->{'lgfmat'}->[3*$atom1+$idx]->[3*$atom2+$jdx];
				}
			}
			return $lgf3x3;
		}
	}

	die "Could not find LGF record: $$rvec[0] $$rvec[1] $$rvec[2]\n";

}

die "Usage: $0 lgf.out cellfile lgfpairing.out LGF_Description\n" if(scalar @ARGV != 4);

my $lgffilename = shift @ARGV;
my $cellfilename = shift @ARGV;
my $lgfpfilename = shift @ARGV;
my $lgfdesc = shift @ARGV;

my %cell = loadcell($cellfilename);
my @Rdiff = ();
for(my $idx = 0; $idx < $cell{'natoms'}; $idx++) {
	$Rdiff[$idx] = [];
	for(my $jdx = 0; $jdx < $cell{'natoms'}; $jdx++) {
		$Rdiff[$idx]->[$jdx] = vecdiff($cell{'aposc'}->[$idx],$cell{'aposc'}->[$jdx]);
	}
}

my $lgfdata = loadlgf($lgffilename, $cell{'natoms'});

my @lgfcar = ();

open(LGFPFILE, "<$lgfpfilename") or die "Could not open LGF pairing file: $lgfpfilename\n";

#loop through all atom pairs
foreach my $lgfp (<LGFPFILE>) {
	my @lgfpsp = split(/\s+/, $lgfp);
	my $atom2 = $lgfpsp[0];
	my $atom1 = $lgfpsp[1];
	my $Rvec = [$lgfpsp[2],$lgfpsp[3],$lgfpsp[4]];
	my $atomtype2 = $lgfpsp[5];
	my $atomtype1 = $lgfpsp[6];

	my $lgf3x3 = searchlgf($lgfdata, $Rvec, $atomtype1, $atomtype2, \@Rdiff);

	my $lgfline = [$atom2,$atom1, $lgf3x3];
	push(@lgfcar, $lgfline);
}

close(LGFPFILE);

@lgfcar = sort { return ($a->[0] == $b->[0] ? $a->[1] <=> $b->[1] : $a->[0] <=> $b->[0]) } @lgfcar;

my $begrange = $lgfcar[0];
my $endrange = $lgfcar[$#lgfcar];
my $lgfelems = scalar @lgfcar;
my $r2beg = $$begrange[0];
my $r2end = $$endrange[0];
my $dislbeg = $$begrange[1];
my $dislend = $$endrange[1];

print "$lgfdesc\n";
print "$r2beg $r2end $dislbeg $dislend $lgfelems\n";

foreach my $elem (@lgfcar) {
	print "$$elem[0] $$elem[1]";
	for(my $idx = 0; $idx < 9; $idx++) {
		printf(" %+.14e", $$elem[2]->[$idx]);
	}
	print "\n";
}
