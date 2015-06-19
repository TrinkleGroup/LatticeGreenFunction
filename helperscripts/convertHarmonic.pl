#!/usr/bin/env perl

use strict;
use warnings;
use Utils::PhonUtil;

open(INHARM, '<HARMONIC') or die "Error opening force matrix file 'HARMONIC'\n";
open(INPHON, '<INPHON') or die "Error opening 'INPHON'\n";
open(INPOSCAR, '<POSCAR') or die "Error opening 'POSCAR'\n";

# strip blanks and comments
my @inphondata = grep(!/(?:^#)|(?:^\s*$)/, <INPHON>);
close(INPHON);

my @ndims = ();

#get supercell dimensions
foreach my $line (@inphondata) {
	next if($line !~ /NDIM\s*=\s*(\d+)\s+(\d+)\s+(\d+)/);
	$ndims[0] = $1;
	$ndims[1] = $2;
	$ndims[2] = $3;
	last;
}

die "No supercell dimensions defined in INPHON!\n" if ($#ndims == 0);

print "$ndims[0] $ndims[1] $ndims[2]\n";
my $primcells = $ndims[0]*$ndims[1]*$ndims[2];

my $natomtypes;
#get number of atom types
foreach my $line (@inphondata) {
	next if($line !~ /NTYPES\s*=\s*(\d+)/);
	$natomtypes = $1;
	last;
}

my %latt = loadposcar(<INPOSCAR>);
close(INPOSCAR);

print $latt{'scale'}."\n";
print join(' ', @{$latt{'a1'}})."\n";
print join(' ', @{$latt{'a2'}})."\n";
print join(' ', @{$latt{'a3'}})."\n";
print join(' ', @{$latt{'natoms'}})."\n";

print "atompos:\n";
foreach my $atompos (@{$latt{'atompos'}}) {
	foreach (@{$atompos}) {
		print join(' ', @{$_}) . "\n";
	}
	print "=====\n";
}

$latt{'pnatoms'} = [];
$latt{'totalatoms'} = 0;

print "superatoms:\n";
foreach my $atom (@{$latt{'natoms'}}) {
	print "$atom\n";
	push(@{$latt{'pnatoms'}}, $atom/$primcells);
	$latt{'totalatoms'} += $atom/$primcells;
}
print "primatoms:\n";
foreach my $atom (@{$latt{'pnatoms'}}) {
	print "$atom\n";
}
print "total atoms: $latt{'totalatoms'}\n";

my %harmonic = loadharmonic(<INHARM>);
close(INHARM);

my $fullatoms = 0;
print "HARMONIC:\nRcounts:\n";
foreach (sort {$a <=> $b;} keys %{$harmonic{'rcount'}}) {
	print "$_: $harmonic{'rcount'}->{$_}\n";
	if($harmonic{'rcount'}->{$_} == $latt{'totalatoms'}**2) {
		$fullatoms++;
	}
}
print "Rvalues:\n";
foreach (sort {$a <=> $b;} keys %{$harmonic{'rhash'}}) {
	print "$_: $harmonic{'rhash'}->{$_}[0] $harmonic{'rhash'}->{$_}[1] $harmonic{'rhash'}->{$_}[2]\n";
}
print "Rmags:\n";
foreach (sort {$a <=> $b;} keys %{$harmonic{'rhash'}}) {
	print ("$_: " . sqrt($harmonic{'rhash'}->{$_}[0]**2 + $harmonic{'rhash'}->{$_}[1]**2 + $harmonic{'rhash'}->{$_}[2]**2) . "\n");
}
print "full atoms: $fullatoms\n";

open(DMAT, '>dmat.txt') or die 'Could not open dmat.txt for writing, stoppped';
print DMAT "DMAT for: $latt{'desc'}\n";
print DMAT ((scalar keys %{$harmonic{'rcount'}}) . " $latt{'totalatoms'}\n");
foreach my $rvec (sort {$a <=> $b;} keys %{$harmonic{'rcount'}}) {
	my @dmat = @{${$harmonic{'dmat'}}[$rvec]};
	print DMAT "$harmonic{'rhash'}->{$rvec}->[0] " .
		   "$harmonic{'rhash'}->{$rvec}->[1] " .
		   "$harmonic{'rhash'}->{$rvec}->[2]";
	for(my $at1=0; $at1 < $latt{'totalatoms'}; $at1++) {
		for(my $row=0; $row<3; $row++) {
			for(my $at2=0; $at2 < $latt{'totalatoms'}; $at2++) {
				for(my $col=0; $col<3; $col++) {
					print DMAT " ";
					print DMAT defined($dmat[$at2][$at1][$row][$col]) ? $dmat[$at2][$at1][$row][$col] : 0.0;
				}
			}
		}
	}
	print DMAT "\n";
}

close(DMAT);
