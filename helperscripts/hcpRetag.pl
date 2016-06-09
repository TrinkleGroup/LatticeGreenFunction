#!/usr/bin/env perl

use strict;
use warnings;

sub trimcomm(\$) {
	my $ts = shift;
	$$ts =~ s/\s*\#.*//;
	$$ts =~ s/^\s+//;
	return $$ts;
}

sub vecdiff($$) {
	my $atom1 = shift;
	my $atom2 = shift;

	my $atomdiff = [($$atom1[0] - $$atom2[0]),($$atom1[1] - $$atom2[1]),($$atom1[2] - $$atom2[2])];
	return($atomdiff);
}

sub vecmag2($) {
	my $vec = shift;
	my $mag2 = $$vec[0]**2 + $$vec[1]**2 + $$vec[2]**2;
	return $mag2;
}

sub rveccomp($$) {
	my $arvec = shift;
	my $brvec = shift;
	return (abs($$arvec[0] - $$brvec[0]) < 1e-7 ?
	               ((abs($$arvec[1] - $$brvec[1]) < 1e-7) ?
	                    ((abs($$arvec[2] - $$brvec[2]) < 1e-7) ? 0 : $$arvec[2] <=> $$brvec[2])
			  : $$arvec[1] <=> $$brvec[1])
	               : $$arvec[0] <=> $$brvec[0]);
}

sub rgridadd($$) {
	my $rlist = shift;
	my $rvec = shift;

	my $beg = 0;
	my $end = $#$rlist;

	if($end == -1) {
		unshift(@{$rlist}, $rvec);
	}

	#boundary cases:
	my $begvec = $rlist->[$beg];
	my $begcompval = rveccomp($rvec,$begvec);
	if($begcompval < 0) {
		unshift(@{$rlist}, $rvec);
		return;
	} elsif($begcompval == 0) {
		return;
	}

	my $endvec = $rlist->[$end];
	my $endcompval = rveccomp($rvec,$endvec);
	if($endcompval > 0) {
		push(@{$rlist}, $rvec);
		return;
	} elsif($endcompval == 0) {
		return;
	}

	while($beg < $end) {
		my $idx = int (($end+$beg)/2);
		my $searchvec = $rlist->[$idx];
		my $compval = rveccomp($rvec,$searchvec);
		if($compval == 0) {
			return; #found vector;
		} elsif($compval < 0) {
			$end = $idx;
		} elsif($beg == $idx) {
			last;
		} else {
			$beg = $idx;
		}
	}
	#not found: add new r vector
	splice(@{$rlist}, $end, 0, $rvec);
}

sub loadcell($) {
	my $cellfilename = shift;
	my %cell;
	open(CELLFILE, "<$cellfilename") or die "Couldn't open Cell file: $cellfilename.";

	my $a0 = <CELLFILE>;
	chomp($a0);
	trimcomm($a0);
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
	trimcomm($cc);

	$cell{'cc'} = $cc;

	my $natoms = <CELLFILE>;
	chomp($natoms);
	trimcomm($natoms);

	$cell{'natoms'} = $natoms;

	my $apos = [];
	my $aposc = [];

	for(my $idx = 0; $idx < $natoms; $idx++) {
		my $atline = <CELLFILE>;
		chomp($atline);
		trimcomm($atline);
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

sub vmatmult($$) {
	my $vec = shift;
	my $mat = shift;
	my $retvec = [];
	$$retvec[0] = $$vec[0]*$$mat[0]->[0]
	              + $$vec[1]*$$mat[1]->[0]
	              + $$vec[2]*$$mat[2]->[0];
	$$retvec[1] = $$vec[0]*$$mat[0]->[1]
	              + $$vec[1]*$$mat[1]->[1]
	              + $$vec[2]*$$mat[2]->[1];
	$$retvec[2] = $$vec[0]*$$mat[0]->[2]
	              + $$vec[1]*$$mat[1]->[2]
	              + $$vec[2]*$$mat[2]->[2];
	return $retvec;
}

sub crossprod($$) {
	my $vec1 = shift;
	my $vec2 = shift;
	my $retvec = [];
	$$retvec[0] = $$vec1[1]*$$vec2[2] - $$vec1[2]*$$vec2[1];
	$$retvec[1] = $$vec1[2]*$$vec2[0] - $$vec1[0]*$$vec2[2];
	$$retvec[2] = $$vec1[0]*$$vec2[1] - $$vec1[1]*$$vec2[0];

	return $retvec;
}

sub loaddisl_rot($\%) {
	my $dislfilename = shift;
	my $cell = shift;
	open (DISLFILE, "<$dislfilename") or die "Could not open dislocation file: $dislfilename\n";
	my $tdir = <DISLFILE>;
	chomp($tdir);
	trimcomm($tdir);
	my @tvec_u = split(/\s+/, $tdir);

	my $bdir = <DISLFILE>;
	chomp($bdir);
	trimcomm($bdir);
	my @bvec_u = split(/\s+/, $bdir);

	my $mdir = <DISLFILE>;
	chomp($mdir);
	trimcomm($mdir);
	my @mvec_u = split(/\s+/, $mdir);

	close(DISLFILE);

	my $avec = $cell->{'avec'};

	my $tvec_c = vmatmult(\@tvec_u, $avec);
	my $mvec_c = vmatmult(\@mvec_u, $avec);

	my $tvecmag = sqrt(vecmag2($tvec_c));
	my $mvecmag = sqrt(vecmag2($mvec_c));

	my $tvec_n = [];
	$$tvec_n[0] = $$tvec_c[0]/$tvecmag;
	$$tvec_n[1] = $$tvec_c[1]/$tvecmag;
	$$tvec_n[2] = $$tvec_c[2]/$tvecmag;

	my $mvec_n = [];
	$$mvec_n[0] = $$mvec_c[0]/$mvecmag;
	$$mvec_n[1] = $$mvec_c[1]/$mvecmag;
	$$mvec_n[2] = $$mvec_c[2]/$mvecmag;

	my $nvec_n = crossprod($tvec_n,$mvec_n);

	my $rotmat = [$mvec_n,$nvec_n,$tvec_n];

	return $rotmat;
}


die "Usage: $0 cellfile dislfile bulk.xyz disl.xyz newbulk.xyz newdisl.xyz\n" if(scalar @ARGV != 6);

my $cellfilename = shift @ARGV;
my $dislfilename = shift @ARGV;
my $xyzfilename = shift @ARGV;
my $dislxyzfilename = shift @ARGV;
my $newxyzfilename = shift @ARGV;
my $newdislxyzfilename = shift @ARGV;

my %cell = loadcell($cellfilename);
my $R12diff = vecdiff($cell{'aposc'}->[0],$cell{'aposc'}->[1]);
my $R21diff = vecdiff($cell{'aposc'}->[1],$cell{'aposc'}->[0]);

my $rotmat = loaddisl_rot($dislfilename, %cell);

open(XYZFILE, "<$xyzfilename") or die "Couldn't open bulk.xyz file: $xyzfilename.";
open(DISLXYZFILE, "<$dislxyzfilename") or die "Couldn't open disl.xyz file: $dislxyzfilename.";
open(NEWXYZFILE, ">$newxyzfilename") or die "Couldn't open new bulk.xyz file: $newxyzfilename.";
open(NEWDISLXYZFILE, ">$newdislxyzfilename") or die "Couldn't open new disl.xyz file: $newdislxyzfilename.";

my $natoms = <XYZFILE>;
chomp($natoms);
print NEWXYZFILE "$natoms\n";

my $comment = <XYZFILE>;
chomp($comment);
print NEWXYZFILE "$comment\n";

my $d_natoms = <DISLXYZFILE>;
chomp($d_natoms);
print NEWDISLXYZFILE "$d_natoms\n";

my $d_comment = <DISLXYZFILE>;
chomp($d_comment);
print NEWDISLXYZFILE "$d_comment\n";



my @atoms = ( [[],[],[]] , [[],[],[]] );

# Hexagonal specific for atom segregation:
my $cval = $cell{'avec'}->[2]->[2] * $cell{'a0'};

while(my $line = <XYZFILE>) {
	chomp($line);
	my @atompos = split(/\s+/,$line);
	my $regiontype = shift(@atompos);
	my $atomtype;
	#absolute value of c axis unit cell coordinate
	my $c_unit = abs($atompos[2]/$cval - 0.25);
	if(abs($c_unit - int($c_unit+0.5)) < 1e-5) {
		$atomtype = 1;
	} elsif(abs($c_unit - int($c_unit) - 0.5) < 1e-5) {
		$atomtype = 2;
	} else {
		die "Unknown atom type for segregation: " . ($c_unit - int($c_unit+0.5)) . "\n";
	}

	my $d_line = <DISLXYZFILE>;
	chomp($d_line);
	my @d_atompos = split(/\s+/,$d_line);
	my $d_regiontype = shift(@d_atompos);

	print NEWXYZFILE "${regiontype}_${atomtype} " . join(' ', @atompos) . "\n";
	print NEWDISLXYZFILE "${d_regiontype}_${atomtype} " . join(' ', @d_atompos) . "\n";
}

close(XYZFILE);
close(DISLXYZFILE);
close(NEWXYZFILE);
close(NEWDISLXYZFILE);
