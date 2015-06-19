#!/usr/bin/env perl

use strict;
use warnings;

use Math::Round qw(round);

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

sub dotprod($$) {
	my $vec1 = shift;
	my $vec2 = shift;
	my $retval = $$vec1[0]*$$vec2[0] + $$vec1[1]*$$vec2[1] + $$vec1[2]*$$vec2[2];
	return $retval;
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

sub rveccomp($$) {
	my $arvec = shift;
	my $brvec = shift;
	return (abs($$arvec[0] - $$brvec[0]) < 1e-7 ?
	               ((abs($$arvec[1] - $$brvec[1]) < 1e-7) ?
	                    ((abs($$arvec[2] - $$brvec[2]) < 1e-7) ? 0 : $$arvec[2] <=> $$brvec[2])
			  : $$arvec[1] <=> $$brvec[1])
	               : $$arvec[0] <=> $$brvec[0]);
}

sub rmags {
	my $amag2 = vecmag2($a);
	my $bmag2 = vecmag2($b);
	return ( abs($amag2 - $bmag2) < 1e-7 ? 0 : $amag2 <=> $bmag2 );
}

sub searchnearest($$) {
	my $cell = shift;
	my $rvec = shift;

	my $guess = vmatmult($rvec, $cell->{'bvec'});

	my $roundguess = [];
	$$roundguess[0] = round($$guess[0]);
	$$roundguess[1] = round($$guess[1]);
	$$roundguess[2] = round($$guess[2]);

	my $cartguess = vmatmult($roundguess, $cell->{'avec'});
	$$cartguess[0] *= $cell->{'a0'};
	$$cartguess[1] *= $cell->{'a0'};
	$$cartguess[2] *= $cell->{'a0'};

	my $minvec = $cartguess;
	my $mindiff2 = vecmag2(vecdiff($rvec, $minvec));

	foreach my $vec (@{$cell->{'nnvecs'}}) {
		my $vecadd = [];
		$$vecadd[0] = $$vec[0] + $$cartguess[0];
		$$vecadd[1] = $$vec[1] + $$cartguess[1];
		$$vecadd[2] = $$vec[2] + $$cartguess[2];
		my $vecdiff2 = vecmag2(vecdiff($rvec, $vecadd));
		if($vecdiff2 < $mindiff2) {
			$mindiff2 = $vecdiff2;
			$minvec = $vecadd;
		}
	}

	return $minvec;
}

sub loadcell($) {
	my $cellfilename = shift;
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

	# reciprocal without the two pi:

	my $cellvol = $a0*dotprod($$avec[0], crossprod($$avec[1], $$avec[2]));
	my $b1 = crossprod($$avec[1],$$avec[2]);
	my $b2 = crossprod($$avec[2],$$avec[0]);
	my $b3 = crossprod($$avec[0],$$avec[1]);
	$$b1[0] /= $cellvol; $$b1[1] /= $cellvol; $$b1[2] /= $cellvol;
	$$b2[0] /= $cellvol; $$b2[1] /= $cellvol; $$b2[2] /= $cellvol;
	$$b3[0] /= $cellvol; $$b3[1] /= $cellvol; $$b3[2] /= $cellvol;

	my $bvec = [[],[],[]];
	$$bvec[0]->[0] = $$b1[0]; $$bvec[0]->[1] = $$b2[0]; $$bvec[0]->[2] = $$b3[0];
	$$bvec[1]->[0] = $$b1[1]; $$bvec[1]->[1] = $$b2[1]; $$bvec[1]->[2] = $$b3[1];
	$$bvec[2]->[0] = $$b1[2]; $$bvec[2]->[1] = $$b2[2]; $$bvec[2]->[2] = $$b3[2];

	$cell{'bvec'} = $bvec;


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

	my @nnvecs = ([1,0,0],
		      [1,1,0],
		      [0,1,0],
		      [-1,0,0],
		      [-1,-1,0],
		      [0,-1,0],
	              [0,0,1],
	              [1,0,1],
		      [1,1,1],
		      [0,1,1],
		      [-1,0,1],
		      [-1,-1,1],
		      [0,-1,1],
	              [0,0,-1],
	              [1,0,-1],
		      [1,1,-1],
		      [0,1,-1],
		      [-1,0,-1],
		      [-1,-1,-1],
		      [0,-1,-1]);

#		      [1,2,-1],
#		      [-1,2,-1],
#		      [-2,-1,-1],
#		      [-1,-2,-1],
#		      [1,-1,-1],
#		      [2,1,-1]);
#		      [1,2,0],
#		      [-1,2,0],
#		      [-2,-1,0],
#		      [-1,-2,0],
#		      [1,-1,0],
#		      [2,1,0],
#		      [1,2,1],
#		      [-1,2,1],
#		      [-2,-1,1],
#		      [-1,-2,1],
#		      [1,-1,1],
#		      [2,1,1],
	
	my $cnnvecs = [];
	foreach my $dnnvec (@nnvecs) {
		my $nvec = vmatmult($dnnvec, $avec);
		$$nvec[0] *= $a0;
		$$nvec[1] *= $a0;
		$$nvec[2] *= $a0;
		push(@{$cnnvecs}, $nvec);
	}
	
	$cell{'nnvecs'} = $cnnvecs;


	return %cell;
}

my @tvec = ();

sub loaddisl_rot($\%) {
	my $dislfilename = shift;
	my $cell = shift;
	open (DISLFILE, "<$dislfilename") or die "Could not open dislocation file: $dislfilename\n";
	my $tdir = <DISLFILE>;
	chomp($tdir);
	$tdir =~ s/\s*\#.*//;
	my @tvec_u = split(/\s+/, $tdir);

	my $bdir = <DISLFILE>;
	chomp($bdir);
	$bdir =~ s/\s*\#.*//;
	my @bvec_u = split(/\s+/, $bdir);

	my $mdir = <DISLFILE>;
	chomp($mdir);
	$mdir =~ s/\s*\#.*//;
	my @mvec_u = split(/\s+/, $mdir);

	close(DISLFILE);

	my $avec = $cell->{'avec'};

	my $tvec_c = vmatmult(\@tvec_u, $avec);
	$tvec[0] = $$tvec_c[0] * $cell->{'a0'};
	$tvec[1] = $$tvec_c[1] * $cell->{'a0'};
	$tvec[2] = $$tvec_c[2] * $cell->{'a0'};

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


die "Usage: $0 cellfile dislfile disl.xyz Region_1_element Region_2_element Region_3_element Rvec_output LGF_pairing_output.\n" if(scalar @ARGV != 8);

my $cellfilename = shift @ARGV;
my $dislfilename = shift @ARGV;
my $dislxyzfilename = shift @ARGV;
my $lgfpfilename = pop @ARGV;
my $rvecfilename = pop @ARGV;
my @elems = @ARGV;

my %cell = loadcell($cellfilename);
my @Rdiff = ();
for(my $idx = 0; $idx < $cell{'natoms'}; $idx++) {
	$Rdiff[$idx] = [];
	for(my $jdx = 0; $jdx < $cell{'natoms'}; $jdx++) {
		$Rdiff[$idx]->[$jdx] = vecdiff($cell{'aposc'}->[$idx],$cell{'aposc'}->[$jdx]);
	}
}

my $rotmat = loaddisl_rot($dislfilename, %cell);

open(DISLXYZFILE, "<$dislxyzfilename") or die "Couldn't open disl.xyz file: $dislxyzfilename.";

my $d_natoms = <DISLXYZFILE>;
chomp($d_natoms);

my $d_comment = <DISLXYZFILE>;
chomp($d_comment);

my @d_atoms = ([],[],[]);

my $d_atomnum = 1;

# Hexagonal specific for atom segregation:
# my $cval = $cell{'avec'}->[2]->[2] * $cell{'a0'};

while(my $line = <DISLXYZFILE>) {
	chomp($line);
	my @aposrot = split(/\s+/,$line);
	my $regiontype = shift(@aposrot);
	my $atompos = vmatmult(\@aposrot,$rotmat);
	push(@$atompos, $d_atomnum);

	my $regionnum;
	if($regiontype eq $elems[0]) {
		$regionnum = 0;
	} elsif($regiontype eq $elems[1]) {
		$regionnum = 1;
	} elsif($regiontype eq $elems[2]) {
		$regionnum = 2;
	} else {
		die "Error: Invalid atom in xyz file: $dislxyzfilename!\n";
	}
	push(@{$d_atoms[$regionnum]}, $atompos);
	$d_atomnum++;
}

close(DISLXYZFILE);

open(LGFPFILE, ">$lgfpfilename") or die "Could not open LGF pairing output file: $lgfpfilename.\n";
my @newrgrid = ();
foreach my $atompos2 (@{$d_atoms[1]}) {
	for(my $regionnum = 0; $regionnum < 3; $regionnum++) {
		foreach my $atompos (@{$d_atoms[$regionnum]}) {
			my $adiff = vecdiff($atompos,$atompos2);
			my $nearestlvec0 = searchnearest(\%cell, $adiff);
			my $mag0 = vecmag2(vecdiff($adiff,$nearestlvec0));
			my $rdiff;
			$rdiff = vecdiff($adiff, $Rdiff[0]->[1]);
			my $nearestlvec12 = searchnearest(\%cell, $rdiff);
			my $mag12 = vecmag2(vecdiff($rdiff,$nearestlvec12));
			$rdiff = vecdiff($adiff, $Rdiff[1]->[0]);
			my $nearestlvec21 = searchnearest(\%cell, $rdiff);
			my $mag21 = vecmag2(vecdiff($rdiff,$nearestlvec21));
			if($mag0 <= $mag12 && $mag0 <= $mag21) {
				push(@newrgrid, $nearestlvec0);
				printf LGFPFILE "%d %d %+.12f %+.12f %+.12f 1 1\n", $$atompos2[3], $$atompos[3], $$nearestlvec0[0], $$nearestlvec0[1], $$nearestlvec0[2];
				printf STDERR "%.4f %d %d 1 1\n", sqrt($mag0), $$atompos2[3], $$atompos[3];
			} elsif($mag12 <= $mag21) {
				push(@newrgrid, $nearestlvec12);
				printf LGFPFILE "%d %d %+.12f %+.12f %+.12f 2 1\n", $$atompos2[3], $$atompos[3], $$nearestlvec12[0], $$nearestlvec12[1], $$nearestlvec12[2];
				printf STDERR "%.4f %d %d 2 1\n", sqrt($mag12), $$atompos2[3], $$atompos[3];
			} else {
				push(@newrgrid, $nearestlvec21);
				printf LGFPFILE "%d %d %+.12f %+.12f %+.12f 1 2\n", $$atompos2[3], $$atompos[3], $$nearestlvec21[0], $$nearestlvec21[1], $$nearestlvec21[2];
				printf STDERR "%.4f %d %d 1 2\n", sqrt($mag21), $$atompos2[3], $$atompos[3];
			}
		}
	}
}
close(LGFPFILE);

@newrgrid = sort rveccomp @newrgrid;

my @uniqrgrid = ();
push(@uniqrgrid, $newrgrid[0]);

foreach my $rvec (@newrgrid) {
	if(rveccomp($rvec, $uniqrgrid[$#uniqrgrid]) != 0) {
		push(@uniqrgrid, $rvec);
	}
}

open(RVECFILE, ">$rvecfilename") or die "Could not open rvec output file: $rvecfilename.\n";

print RVECFILE scalar @uniqrgrid . " 1.0 C\n";
foreach my $rpoint (@uniqrgrid) {
	printf(RVECFILE "%+.12f %+.12f %+.12f\n", $$rpoint[0], $$rpoint[1], $$rpoint[2]);
}
close(RVECFILE);
