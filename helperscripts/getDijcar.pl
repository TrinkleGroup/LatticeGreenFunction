#!/usr/bin/env perl

use strict;
use warnings;

use POSIX qw(floor);

sub trimcomm(\$) {
	my $ts = shift;
	$$ts =~ s/\s*\#.*//;
	$$ts =~ s/^\s+//;
	return $$ts;
}

sub vecdiff($$) {
	my $atom1 = $_[0];
	my $atom2 = $_[1];

	my $atomdiff = [($$atom1[0] - $$atom2[0]),($$atom1[1] - $$atom2[1]),($$atom1[2] - $$atom2[2])];
	return($atomdiff);
}

sub vecmag2($) {
	my $vec = $_[0];
	my $mag2 = $$vec[0]**2 + $$vec[1]**2 + $$vec[2]**2;
	return $mag2;
}

sub loadcell($) {
	my $cellfilename = $_[0];
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

sub dmat2dsort {
	return rveccomp2d($a->{'rvec'},$b->{'rvec'});
}

sub rveccomp2d($$) {
	my $arvec = $_[0];
	my $brvec = $_[1];
	return ((abs($$arvec[1] - $$brvec[1]) < 1e-7) ?
	               ((abs($$arvec[2] - $$brvec[2]) < 1e-7) ? 0 : $$arvec[2] <=> $$brvec[2])
	             : $$arvec[1] <=> $$brvec[1]);
}

sub loaddmat2d($$) {
	my $dmatfilename = shift;
	my $natoms = shift;

	open(DMATFILE, "<$dmatfilename") or die "Couldn't open 2d Dmat file: $dmatfilename.";

	my $comm = <DMATFILE>;
	my $numvecs = <DMATFILE>;

	my $dmat = [];

	while(my $line = <DMATFILE>) {
		chomp($line);
		my @linespl = split(/\s+/,$line);
		my $rvec = [];
		$$rvec[0] = shift(@linespl);
		$$rvec[1] = shift(@linespl);
		$$rvec[2] = shift(@linespl);
		my $dmatelem = [];

		for(my $idx = 0; $idx < 3*$natoms; $idx++) {
			$$dmatelem[$idx] = [];
			for(my $jdx = 0; $jdx < 3*$natoms; $jdx++) {
				$$dmatelem[$idx]->[$jdx] = $linespl[$idx*3*$natoms + $jdx];
			}
		}
		my $dmatl;
		$dmatl->{'rvec'} = $rvec;
		$dmatl->{'dmat'} = $dmatelem;
		push(@{$dmat}, $dmatl);
	}
	close(DMATFILE);
	@{$dmat} = sort dmat2dsort @{$dmat};
	return $dmat;
}

sub searchdmat($$$$) {
	my $dmat = shift;
	my $rvec = shift;
	my $atom1 = shift;
	my $atom2 = shift;

	my $beg = 0;
	my $end = $#$dmat;

	my $dmatrec;
	my $found = 0;

	while($beg < $end) {
		my $idx = int (($end+$beg)/2);
		my $searchvec = $$dmat[$idx]->{'rvec'};
		my $compval = rveccomp2d($rvec,$searchvec);
		if($compval == 0) {
			$dmatrec = $$dmat[$idx];
			$found = 1;
			last;
		} elsif($compval < 0) {
			$end = $idx;
		} elsif($beg == $idx) {
			$searchvec = $$dmat[$end]->{'rvec'};
			$compval = rveccomp2d($rvec,$searchvec);
			if($compval == 0) {
				$dmatrec = $$dmat[$idx];
				$found = 1;
			}
			last;
		} else {
			$beg = $idx;
		}
	}

	my $dmat3x3 = [];

	if($found == 0) {
		#die "Could not find 2d DMAT element for rvec: $$rvec[0], $$rvec[1], $$rvec[2]!\n";
		# Outside cut-off
		for(my $idx = 0; $idx < 3; $idx++) {
			for(my $jdx = 0; $jdx < 3; $jdx++) {
				$$dmat3x3[3*$idx+$jdx] = 0.0;
			}
		}
	} else {
		for(my $idx = 0; $idx < 3; $idx++) {
			for(my $jdx = 0; $jdx < 3; $jdx++) {
				$$dmat3x3[3*$idx+$jdx] = $dmatrec->{'dmat'}->[3*$atom1+$idx]->[3*$atom2+$jdx];
			}
		}
	}
	return $dmat3x3;
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
	my @tvec_u = split(/\s/, $tdir);

	my $bdir = <DISLFILE>;
	chomp($bdir);
	trimcomm($bdir);
	my @bvec_u = split(/\s/, $bdir);

	my $mdir = <DISLFILE>;
	chomp($mdir);
	trimcomm($mdir);
	my @mvec_u = split(/\s/, $mdir);

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

die "Usage: $0 2ddmatfile cellfile dislfile bulk.xyz Region_1_element Region_2_element Region_3_element.\n" if(scalar @ARGV != 7);

my $dmatfilename = shift @ARGV;
my $cellfilename = shift @ARGV;
my $dislfilename = shift @ARGV;
my $xyzfilename = shift @ARGV;
my @elems = @ARGV;

my %cell = loadcell($cellfilename);
my $R12diff = vecdiff($cell{'aposc'}->[0],$cell{'aposc'}->[1]);
my $R21diff = vecdiff($cell{'aposc'}->[1],$cell{'aposc'}->[0]);
printf(STDERR "1-2: %+.12f %+.12f %+.12f\n", $$R12diff[0], $$R12diff[1], $$R12diff[2]);
printf(STDERR "2-1: %+.12f %+.12f %+.12f\n", $$R21diff[0], $$R21diff[1], $$R21diff[2]);

my $dmatdata = loaddmat2d($dmatfilename, $cell{'natoms'});

my $rotmat = loaddisl_rot($dislfilename, %cell);

open(XYZFILE, "<$xyzfilename") or die "Couldn't open xyz file: $xyzfilename.";

my $natoms = <XYZFILE>;
chomp($natoms);

my $comment = <XYZFILE>;
chomp($comment);

my @atoms = ( [[],[],[]] , [[],[],[]] );

my $atomnum = 1;

# Hexagonal specific for atom segregation:
my $cval = $cell{'avec'}->[2]->[2] * $cell{'a0'};

while(my $line = <XYZFILE>) {
	chomp($line);
	my @aposrot = split(/\s+/,$line);
	my $regiontype = shift(@aposrot);
	my $atompos = vmatmult(\@aposrot,$rotmat);
	push(@$atompos, $atomnum);
	my $atomtype;
	#absolute value of c axis unit cell coordinate
	my $c_unit = abs($$atompos[2]/$cval - 0.25);
	if(abs($c_unit - int($c_unit+0.5)) < 1e-5) {
		$atomtype = 0;
	} elsif(abs($c_unit - int($c_unit) - 0.5) < 1e-5) {
		$atomtype = 1;
	} else {
		die "Unknown atom type for segregation: " . ($c_unit - int($c_unit+0.5)) . "\n";
	}
	my $regionnum;
	if($elems[0] eq $regiontype) {
		$regionnum = 0;
	} elsif($elems[1] eq $regiontype) {
		$regionnum = 1;
	} elsif($elems[2] eq $regiontype) {
		$regionnum = 2;
	} else {
		die "Error: Invalid atom in xyz file: $xyzfilename!\n";
	}
	push(@{$atoms[$atomtype]->[$regionnum]}, $atompos);
	$atomnum++;
}

close(XYZFILE);

my @dijcar = ();

#loop through all region 2 atoms
for(my $atomtype2 = 0; $atomtype2 < 2; $atomtype2++) {
	foreach my $atompos2 (@{$atoms[$atomtype2]->[1]}) {
		#loop through all atoms.
		for(my $atomtype = 0; $atomtype < 2; $atomtype++) {
			for(my $regionnum = 0; $regionnum < 3; $regionnum++) {
				foreach my $atompos (@{$atoms[$atomtype]->[$regionnum]}) {
					my $adiff = vecdiff($atompos,$atompos2);
					my $rdiff;
					if($atomtype != $atomtype2) {
						if($atomtype == 0) {
							$rdiff = vecdiff($adiff, $R12diff);
						} else {
							$rdiff = vecdiff($adiff, $R21diff);
						}
					} else {
						$rdiff = $adiff;
					}
					if($$atompos[3] == 85 && $$atompos2[3] == 108) {
						printf(STDERR "85: %d, 108: %d\n", $atomtype, $atomtype2);
						printf(STDERR "85,108: %+.12f %+.12f %+.12f\n", $$adiff[0], $$adiff[1], $$adiff[2]);
						printf(STDERR "85,108: %+.12f %+.12f %+.12f\n", $$rdiff[0], $$rdiff[1], $$rdiff[2]);
						printf(STDERR "85: %+.12f %+.12f %+.12f\n", $$atompos[0], $$atompos[1], $$atompos[2]);
						printf(STDERR "108: %+.12f %+.12f %+.12f\n", $$atompos2[0], $$atompos2[1], $$atompos2[2]);
					}
					my $dmat3x3 = searchdmat($dmatdata, $rdiff, $atomtype, $atomtype2);

					my $dmatline = [$atompos2->[3],$atompos->[3], $dmat3x3];
					push(@dijcar, $dmatline);
				}
			}
		}
	}
}

@dijcar = sort { return ($a->[0] == $b->[0] ? $a->[1] <=> $b->[1] : $a->[0] <=> $b->[0]) } @dijcar;

foreach my $elem (@dijcar) {
	print "$$elem[0] $$elem[1]";
	for(my $idx = 0; $idx < 9; $idx++) {
		printf(" %+.14f", $$elem[2]->[$idx]);
	}
	print "\n";
}
