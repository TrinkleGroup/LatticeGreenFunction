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

sub rveccomp2d($$) {
	my $arvec = shift;
	my $brvec = shift;
	return (abs($$arvec[0] - $$brvec[0]) < 1e-7 ?
	               ((abs($$arvec[1] - $$brvec[1]) < 1e-7) ?
	                    0 : $$arvec[1] <=> $$brvec[1])
	               : $$arvec[0] <=> $$brvec[0]);
}

sub dmat2dsort {
	return rveccomp2d($a->{'rvec'},$b->{'rvec'});
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
		my $dmatl = {};
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
		for(my $idx = 0; $idx < 9; $idx++) {
			$$dmat3x3[$idx] = 0.0;
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

die "Usage: $0 dmat2d.out cellfile lgfpairing.out.\n" if(scalar @ARGV != 3);

my $dmatfilename = shift @ARGV;
my $cellfilename = shift @ARGV;
my $lgfpfilename = shift @ARGV;

my %cell = loadcell($cellfilename);
my @Rdiff = ();
for(my $idx = 0; $idx < $cell{'natoms'}; $idx++) {
	$Rdiff[$idx] = [];
	for(my $jdx = 0; $jdx < $cell{'natoms'}; $jdx++) {
		$Rdiff[$idx]->[$jdx] = vecdiff($cell{'aposc'}->[$idx],$cell{'aposc'}->[$jdx]);
	}
}

my $dmat2ddata = loaddmat2d($dmatfilename, $cell{'natoms'});

my @dijcar = ();

open(LGFPFILE, "<$lgfpfilename") or die "Could not open LGF pairing file: $lgfpfilename\n";

#loop through all atom pairs
foreach my $lgfp (<LGFPFILE>) {
	my @lgfpsp = split(/\s+/, $lgfp);
	# R -> -R for the dmat
	my $atom1 = $lgfpsp[0];
	my $atom2 = $lgfpsp[1];
	my $Rvec = [-$lgfpsp[2],-$lgfpsp[3],-$lgfpsp[4]];
	my $atomtype1 = $lgfpsp[5];
	my $atomtype2 = $lgfpsp[6];

	my $dmat3x3 = searchdmat($dmat2ddata, $Rvec, $atomtype1-1, $atomtype2-1);

	my $dmatline = [$atom2,$atom1,$dmat3x3];
	push(@dijcar, $dmatline);
}

close(LGFPFILE);

@dijcar = sort { return ($a->[0] == $b->[0] ? $a->[1] <=> $b->[1] : $a->[0] <=> $b->[0]) } @dijcar;

foreach my $elem (@dijcar) {
	print "$$elem[0] $$elem[1]";
	for(my $idx = 0; $idx < 9; $idx++) {
		printf(" %+.14f", ${$$elem[2]}[$idx]);
	}
	print "\n";
}
