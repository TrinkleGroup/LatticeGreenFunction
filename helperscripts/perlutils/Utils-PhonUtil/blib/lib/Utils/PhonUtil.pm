package Utils::PhonUtil;

use 5.008006;
use strict;
use warnings;

require Exporter;
use AutoLoader qw(AUTOLOAD);

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use Utils::PhonUtil ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	trim
	vecdiff
	vecmag2
	vmatmult
	dotprod
	crossprod
	cellvol
	loadposcar
	loadharmonic
);

our $VERSION = '0.02';


# Preloaded methods go here.

sub trim {
	$_[0]=~s/^\s+//;
	$_[0]=~s/\s+$//;
	return $_[0];
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

sub cellvol(\@\@\@) {
	my $a1 = shift;
	my $a2 = shift;
	my $a3 = shift;
	my $vol = $a1->[0]*($a2->[1]*$a3->[2] - $a2->[2]*$a3->[1])
	        + $a1->[1]*($a2->[2]*$a3->[0] - $a2->[0]*$a3->[2])
	        + $a1->[2]*($a2->[0]*$a3->[1] - $a2->[1]*$a3->[0]);
	return abs($vol);
}

sub loadposcar(@) {
	my %latt;
	$latt{'desc'} = shift;
	chomp($latt{'desc'});
	$latt{'scale'} = &trim(shift);
	$latt{'a1'} = [shift =~ /([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)/g];
	$latt{'a2'} = [shift =~ /([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)/g];
	$latt{'a3'} = [shift =~ /([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)/g];
	$latt{'amat'} = [$latt{'a1'},$latt{'a2'},$latt{'a3'}];
	die "Lattice vectors are not 3 elements, stopped"
		unless ($#{$latt{'a1'}} == 2 && $#{$latt{'a2'}} == 2 && $#{$latt{'a3'}} == 2);
	$latt{'lattvol'} = cellvol(@{$latt{'a1'}},@{$latt{'a2'}},@{$latt{'a3'}});
	if($latt{'scale'} < 0) {
		my $tmpscale = (-$latt{'scale'}/$latt{'lattvol'})**(1/3);
		$latt{'scale'} = $tmpscale;
	}
	$latt{'natoms'} = [shift =~ /(\d+)/g];
	my $coordstr = shift;
	if($coordstr =~ /^\s*[Ss]/) { #Selective Dynamics
		$latt{'sd'} = 1;
		$coordstr = shift;
	} else {
		$latt{'sd'} = 0;
	}
	# Cartesian or Direct coordinates
	$latt{'coord'} = ($coordstr =~ /^\s*[CKck]/) ? 'C' : 'D';

	my $idx = 0;
	$latt{'atompos'} = [];
	$latt{'cartpos'} = [];
	my $atomsum = 0;
	for(my $atom = 0; $atom <= $#{$latt{'natoms'}}; $atom++) {
		my $atomposlist = [];
		my $cartposlist = [];
		$atomsum += ${$latt{'natoms'}}[$atom];
		while($idx<$atomsum) {
			my $atompos = [$_[$idx] =~ /([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)/g];
			push(@{$atomposlist}, $atompos);
			my $cartpos;
			if($latt{'coord'} eq 'D') {
				$cartpos = vmatmult($atompos, $latt{'amat'});
			} else {
				$cartpos = [$$atompos[0], $$atompos[1], $$atompos[2]];
			}
			$$cartpos[0] *= $latt{'scale'};
			$$cartpos[1] *= $latt{'scale'};
			$$cartpos[2] *= $latt{'scale'};
			push(@{$cartposlist}, $cartpos);
			$idx++;
		}
		push(@{$latt{'atompos'}}, $atomposlist);
		push(@{$latt{'cartpos'}}, $cartposlist);
	}
	return %latt;
}

sub loadharmonic(@) {
	my %harm;
	if(shift !~ /([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s+cutoff radius/) {
		die 'HARMONIC is not the correct format, stopped';
	}
	$harm{'cutoff'} = $1;
	if(shift !~ /(\d+)\s+number of vectors/) {
		die 'HARMONIC is not the correct format, stopped';
	}
	$harm{'nvec'} = $1;
	$harm{'vecs'} = [];
	$harm{'dmat'} = [];
	my $maxacount = 0;
	for(my $idx = 0; $idx < $harm{'nvec'}; $idx++) {
		my $s = $idx*6;
		my $row;
		($row->{'a1'}, $row->{'a2'}, $row->{'r'}, $row->{'ws'}) = ($_[$s] =~ /vector:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/);
		#next if($row->{'ws'} == 0);
		$row->{'R+t2-t1'} = [$_[$s+1] =~ /([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)/g];
		die 'HARMONIC is not the correct format, stopped' if $_[$s+2] !~ /force constant matrix:/;
		$row->{'dmat'} = [];
		push(@{$row->{'dmat'}}, [$_[$s+3] =~ /([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)/g]);
		push(@{$row->{'dmat'}}, [$_[$s+4] =~ /([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)/g]);
		push(@{$row->{'dmat'}}, [$_[$s+5] =~ /([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)/g]);
		${$harm{'dmat'}}[$row->{'r'}][$row->{'a1'}-1][$row->{'a2'}-1] = $row->{'dmat'};
		$harm{'rcount'}->{$row->{'r'}}++;
		if($row->{'a1'} > $maxacount) {
			$maxacount = $row->{'a1'};
		}
		if($row->{'a2'} > $maxacount) {
			$maxacount = $row->{'a2'};
		}
		if($row->{'a1'} == $row->{'a2'}) {
			$harm{'rhash'}->{$row->{'r'}} = $row->{'R+t2-t1'};
		}
		push(@{$harm{'vecs'}}, $row);
	}
	$harm{'t2-t1'} = [];
	while(my ($key, $val) = each(%{$harm{'rcount'}})) {
		if($val == $maxacount**2) {
			my $R = $harm{'rhash'}->{$key};
			foreach my $vec (@{$harm{'vecs'}}) {
				next unless $vec->{'r'} == $key;
				$harm{'t2-t1'}->[$vec->{'a1'}][$vec->{'a2'}][0] = $vec->{'R+t2-t1'}->[0] - $R->[0];
				$harm{'t2-t1'}->[$vec->{'a1'}][$vec->{'a2'}][1] = $vec->{'R+t2-t1'}->[1] - $R->[1];
				$harm{'t2-t1'}->[$vec->{'a1'}][$vec->{'a2'}][2] = $vec->{'R+t2-t1'}->[2] - $R->[2];
			}
			last;
		}
	}
	die 'Could not get all t2-t1 values, stopped' if($#{$harm{'t2-t1'}} == -1);
	foreach my $key (keys %{$harm{'rcount'}}) {
		next if(defined($harm{'rhash'}->{$key}));
		foreach my $vec (@{$harm{'vecs'}}) {
			next unless $vec->{'r'} == $key;
			$harm{'rhash'}->{$key}->[0] = $vec->{'R+t2-t1'}->[0] - $harm{'t2-t1'}->[$vec->{'a1'}][$vec->{'a2'}][0];
			$harm{'rhash'}->{$key}->[1] = $vec->{'R+t2-t1'}->[1] - $harm{'t2-t1'}->[$vec->{'a1'}][$vec->{'a2'}][1];
			$harm{'rhash'}->{$key}->[2] = $vec->{'R+t2-t1'}->[2] - $harm{'t2-t1'}->[$vec->{'a1'}][$vec->{'a2'}][2];
			last;
		}
	}
	return %harm;
}

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Utils::PhonUtil - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Utils::PhonUtil;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for Utils::PhonUtil, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Joseph Yasi, E<lt>yasi@uiuc.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2007 by Joseph Yasi

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.6 or,
at your option, any later version of Perl 5 you may have available.


=cut
