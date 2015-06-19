#!/usr/bin/env perl

use strict;
use warnings;

my $text = <>;
chomp($text);
my $numstr = <>;
$numstr =~ s/^\s+//;
my @numarr = split(/\s+/,$numstr);
my $numrvecs = $numarr[0];
my $numatoms = $numarr[1];

my @rvecs = ();
my @dmat2d = ();
my @rvec2dmap = ();
my @rvecs2d = ();

while(my $line = <>) {
	$line =~ s/^\s+//;
	chomp($line);
	my @dline = split(/\s+/, $line);
	my $rv2d = [$dline[0],$dline[1],$dline[2]];
	my $dm2d = [];
	for(my $jdx=3; $jdx<3+9*$numatoms**2; $jdx++) {
		push(@{$dm2d}, $dline[$jdx]);
	}
	push(@rvecs,$rv2d);
	my $vecfound = 0;
	my $idx = 0;
	foreach my $rvec2d (@rvecs2d) {
		my $vecdiff = ($$rvec2d[0]-$$rv2d[0])**2 + ($$rvec2d[1]-$$rv2d[1])**2;
		if($vecdiff < 1.0e-7) {
			$vecfound = 1;
			push(@rvec2dmap, $idx);
			for(my $jdx=0;$jdx<9*$numatoms**2;$jdx++) {
				${$dmat2d[$idx]}[$jdx] += $$dm2d[$jdx];
			}
			if(abs($$rvec2d[2]) > abs($$rv2d[2])) {
				$rvecs2d[$idx] = $rv2d;
			}
			last;
		}
		$idx++;
	}
	if($vecfound == 0) {
		push(@rvecs2d, $rv2d);
		push(@dmat2d, $dm2d);
		push(@rvec2dmap, $idx);
	}
}

print "$text 2d\n";
print(scalar @rvecs2d . " $numatoms\n");
for(my $idx=0; $idx<scalar @rvecs2d; $idx++) {
	printf("%+.12f %+.12f %+.12f", ${$rvecs2d[$idx]}[0], ${$rvecs2d[$idx]}[1], ${$rvecs2d[$idx]}[2]);  
	for(my $jdx=0;$jdx<9*$numatoms**2;$jdx++) {
		printf(" %+.12e", ${$dmat2d[$idx]}[$jdx]);
	}
	print "\n";
}
