#!/usr/bin/env perl

use strict;
use warnings;
use Utils::PhonUtil;

open(INHARM, '<HARMONIC') or die 'Could not open HARMONIC, stopped';
my %harm = loadharmonic(<INHARM>);
close(INHARM);
