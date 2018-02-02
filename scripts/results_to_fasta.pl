#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

cat results.txt | results_to_fasta directory

Description:

This script parses 'Matlock cluster' results. The output directory should
be empty, the script just appends to files if they exist.

";


my ($help);
my $opt_success = GetOptions('help'    => \$help,
			      );

die $usage if $help || ! $opt_success;


my $directory = shift;

LINE: while (<STDIN>) {

    next LINE if($_ =~ /\#/);

    my @l = split /\t/, $_;

    open (my $IN, '>>', "$directory/$l[3].cluster.fasta") or die "Can't open $l[3].cluster.fasta for reading\n$!\n";

    print $IN ">$l[0] $l[1] $l[2] $l[3]\n$l[4]";

    close($IN);

}

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

sub {

}
