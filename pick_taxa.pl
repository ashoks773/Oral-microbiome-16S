#!usr/bin/perl
use strict;

use Data::Dumper;
my $hash={};

my $file = $ARGV[0];
open (IF, $file);
while (chomp (my $line=<IF>))
{
	my @arr = split ("\t", $line);
	$hash->{$arr[0]}=$arr[1];
}
#print Dumper $hash;

my $file1 = $ARGV[1];
open (IF1, $file1);
open (OF, ">$file1.taxonomy.txt");
while (chomp (my $line1=<IF1>))
{
        my @arr1 = split ("\t", $line1);
	if (exists ($hash->{$arr1[0]}))
	{
		print OF "$hash->{$arr1[0]}.\t$line1\n";
	}
}
