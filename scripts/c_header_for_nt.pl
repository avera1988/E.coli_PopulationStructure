#!/usr/bin/perl
#
#This script change the header of a geneMark translated ang get_homologues cluster and gives a header like: 
#	>gene_1_genome_1
#
#Author Arturo Vera
# avera@ccg.unam.mx
#############################################
use warnings;
use strict;
my(@header,@t,$headerlast,$seq);

while(<>){
	chomp;
	if($_ =~ /^>/){
	chomp;
	@header=split(/\t/); 
	foreach($header[0]){
		@t=split(/\|/);
	}
	$headerlast = $header[@header -1];
	$headerlast =~ s/\[//g;
	$headerlast =~ s/\]//g;
	}else{
	$seq=$_;
	print "$t[0]_genome_$headerlast\n$seq\n"
	}
}
