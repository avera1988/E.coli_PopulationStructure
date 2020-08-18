#!/usr/bin/perl

while(<>){
	chomp;
	if($_ ne "--"){
		print "$_\n";
	}
}
