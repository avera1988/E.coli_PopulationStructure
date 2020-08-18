#!/usr/bin/bash
#############################################################
#	This scritp pars the faa clusters and look for the nt sequences in a concatenated file all.ecoli.fasta
#	Dpendencies c_header_for_nt.pl, removing.hyphen.pl all.ecoli.fasta (index file with all genome ffn seq)
#Author Arturo Vera
############################################################

for i in *.faa; do
	name=$(basename $i .faa);
	perl c_header_for_nt.pl $i|\
	awk '{if($0 ~/^>/) print $0}'|\
	sed 's/>//'|\
	fgrep -w -f - -A 1 \
	all.ecoli.fasta |\
	perl removing.hyphen.pl > $name.ffn;
done

for i in *.ffn ;
	do
	perl -e 'while(<>){chomp;if($_ =~/^>/){@h=split(/\_/);}else{$seq=$_;print ">$h[2]_$h[3]\n$seq\n";}}' $i >$i.mod;
done
