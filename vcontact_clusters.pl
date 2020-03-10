```perl
#!/usr/local/bin/perl
# Retrieves a list of contigs and clusters from vCONTACT:
# awk -F',' '{print $5 "\t" $14}' ../cytophage_METADATA3.csv | grep -v '_pro' | sort > contigClusters.txt
# to recover FASTA files per cluster (assuming that is a FASTA single-line)
#usage: perl vcontact_clusters.pl <contig_list> <fasta_file> <outfile_prefix>

use strict 'vars';
use List::MoreUtils qw(first_index);
my $infile = @ARGV[0];
my $fafile = @ARGV[1];
my $oufile = @ARGV[2];

my %hashr = ();

# Read list file and make a hash
open (INFILE, "$infile")   || die "can't open file $infile\n";
while (<INFILE>) {
	my @lyne = split(/\t/,$_);
	$hashr{$lyne[1]} = $lyne[0];
}
close INFILE;

# Retrieve sequences
open (FASTAFILE, "$fafile") || die "can't open file $fafile\n";
chomp (my @fasta = <FASTAFILE>);
close FASTAFILE;

my $index;
my $m;
my $contig;
my $cluster;
my @notfound;
my $scalar;

open (OUTFILE,">> $oufile")   || die "can't open file $oufile\n";

while (my ($contig,$cluster) = each %hashr) {
	#print "$contig is $cluster \n";
	$contig =~ s/VIRSorter_//;
	$contig =~ s/_len.*/_/;
	chomp $contig;
	print ("grep '$contig' $fafile -A1 >> $cluster \n");
        system("grep '$contig' $fafile -A1 >> $cluster");
  print ("grep '$contig' ~/Documents/phages/VirSorter/infiles/virome_source.fasta -A1 >> $cluster \n");
  system("grep '$contig' ~/Documents/phages/VirSorter/infiles/virome_source.fasta -A1 >> $cluster");
	my $index = first_index {/$contig/} @fasta;
	my ($index) = grep { $fasta[$_] eq $contig } (0 .. @fasta-1);
	if ($index == -1) {
		$scalar = $cluster."\t".$contig;
		push (@notfound,$scalar);
		undef $index;
		next
	} else {
		$m = $index - 1;
		chomp $contig;
		print OUTFILE "\n>".$cluster."\t".$contig."\n".$fasta[$m];
		undef $index;
	}
}
close OUTFILE;
while (@notfound) {
	print "$_"
}
```
