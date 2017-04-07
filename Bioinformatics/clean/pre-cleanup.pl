
use warnings;
use strict;

my @original_files = </moritz/jgb/readsFeb2016/LN01/all_raw_reads/*.fastq.gz>;         ## all original fastq sequence files should be saved in one folder;

my $Result_dir1 = "/moritz/jgb/readsFeb2016/LN01/combined/"; ## save combined fastq files in this folder (archival copy);
my $Result_dir2 = "/moritz/jgb/readsFeb2016/LN01/preclean/"; ## pre-clean the reads and save them in this folder;

foreach (<@original_files>) {
	my $ori_file = $_ =~ /(LN01_indexing\S+_R[1|2])_\d+.fastq.gz/; 
	my $lib_name = $1;
	my $file2 = "/moritz/jgb/readsFeb2016/LN01/all_raw_reads/" . $1 . '*' . 'fastq.gz';
	my $combined = $Result_dir1 . $1 . '.fastq.gz'; 
	system ("cat $file2 > $combined ") unless -e $combined;
	print "merging", $combined, "!", "\n";
	}

my @merged_files = < $Result_dir1*.fastq.gz> ;
foreach my $file (<@merged_files>) {
	my $out = $Result_dir2 .  $1 . $2 . ".fq" if $file  =~ /(LN01_indexing\d+)_\S+(_R[1|2]).fastq.gz/; 
	my $redundancy = '^--$' ;
	print "cleaning","\t",$file,"\n";
		if ($file =~ m/R(\d+)/) {
			if ($1 == 1) {
				system ("zcat $file | grep -A 3 \'^@.* [^:]*:N:[^:]*:\' | grep -v $redundancy | sed \'s/ 1:N:0:.*/\\/1/g\' >> $out");
			}
			if ($1 == 2) {
				system ("zcat $file | grep -A 3 \'^@.* [^:]*:N:[^:]*:\' | grep -v $redundancy | sed \'s/ 2:N:0:.*/\\/2/g\' >> $out");
			}

		}
}


