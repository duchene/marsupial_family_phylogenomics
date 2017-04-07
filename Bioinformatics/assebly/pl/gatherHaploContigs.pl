#!/usr/bin/perl -w

# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use warnings;
use strict;

my ($assemdir, $samples_list, $exonlist) = @ARGV;

my $gatherhaplocontigs_dir = "$assemdir/gatherHaplocontigs/";
unless(-e $gatherhaplocontigs_dir or mkdir $gatherhaplocontigs_dir) { die "[ERROR gatherhaplocontigs] Could not make $gatherhaplocontigs_dir\n"; }    


open LIBS, "<$samples_list" or die "[ERROR gatherHaplocontigs] Could not open the sample names file $samples_list\n";
my @libs  = <LIBS>;
close(LIBS);
chomp(@libs);

open EXONS, "<$exonlist" or die "[ERROR gatherHaplocontigs] Could not open the target exon IDs file $exonlist\n";
my @exons = <EXONS>;
close(EXONS);
chomp(@exons);

foreach my $exon_name (@exons) {
    
    # Output file for each target exon
    my $exon_best_h0_contigs    = "$gatherhaplocontigs_dir/${exon_name}_best_h0_contigs.fasta";
    my $exon_best_h1_contigs    = "$gatherhaplocontigs_dir/${exon_name}_best_h1_contigs.fasta";

    foreach my $lib (@libs) {

        # For each sample, extract best contig sequence and change sequence ID to sample name
        my $vcf2haplofasta_refs_h0_fil = "$assemdir/$lib/${lib}_vcf2haplofasta/${lib}_best2refs.vcf2haplofasta_refs_h0.fasta";
        my $vcf2haplofasta_refs_h1_fil = "$assemdir/$lib/${lib}_vcf2haplofasta/${lib}_best2refs.vcf2haplofasta_refs_h1.fasta";

        open H0_OUT_FIL, ">>$exon_best_h0_contigs", or die "could not open H0 file";
        open H1_OUT_FIL, ">>$exon_best_h1_contigs", or die "could not open H1 file";

        open H0_IN_FIL, "<$vcf2haplofasta_refs_h0_fil", or die "could not open H0 file";
        open H1_IN_FIL, "<$vcf2haplofasta_refs_h1_fil", or die "could not open H1 file";

        my @liblines_h0 = <H0_IN_FIL>;
        my $allliblines_h0 = join('',@liblines_h0);
        my @fasta_h0 = split(/\n>/,$allliblines_h0);

        foreach my $line (@fasta_h0) {
                 
             $line =~ m/^(\S+)_h0/;
             if ($exon_name eq $1) {

                  my $newline = $line;
                  $newline =~ s/$exon_name/$lib/;
                  print H0_OUT_FIL ">$newline\n";
             }
        }


        my @liblines_h1 = <H1_IN_FIL>;
        my $allliblines_h1 = join('',@liblines_h1);
        my @fasta_h1 = split(/\n>/,$allliblines_h1);

        foreach my $line (@fasta_h1) {
                 
             $line =~ m/^(\S+)_h1/;
             if ($exon_name eq $1) {
                  my $newline = $line;
                  $newline =~ s/$exon_name/$lib/;
                  print H1_OUT_FIL ">$newline\n";
             }
        }





    }
}


__END__

=head1 NAME

gatherHaplocontigs - For each target, from each sample, gather the best contig sequence that was amended with IUPAC codes in vcf2haplofasta.pl

=head1 USAGE

=over

=item B<perl gathercontigs.pl [ARGS]>

All arguments are required, in order.

=back

=head1 ARGUMENTS

=over

=item $assemdir

Output directory.

=item $samples_list

List of sample names.

=item $exonlist

Text file with target exon IDs.

=back

=head1 SUBROUTINES

None.

=head1 DIAGNOSTICS

=over 

=item [ERROR gathercontigs] Could not open/make ...

Required input/output file not openable. Check that previous pipeline stages have completed successfully, and the specified files exist. Check that you have adequate permissions in the output directory.

=back

=head1 DEPENDENCIES

None.

=head1 KNOWN BUGS

None.

=head1 NOTES

None.

=head1 AUTHORS

Ben Bai, ANU (u5205339@anu.edu.au)

Dr. Jason Bragg, ANU (jason.bragg@anu.edu.au)

=head1 SEE ALSO

exon-capture-phylo Design and Usage manual

=head1 LICENCE AND COPYRIGHT

Copyleft 2013-2014 by the authors.

This script is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 

=cut
