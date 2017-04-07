#!/usr/bin/perl -w

# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use warnings;
use strict;
use Bio::SeqIO;

my ($lib, $assemdir) = @ARGV;

# Open final GATK output
my $vcf = "$assemdir/$lib/${lib}_gatkSNPcalls/$lib.ReadGrouped.ReadBackedPhased.VariantAnnotated.VariantFiltered.vcf";
#my $vcf = "test.vcf";
open(VCF, "<$vcf") or die "[ERROR vcf2haplofasta $lib] Could not open .vcf input file $vcf\n";
my @vcf = <VCF>;
close VCF;

my $vcf2haplofasta_dir = "$assemdir/$lib/${lib}_vcf2haplofasta/";
unless(-e $vcf2haplofasta_dir or mkdir $vcf2haplofasta_dir) { die "[ERROR vcf2haplofasta $lib] Could not make $vcf2haplofasta_dir\n"; }

my $reffil = "$assemdir/$lib/${lib}_best2refs/${lib}_best2refs.fasta";
#my $reffil = "references.fa";
my $ref_IN = Bio::SeqIO->new(-file => "<$reffil",
                             -format => "fasta",
                             -alphabet => "dna");

my $output_reffil_h0 = "$vcf2haplofasta_dir/${lib}_best2refs.vcf2haplofasta_refs_h0.fasta";
#my $output_reffil = "out.fasta";
my $ref_OUT_h0 = Bio::SeqIO->new(-file => ">$output_reffil_h0",
                              -format => "fasta",
                              -alphabet => "dna");

my $output_reffil_h1 = "$vcf2haplofasta_dir/${lib}_best2refs.vcf2haplofasta_refs_h1.fasta";
my $ref_OUT_h1 = Bio::SeqIO->new(-file => ">$output_reffil_h1",
                              -format => "fasta",
                              -alphabet => "dna");



# Tally types of variants
my $output_statsfil = "$vcf2haplofasta_dir/${lib}_best2refs.vcf2haplofasta_refs.stats";
#my $output_statsfil = "out.stat";
open(STATS_OUT, ">$output_statsfil"); #or die "[ERROR vcf2haplofasta $lib] Could not open output reference file $output_statsfil\n";
print STATS_OUT "REF\tHET\tALT_HOMO\tINDEL\tSEQ_LEN\tPHASE\tBREAKS\n";


# Loop through reference sequences
while (my $ref = $ref_IN->next_seq) {
    my $ref_name = $ref->display_id();

    my $h0_name    = $ref_name . "_h0";
    my $h1_name    = $ref_name . "_h1";
    my $ref_seq_h0 = $ref->seq();
    my $ref_seq_h1 = $ref->seq();

    # Whilst tallying detected variant types
    my $hets  = 0;
    my $alts  = 0;
    my $indel = 0;
    my $other = 0;
    my $phase = 0;

 
    foreach my $vcline (@vcf) {
        # a0 and a1 are the reference and alternate bases
        my ($vcf_chrom_name, $snp_pos, undef, $a0, $a1, undef, $filter_result, undef, undef, $phasing_info) = split(/\t/, $vcline);
#print "$vcf_chrom_name\t$ref_name\t$filter_result\n";
        # Check the variant passed filtering and characterise it based on phasing information
        if ($vcf_chrom_name eq $ref_name && $filter_result eq "PASS") {
#print "passed\n";
            # if indel
            if (length($a0) > 1 || length($a1) > 1) {
                $indel++;
            # if alternate homozygote, replace the reference with the alternate
            } elsif ($phasing_info =~ /^1[\|\/]1/) { 
                substr($ref_seq_h0, $snp_pos-1, 1) = $a1;
                substr($ref_seq_h1, $snp_pos-1, 1) = $a1;  
                $alts++;

            # if heterozygote, assign alleles to haplotypes
            } elsif ($phasing_info =~ /^0[\|\/]1/) {
                substr($ref_seq_h0, $snp_pos-1, 1) = $a0;
                substr($ref_seq_h1, $snp_pos-1, 1) = $a1;  
                if ($phasing_info =~ /\|/) {$phase++;}
                $hets++;

            } elsif ($phasing_info =~ /^1[\|\/]0/) {  
                substr($ref_seq_h0, $snp_pos-1, 1) = $a1;
                substr($ref_seq_h1, $snp_pos-1, 1) = $a0;  
                if ($phasing_info =~ /\|/) {$phase++;}
                $hets++;

            } else {
                $other++;
            }
        }
    }

#print "$ref_seq_h0\t$ref_seq_h1\n";

    my $ref_h0 = $ref;
    $ref_h0->id("$h0_name");
    $ref_h0->seq("$ref_seq_h0");
    my $seqlen_h0 = $ref_h0->length();
    $ref_OUT_h0->write_seq($ref_h0);

    my $ref_h1 = $ref;
    $ref_h1->id("$h1_name");
    $ref_h1->seq("$ref_seq_h1");
    my $seqlen_h1 = $ref_h1->length();
    $ref_OUT_h1->write_seq($ref_h1);

#    $ref_OUT->width($seqlen_h0);
#    $ref_OUT->width($seqlen_h1);

    my $breaks = $hets - $phase - 1;
    print STATS_OUT "$ref_name\t$hets\t$alts\t$indel\t$seqlen_h0\t$phase\t$breaks\n";
}

close STATS_OUT;

__END__

=head1 NAME

vcf2haplofasta - Write SNPs to best contigs reference.

=head1 USAGE

=over

=item B<perl vcf2haplofasta.pl [ARGS]>

All arguments are required, in order.

=back

=head1 ARGUMENTS

=over

=item $lib

Sample name.

=item $assemdir

Output directory.

=back

=head1 SUBROUTINES

None.

=head1 DIAGNOSTICS

=over 

=item [ERROR vcf2haplofasta $lib] Could not open/make ...

Required input/output file not openable. Check that previous pipeline stages have completed successfully, and the specified files exist. Check that you have adequate permissions in the output directory.

=item [ERROR vcf2haplofasta $lib] Invalid base detected ...

Base combination does not have a corresponding IUPAC code. Check for .vcf file corruption.

=back

=head1 DEPENDENCIES

BioPerl

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
