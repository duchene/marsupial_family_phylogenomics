#!/usr/bin/perl -w

# Documentation @ __END__
# WARNING: not designed for use as a standalone module.

use warnings;
use strict;

my ($assemdir, $samples_list, $exonlist) = @ARGV;


open LIBS, "<$samples_list" or die "[ERROR gathercontigs] Could not open the sample names file $samples_list\n";
open EXONS, "<$exonlist" or die "[ERROR gathercontigs] Could not open the target exon IDs file $exonlist\n";
my @libs  = <LIBS>;
my @exons = <EXONS>;
chomp(@libs); chomp(@exons);
close(LIBS); close(EXONS);

foreach my $lib (@libs) {
    my $assemlib = "$assemdir/$lib/";
    my $contigcountfile = "$assemlib/$lib" . "_contigcount.txt";
    open COUNT, ">$contigcountfile" or die "[ERROR countcontigs] Could not open the contig count file $contigcountfile\n";

    my %contignum=();
    foreach my $exon_name (@exons) {

        if ($exon_name =~ /^(\S+?)_\S+/) {

            # Non-whitespace portion of exon ID before underscore should be the orthologous protein ID
            my $prot = $1;

            my $bestcontig_distrib_dir = "$assemlib/$prot/${prot}_bestcontig_distrib/";

            my $contigfilall = "$bestcontig_distrib_dir/${exon_name}_velvet_contigs.cap3ed.exonerated.filtered.fasta";

            # Collate all contigs for the target that passed filtering
            if (-e $contigfilall) {

                open CTGA, "<$contigfilall" or die "[ERROR gathercontigs] Could not open the filtered contigs file $contigfilall\n";
                my @confilalllines = <CTGA>;

                # Eliminate multiline sequences
                my $alltmpstring = join('', @confilalllines);
                my @tmplines = split(/\n>/, $alltmpstring);

                # Count contigs
                my $tmpnum = ( scalar @tmplines );
                $contignum{$exon_name} = $tmpnum;

                close(CTGA);
            
            } else {
                $contignum{$exon_name} = 0;
            }

            print COUNT $exon_name ."\t". $contignum{$exon_name} ."\n";
        }

    }
    close(COUNT);  
}




__END__

=head1 NAME

gathercontigs - For each target, from each sample, gather all contigs that passed filtering, and the best contig. Also count numbers of filtered contigs.

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

The regex used to parse target exon names restricts the sequence ID format:
    if ($exon_name =~ /^(\S+?)_\S+/) 

=head1 AUTHORS

Ben Bai, ANU (u5205339@anu.edu.au)

Dr. Jason Bragg, ANU (jason.bragg@anu.edu.au)

=head1 SEE ALSO

exon-capture-phylo Design and Usage manual

=head1 LICENCE AND COPYRIGHT

Copyleft 2013-2014 by the authors.

This script is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 

=cut
