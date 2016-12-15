use strict;
use warnings;

use Test::More;

use AlignDB::Codon;
use Bio::Align::DNAStatistics;

{    # syn_sites
    my $codon_obj = AlignDB::Codon->new( table_id => 1 );
    my $comp_obj = Bio::Align::DNAStatistics->new;

    my $result1 = $codon_obj->syn_sites;
    my $result2 = $comp_obj->get_syn_sites;

    is_deeply( $result1, $result2, "syn_sites" );
}

{    # syn_changes
    my $codon_obj = AlignDB::Codon->new( table_id => 1 );
    my $comp_obj = Bio::Align::DNAStatistics->new;

    my $result1 = $codon_obj->syn_changes;
    my %result2 = $comp_obj->get_syn_changes;

    is_deeply( $result1, \%result2, "syn_changes" );
}

done_testing();
