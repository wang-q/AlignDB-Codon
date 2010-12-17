#!/usr/bin/perl
use strict;
use warnings;

use Test::More qw(no_plan);

BEGIN {
    use_ok('AlignDB::Codon');
}

# compare normal codons
{
    my $codon_obj = AlignDB::Codon->new( table_id => 1 );

    #    my $seq1 = <<ENDLINE;
    #GTCTGTTCCAAGGGCCTTTGCGTCAGG-TGGGC-TCAGGGTT---------------CCAGGGTGGCTGG
    #ACCCCAGGCCCCAGCTCTGCAGCAGGGAGGACGTGGCTGGGCTCGTGAAGCATGTGGGGGTGAGCCCAGG
    #GGCCCCAAGGCAGGGCACCTGGCCTTCAGCCTGCCTCAGCCCTGCCTGTCTCCCA
    #ENDLINE
    #    my $seq2 = <<ENDLINE ;
    #GTCTGTTCCAAGGGCCTTCGAGCCAGTCTGGGCCCCAGGGCTGCCCCACTCGGGGTTCCAGAGCAGTTGG
    #ACCCCAGGTCTCAGC---------GGGAGGGTGTGGCTGGGCTC-TGAAGCATTT--GGGTGAGCCCAGG
    #GGCTC-AGGGCAGGGCACCTG-CCTTCAGC-GGCCTCAGC-CTGCCTGTCTCCCA
    #ENDLINE

    my @compare = (
        [qw{ TTT TTA }], [qw{ TTT TTC }], [qw{ TTT GTA }], [qw{ TTT GTG }],
        [qw{ TTG AGA }],
    );

    my @expect
        = ( [ 0, 1 ], [ 1, 0 ], [ 0.5, 1.5 ], [ 0.5, 1.5 ], [ 0.75, 2.25 ], );

    for my $i ( 0 .. $#compare ) {
        my ( $cod1, $cod2 ) = @{ $compare[$i] };
        my ( $exp1, $exp2 ) = @{ $expect[$i] };
        my ( $syn,  $nsy )  = $codon_obj->comp_codons( $cod1, $cod2 );
        is( $syn, $exp1, "syn|$i" );
        is( $nsy, $exp2, "nsy|$i" );
    }
}

# compare normal codons
{
    my $codon_obj = AlignDB::Codon->new( table_id => 1 );

    my @compare = (
        [qw{ TTT TTA 2 }], [qw{ TTT TTA 0 }],
        [qw{ TTT TTC 2 }], [qw{ TTT TTC 1 }],
        [qw{ TTT GTA 0 }], [qw{ TTT GTA 1 }],
        [qw{ TTT GTA 2 }], [qw{ TTG AGA 0 }],
        [qw{ TTG AGA 1 }], [qw{ TTG AGA 2 }],
    );

    my @expect = (
        [ 0,    1 ],
        [ 0,    0 ],
        [ 1,    0 ],
        [ 0,    0 ],
        [ 0,    1 ],
        [ 0,    0 ],
        [ 0.5,  0.5 ],
        [ 0,    1 ],
        [ 0,    1 ],
        [ 0.75, 0.25 ],

    );

    for my $i ( 0 .. $#compare ) {
        my ( $cod1, $cod2, $pos ) = @{ $compare[$i] };
        my ( $exp1, $exp2 ) = @{ $expect[$i] };
        my ( $syn, $nsy ) = $codon_obj->comp_codons( $cod1, $cod2, $pos );
        is( $syn, $exp1, "syn_pos|$i" );
        is( $nsy, $exp2, "nsy_pos|$i" );
    }

}
