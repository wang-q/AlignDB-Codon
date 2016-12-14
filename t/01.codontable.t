#!/usr/bin/perl
use strict;
use warnings;

use Test::More tests => 66;

BEGIN {
    use_ok('AlignDB::Codon');
}

# id <=> name
{

    # all translation tables in Bio::Tools::CodonTable
    my @NAMES =    #id
        (
        'Standard',                                                                      # 1
        'Vertebrate Mitochondrial',                                                      # 2
        'Yeast Mitochondrial',                                                           # 3
        'Mold, Protozoan, and Coelenterate Mitochondrial and Mycoplasma/Spiroplasma',    # 4
        'Invertebrate Mitochondrial',                                                    # 5
        'Ciliate, Dasycladacean and Hexamita Nuclear',                                   # 6
        '', '',
        'Echinoderm and Flatworm Mitochondrial',                                         # 9
        'Euplotid Nuclear',                                                              # 10
        'Bacterial, Archaeal and Plant Plastid',                                         # 11
        'Alternative Yeast Nuclear',                                                     # 12
        'Ascidian Mitochondrial',                                                        # 13
        'Alternative Flatworm Mitochondrial',                                            # 14
        'Blepharisma Nuclear',                                                           # 15
        'Chlorophycean Mitochondrial',                                                   # 16
        '', '', '', '',
        'Trematode Mitochondrial',    # 21
                                      #'Scenedesmus obliquus Mitochondrial',             # 22
                                      #'Thraustochytrium Mitochondrial'                  # 23
        );

    for ( 0 .. $#NAMES ) {
        my $table_id   = $_ + 1;
        my $table_name = $NAMES[$_];
        next if length $table_name < 1;
        my $codon_obj = AlignDB::Codon->new( table_id => $table_id );
        ok( defined $codon_obj,                "Init object $table_id" );
        ok( $codon_obj->isa('AlignDB::Codon'), "ISA $table_id" );
        is( $codon_obj->table_id(),   $table_id,   "table_id $table_id" );
        is( $codon_obj->table_name(), $table_name, "table_name $table_id" );
    }

    print "\n";
}

# wrong codon table id
{

    # not a number
    {
        my @ids = qw{ a b c };
        for my $table_id (@ids) {
            eval { my $codon_obj = AlignDB::Codon->new( table_id => $table_id ); };
            ok( $@ =~ /Int/, "not a number" );
        }
    }

    # out of range
    {
        my @ids = qw{ 55 100 };
        for my $table_id (@ids) {
            eval { my $codon_obj = AlignDB::Codon->new( table_id => $table_id ); };
            ok( $@ =~ /range/, "out of range" );
        }
    }
}
