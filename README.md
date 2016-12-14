[![Build Status](https://travis-ci.org/wang-q/AlignDB-Codon.svg?branch=master)](https://travis-ci.org/wang-q/AlignDB-Codon) [![Coverage Status](http://codecov.io/github/wang-q/AlignDB-Codon/coverage.svg?branch=master)](https://codecov.io/github/wang-q/AlignDB-Codon?branch=master) [![MetaCPAN Release](https://badge.fury.io/pl/AlignDB-Codon.svg)](https://metacpan.org/release/AlignDB-Codon)
# NAME

AlignDB::Codon - translate sequences and calculate Dn/Ds

# DESCRIPTION

AlignDB::Codon provides methods to translate sequences and calculate Dn/Ds with different codon tables.

For more information, see [AlignDB](https://metacpan.org/pod/AlignDB).

# ATTRIBUTES

## one2three

lookup hash for one-letter aa names to three-letter ones, isa HashRef

## three2one

lookup hash for three-letter aa names to one-letter ones, isa HashRef

## codons

all codons, isa ArrayRef

## codon2aa

lookup hash for codons to aa, isa HashRef

## table\_id

codon table id, in Bio::Tools::CodonTable

## table\_name

codon table name, in Bio::Tools::CodonTable

## codon\_table

isa Bio::Tools::CodonTable Object

## syn\_sites

lookup hash for the number of synonymous changes per codon, isa HashRef

## syn\_changes

lookup hash of all pairwise combinations of codons differing by 1
1 = synonymous, 0 = non-synonymous, -1 = stop,
isa HashRef

# METHODS

## change\_codon\_table

         Usage : $obj->change_codon_table(2);
       Purpose : change used codon table and recalc all attributes
       Returns : none
    Parameters : a legal codon table id
        Throws : codon table id is not defined
               : or
               : codon table id should be in range of 1-6,9-16,21
      Comments : none
      See Also : n/a

## convert\_123

         Usage : my $three_format = $obj->convert_123('ARN');
       Purpose : convert aa code from one-letter to three-letter
       Returns : Str
    Parameters : IUPAC one-letter amino acid string
        Throws : Given characters not in IUPAC table
      Comments : none
      See Also : convert_321

## convert\_321

         Usage : my $one_format = $obj->convert_321('AlaArgAsn');
       Purpose : convert aa code from three-letter to one-letter
       Returns : Str
    Parameters : IUPAC three-letter amino acid string
        Throws : Given characters not in IUPAC table
      Comments : none
      See Also : convert_123

## comp\_codons

         Usage : my ($syn, $nsy) = $obj->comp_codons('TTT', 'GTA');
               : or
               : my ($syn1, $nsy1) = $obj->comp_codons('TTT', 'GTA', 1);
       Purpose : compares 2 codons to find the number of synonymous and
               :   non-synonymous mutations between them
               : if the third parameter is given, this method will
               :   return syn&nsy at this position
       Returns : (Num, Num)
    Parameters : Codon, Codon, Codon Position (optional, in 0 .. 2)
        Throws : Wrong codon
               : Wrong codon position
      Comments : none
      See Also : n/a

## is\_start\_codon

         Usage : my $bool = $obj->is_start_codon('ATG')
       Purpose : returns true (1) for codons that can be used as a
               :   translation start, false (0) for others.
       Returns : boolean
    Parameters : Codon
        Throws : no exceptions
      Comments : none
      See Also : n/a

## is\_ter\_codon

         Usage : my $bool = $obj->is_ter_codon('GAA')
       Purpose : returns true (1) for codons that can be used as a
               :   translation terminator, false (0) for others.
       Returns : boolean
    Parameters : Codon
        Throws : no exceptions
      Comments : none
      See Also : n/a

## is\_unknown\_codon

         Usage : my $bool = $obj->is_unknown_codon('GAJ')
       Purpose : returns true (1) for codons that are valid,
               :   true (1) for others.
       Returns : boolean
    Parameters : Codon
        Throws : no exceptions
      Comments : none
      See Also : n/a

# AUTHOR

Qiang Wang &lt;wang-q@outlook.com>

# COPYRIGHT AND LICENSE

This software is copyright (c) 2008 by Qiang Wang.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.
