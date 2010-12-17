package AlignDB::Codon;
# ABSTRACT: translate sequences and calculate Dn/Ds

use Moose;
use Carp;

use Bio::Tools::CodonTable;
use List::MoreUtils qw(none);

use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;

=attr one2three

lookup hash for one-letter aa names to three-letter ones, isa HashRef

=cut

has 'one2three' => ( is => 'ro', isa => 'HashRef', );

=attr three2one

lookup hash for three-letter aa names to one-letter ones, isa HashRef

=cut

has 'three2one' => ( is => 'ro', isa => 'HashRef', );

=attr codons

all codons, isa ArrayRef

=cut

has 'codons' => ( is => 'ro', isa => 'ArrayRef', );

=attr codon2aa

lookup hash for codons to aa, isa HashRef

=cut

has 'codon2aa' => ( is => 'ro', isa => 'HashRef', );

=attr table_id

codon table id, in Bio::Tools::CodonTable

=cut

has 'table_id' => ( is => 'ro', isa => 'Int', );

=attr table_name

codon table name, in Bio::Tools::CodonTable

=cut

has 'table_name' => ( is => 'ro', isa => 'Str', );

=attr codon_table

isa Bio::Tools::CodonTable Object

=cut

has 'codon_table' => ( is => 'ro', isa => 'Object', );

=attr syn_sites

lookup hash for the number of synonymous changes per codon, isa HashRef

=cut

has 'syn_sites' => ( is => 'ro', isa => 'HashRef', );

=attr syn_changes

lookup hash of all pairwise combinations of codons differing by 1
1 = synonymous, 0 = non-synonymous, -1 = stop,
isa HashRef
 
=cut

has 'syn_changes' => ( is => 'ro', isa => 'HashRef', );

sub BUILD {
    my $self = shift;

    # Load aa code
    my ( $one2three, $three2one ) = $self->_load_aa_code();
    $self->{one2three} = $one2three;
    $self->{three2one} = $three2one;

    my @codons = $self->_make_codons();
    $self->{codons} = \@codons;

    my $table_id = $self->table_id() ? $self->table_id() : 1;
    $self->change_codon_table($table_id);

    return;
}

sub _load_aa_code {
    my $self = shift;

    my %one2three = (
        A   => 'Ala',    # Alanine
        R   => 'Arg',    # Arginine
        N   => 'Asn',    # Asparagine
        D   => 'Asp',    # Aspartic acid
        C   => 'Cys',    # Cysteine
        Q   => 'Gln',    # Glutamine
        E   => 'Glu',    # Glutamic acid
        G   => 'Gly',    # Glycine
        H   => 'His',    # Histidine
        I   => 'Ile',    # Isoleucine
        L   => 'Leu',    # Leucine
        K   => 'Lys',    # Lysine
        M   => 'Met',    # Methionine
        F   => 'Phe',    # Phenylalanine
        P   => 'Pro',    # Proline
        S   => 'Ser',    # Serine
        T   => 'Thr',    # Threonine
        W   => 'Trp',    # Tryptophan
        Y   => 'Tyr',    # Tyrosine
        V   => 'Val',    # Valine
        B   => 'Asx',    # Aspartic acid or Asparagine
        Z   => 'Glx',    # Glutamine or Glutamic acid
        X   => 'Xaa',    # Any or unknown amino acid
        '*' => '***',    # Stop codon
    );
    my %three2one = reverse %one2three;

    return ( \%one2three, \%three2one );
}

=method change_codon_table

      Usage : $obj->change_codon_table(2);
    Purpose : change used codon table and recalc all attributes
    Returns : none
 Parameters : a legal codon table id
     Throws : codon table id is not defined
            : or
            : codon table id should be in range of 1-6,9-16,21
   Comments : none
   See Also : n/a

=cut

sub change_codon_table {
    my $self = shift;
    my $id   = shift;

    # all codon table ids in Bio::Tools::CodonTable,
    # except the following two
    # 	 'Scenedesmus obliquus Mitochondrial', #22
    #    'Thraustochytrium Mitochondrial' #23
    my $id_set = AlignDB::IntSpan->new("1-6,9-16,21");

    if ( not defined $id ) {
        croak "codon table id is not defined\n";
    }
    elsif ( $id_set->contain($id) ) {
        my $codon_table = Bio::Tools::CodonTable->new( -id => $id );

        $self->{table_id}    = $id;
        $self->{codon_table} = $codon_table;
        $self->{table_name}  = $codon_table->name();

        $self->{codon2aa}    = $self->_codon2aa();
        $self->{syn_sites}   = $self->_syn_sites();
        $self->{syn_changes} = $self->_syn_changes();
    }
    else {
        croak "codon table id should be in range of $id_set\n";
    }

    return;
}

sub _codon2aa {
    my $self = shift;

    my $codons      = $self->codons();
    my $codon_table = $self->codon_table();

    my %codon2aa;
    foreach my $codon (@$codons) {
        my $aa = $codon_table->translate_strict($codon);
        $codon2aa{$codon} = $aa;
    }

    return \%codon2aa;
}

sub _make_codons {
    my $self = shift;

    # makes all codon combinations, returns array of them
    my @nucs = qw(T C A G);
    my @codons;
    for my $i (@nucs) {
        for my $j (@nucs) {
            for my $k (@nucs) {
                push @codons, "$i$j$k";
            }
        }
    }

    return @codons;
}

=method convert_123

      Usage : my $three_format = $obj->convert_123('ARN');
    Purpose : convert aa code from one-letter to three-letter
    Returns : Str
 Parameters : IUPAC one-letter amino acid string
     Throws : Given characters not in IUPAC table
   Comments : none
   See Also : convert_321

=cut

sub convert_123 {
    my $self    = shift;
    my $peptide = shift;

    $peptide = uc $peptide;
    my $three_of = $self->one2three();

    my $converted;
    for my $pos ( 0 .. length($peptide) - 1 ) {
        my $aa_code = substr( $peptide, $pos, 1 );
        if ( $three_of->{$aa_code} ) {
            $converted .= $three_of->{$aa_code};
        }
        else {
            carp "Wrong single-letter amino acid code [$aa_code]!\n";
            $converted .= ' ' x 3;
        }
    }
    return $converted;
}

=method convert_321

      Usage : my $one_format = $obj->convert_321('AlaArgAsn');
    Purpose : convert aa code from three-letter to one-letter
    Returns : Str
 Parameters : IUPAC three-letter amino acid string
     Throws : Given characters not in IUPAC table
   Comments : none
   See Also : convert_123

=cut

sub convert_321 {
    my $self    = shift;
    my $peptide = shift;

    $peptide = lc $peptide;
    my $one_of = $self->three2one();

    my $converted;
    for ( my $pos = 0; $pos < length($peptide); $pos += 3 ) {
        my $aa_code = substr( $peptide, $pos, 3 );
        $aa_code = ucfirst $aa_code;
        if ( $one_of->{$aa_code} ) {
            $converted .= $one_of->{$aa_code};
        }
        else {
            carp "Wrong three-letter amino acid code [$aa_code]!\n";
            $converted .= ' ' x 3;
        }

    }

    return $converted;
}

sub _syn_changes {
    my $self = shift;

    my $codons   = $self->codons();
    my $codon2aa = $self->codon2aa();

    my $arr_len = scalar @$codons;

    my %results;
    for ( my $i = 0; $i < $arr_len - 1; $i++ ) {
        my $cod1 = $codons->[$i];
        for ( my $j = $i + 1; $j < $arr_len; $j++ ) {
            my $cod2     = $codons->[$j];
            my $diff_cnt = 0;
            for my $pos ( 0 .. 2 ) {
                if ( substr( $cod1, $pos, 1 ) ne substr( $cod2, $pos, 1 ) ) {
                    $diff_cnt++;
                }
            }
            next if $diff_cnt != 1;

            # synonymous change
            if ( $codon2aa->{$cod1} eq $codon2aa->{$cod2} ) {
                $results{$cod1}{$cod2} = 1;
                $results{$cod2}{$cod1} = 1;
            }

            # stop codon
            elsif ( $codon2aa->{$cod1} eq '*' or $codon2aa->{$cod2} eq '*' ) {
                $results{$cod1}{$cod2} = -1;
                $results{$cod2}{$cod1} = -1;
            }

            # non-synonymous change
            else {
                $results{$cod1}{$cod2} = 0;
                $results{$cod2}{$cod1} = 0;
            }
        }
    }

    return \%results;
}

sub _syn_sites {
    my $self = shift;

    my $codons   = $self->codons();
    my $codon2aa = $self->codon2aa();

    my %raw_results;
    for my $cod (@$codons) {
        my $aa = $codon2aa->{$cod};

        # calculate number of synonymous mutations vs non-syn mutations
        for my $i ( 0 .. 2 ) {
            my $s = 0;
            my $n = 3;
            for my $nuc (qw(A T C G)) {
                next if substr( $cod, $i, 1 ) eq $nuc;
                my $test = $cod;
                substr( $test, $i, 1, $nuc );
                if ( $codon2aa->{$test} eq $aa ) {
                    $s++;
                }
                if ( $codon2aa->{$test} eq '*' ) {
                    $n--;
                }
            }
            $raw_results{$cod}[$i] = {
                's' => $s,
                'n' => $n
            };
        }
    }

    my %final_results;
    for my $cod ( sort keys %raw_results ) {
        my $t = 0;
        map { $t += ( $_->{'s'} / $_->{'n'} ) } @{ $raw_results{$cod} };
        $final_results{$cod} = { 's' => $t, 'n' => 3 - $t };
    }

    return \%final_results;
}

=method comp_codons

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

=cut

sub comp_codons {
    my $self = shift;
    my $cod1 = shift;
    my $cod2 = shift;
    my $pos  = shift;

    my $syn_changes = $self->syn_changes();
    my $codon2aa    = $self->codon2aa();

    my $syn_cnt = 0;    # total synonymous changes
    my $nsy_cnt = 0;    # total non-synonymous changes

    # ignore codon if beeing compared with gaps!
    if ( $cod1 =~ /\-/ or $cod2 =~ /\-/ ) {
        return ( $syn_cnt, $nsy_cnt );
    }

    # check codons
    for ( $cod1, $cod2 ) {
        if ( !exists $codon2aa->{$_} ) {
            croak Dump( { cod1 => $cod1, cod2 => $cod2 } ), "Wrong codon\n";
        }
    }

    # check codon position
    if ( defined $pos ) {
        if ( none { $_ == $pos } ( 0 .. 2 ) ) {
            croak Dump( { pos => $pos } ), "Wrong codon position\n";
        }
    }

    my %mutator = (
        2 =>    # codon positions to be altered
                # depend on which is the same
            {
            0 => [ [ 1, 2 ], [ 2, 1 ] ],
            1 => [ [ 0, 2 ], [ 2, 0 ] ],
            2 => [ [ 0, 1 ], [ 1, 0 ] ],
            },
        3 =>    # all need to be altered
            [
            [ 0, 1, 2 ],
            [ 1, 0, 2 ],
            [ 0, 2, 1 ],
            [ 1, 2, 0 ],
            [ 2, 0, 1 ],
            [ 2, 1, 0 ],
            ],
    );

    my ( $diff_cnt, $codon_pos ) = $self->count_diffs( $cod1, $cod2 );

    if ( $diff_cnt == 0 ) {    # ignore if codons are identical
    }
    elsif ( $diff_cnt == 1 ) {    # In $codon_pos where bases are different
        if ( !defined $pos or $codon_pos == $pos ) {
            $syn_cnt = $syn_changes->{$cod1}{$cod2};
            $nsy_cnt = 1 - $syn_changes->{$cod1}{$cod2};
        }
    }
    elsif ( $diff_cnt == 2 ) {    # In $codon_pos where bases are the same
        my ( $s_cnt, $n_cnt ) = ( 0, 0 );
        my $pathway = 2;          # will stay 2 unless there are stop codons
                                  #   at intervening point
    PATH: for my $perm ( @{ $mutator{2}{$codon_pos} } ) {
            my $altered = $cod1;
            my $prev    = $cod1;
            my ( $sub_s_cnt, $sub_n_cnt ) = ( 0, 0 );

            for my $mut_i (@$perm) {    # index of codon mutated
                my $mut_base = substr( $cod2, $mut_i, 1 );
                substr( $altered, $mut_i, 1, $mut_base );
                if ( $codon2aa->{$altered} eq '*' ) {
                    $pathway--;
                    next PATH;          # abadon this pathway
                }
                else {
                    if ( !defined $pos or $mut_i == $pos ) {
                        $sub_s_cnt += $syn_changes->{$prev}{$altered};
                        $sub_n_cnt += 1 - $syn_changes->{$prev}{$altered};
                    }
                }
                $prev = $altered;
            }

            $s_cnt += $sub_s_cnt;
            $n_cnt += $sub_n_cnt;
        }
        if ( $pathway != 0 ) {
            $syn_cnt = ( $s_cnt / $pathway );
            $nsy_cnt = ( $n_cnt / $pathway );
        }
    }
    elsif ( $diff_cnt == 3 ) {
        my ( $s_cnt, $n_cnt ) = ( 0, 0 );
        my $pathway = 6;    # will stay 6 unless there are stop codons
                            #   at intervening point
    PATH: for my $perm ( @{ $mutator{'3'} } ) {
            my $altered = $cod1;
            my $prev    = $cod1;
            my ( $sub_s_cnt, $sub_n_cnt ) = ( 0, 0 );

            for my $mut_i (@$perm) {    #index of codon mutated
                my $mut_base = substr( $cod2, $mut_i, 1 );
                substr( $altered, $mut_i, 1, $mut_base );
                if ( $codon2aa->{$altered} eq '*' ) {
                    $pathway--;
                    next PATH;          # abadon this pathway
                }
                else {
                    if ( !defined $pos or $mut_i == $pos ) {
                        $sub_s_cnt += $syn_changes->{$prev}{$altered};
                        $sub_n_cnt += 1 - $syn_changes->{$prev}{$altered};
                    }
                }
                $prev = $altered;
            }

            $s_cnt += $sub_s_cnt;
            $n_cnt += $sub_n_cnt;
        }

        # calculate number of synonymous/non synonymous mutations for that
        # codon and add to total
        if ( $pathway != 0 ) {
            $syn_cnt = ( $s_cnt / $pathway );
            $nsy_cnt = ( $n_cnt / $pathway );
        }
    }    # endif $diffcnt = 3

    return ( $syn_cnt, $nsy_cnt );
}

# counts the number of nucleotide differences between 2 codons
# when 1 nucleotide is different, returns this value plus the codon index
#   of which nucleotide is different
# when 2 nucleotides are different, returns this value plus the codon index
#   of which nucleotide is the same
# So comp_codons() knows which nucleotides to change or not to change
sub count_diffs {
    my $self = shift;
    my $cod1 = shift;
    my $cod2 = shift;

    my $cnt = 0;
    my $return_pos;
    my @sames;    # store same base position
    my @diffs;    # store diff base position

    if ( length $cod1 != 3 or length $cod2 != 3 ) {
        carp Dump( { cod1 => $cod1, cod2 => $cod2 } ), "Codon length error\n";
        return ( $cnt, $return_pos );
    }

    # just for 2 differences
    for ( 0 .. 2 ) {
        if ( substr( $cod1, $_, 1 ) ne substr( $cod2, $_, 1 ) ) {
            $cnt++;
            push @diffs, $_;
        }
        else {
            push @sames, $_;
        }
    }

    if ( $cnt == 1 ) {
        $return_pos = $diffs[0];
    }
    elsif ( $cnt == 2 ) {
        $return_pos = $sames[0];
    }

    return ( $cnt, $return_pos );
}

sub translate {
    my $self  = shift;
    my $seq   = shift;
    my $frame = shift;

    # check $frame
    if ( defined $frame ) {
        if ( none { $_ == $frame } ( 0 .. 2 ) ) {
            croak Dump( { frame => $frame } ), "Wrong frame\n";
        }
    }
    else {
        $frame = 0;
    }

    $seq = substr( $seq, $frame );    # delete first $frame bases from $seq
    my $offset = length($seq) - ( length($seq) % 3 );
    substr( $seq, $offset, length($seq), '' );    # now $seq is 3n bp

    my $codon_table = $self->codon_table();

    my $peptide = $codon_table->translate($seq);

    return $peptide;
}

=method is_start_codon

      Usage : my $bool = $obj->is_start_codon('ATG')
    Purpose : returns true (1) for codons that can be used as a
            :   translation start, false (0) for others.
    Returns : boolean
 Parameters : Codon
     Throws : no exceptions
   Comments : none
   See Also : n/a

=cut

sub is_start_codon {
    my $self = shift;
    my $cod  = shift;

    my $codon_table = $self->codon_table();

    return $codon_table->is_start_codon($cod);
}

=method is_ter_codon

      Usage : my $bool = $obj->is_ter_codon('GAA')
    Purpose : returns true (1) for codons that can be used as a
            :   translation terminator, false (0) for others.
    Returns : boolean
 Parameters : Codon
     Throws : no exceptions
   Comments : none
   See Also : n/a

=cut

sub is_ter_codon {
    my $self = shift;
    my $cod  = shift;

    my $codon_table = $self->codon_table();

    return $codon_table->is_ter_codon($cod);
}

=method is_unknown_codon

      Usage : my $bool = $obj->is_unknown_codon('GAJ')
    Purpose : returns true (1) for codons that are valid,
            :   true (1) for others.
    Returns : boolean
 Parameters : Codon
     Throws : no exceptions
   Comments : none
   See Also : n/a

=cut

sub is_unknown_codon {
    my $self = shift;
    my $cod  = shift;

    my $codon_table = $self->codon_table();

    return $codon_table->is_unknown_codon($cod);
}

1;    # Magic true value required at end of module

__END__

=head1 DESCRIPTION

AlignDB::Codon is a part of AlignDB package. It provides methods to translate
sequences and calculate Dn/Ds with different codon tables.

For more information, see L<AlignDB>.
