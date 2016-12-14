requires 'Moose';
requires 'Bio::Tools::CodonTable';
requires 'List::MoreUtils';
requires 'YAML';
requires 'AlignDB::IntSpan';
requires 'perl', '5.010001';

on test => sub {
    requires 'Test::More', 0.88;
};
