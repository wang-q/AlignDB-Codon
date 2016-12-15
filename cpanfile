requires 'Moose';
requires 'List::MoreUtils';
requires 'YAML::Syck';
requires 'AlignDB::IntSpan';
requires 'perl', '5.010001';

on test => sub {
    requires 'Test::More', 0.88;
    requires 'Bio::Align::DNAStatistics';
};
