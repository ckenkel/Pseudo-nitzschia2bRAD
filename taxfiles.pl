#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without guarantees or restrictions

# -- check for dependencies
use Bio::DB::Taxonomy;

BEGIN { *Bio::DB::Taxonomy::flatfile::DESTROY = sub {} }

my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -nodesfile => 'nodes.dmp', -namesfile => 'names.dmp', -directory => 'DBdir'); 

system("mkdir ./bin");
my $new_dir = "./bin";
my $safe_dir = "./DBdir";
system("cp $safe_dir/* $new_dir/");
