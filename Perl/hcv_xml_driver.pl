#!/usr/bin/perl

# use module
use XML::Simple;
use hcv_xmlEngine;


# create object
$xml = new XML::Simple;

# read XML file
$samples = $xml->XMLin($ARGV[0]);

$engine = new xmlEngine($samples);
$engine->{db_dsn}  ='dbi:Oracle:dev';  # Replace with real dsn info
$engine->{db_user} ='rule_engine'; # Replace with real user
$engine->{db_pass} ='rule_engine'; # Replace with real pwd

$engine->doSamples();  # Calls main function in xmlEngine and kicks things off

$engine->printReport();
#$engine->doSQL();
#$engine->doOracle();

