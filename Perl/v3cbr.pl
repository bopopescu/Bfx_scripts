#!/usr/bin/perl -w
#
# V3CBR
# HIV V3 Tropism Case Based Reasoning
#
# usage: v3cbr.pl -q query_sequence -c case_sequences
#
# Mojgan Haddad
# Dec 2010
#

use strict;

use Getopt::Std;

# Include V3 specific modules
use Scoring;

# Parse arguments 
my %options = ();
getopts("s:q:c:n:d",\%options);

my $query = "";
my $cases = "";
my $display = 0;
my $numHits = 1;
my $scoring = "";
my $totalX4      = 0;
my $totalR5      = 0;
my $numX4Correct = 0;
my $numR5Correct = 0;


if($options{q}) {
	$query = $options{q};
} else { usage(); }

if($options{c}) {
	$cases = $options{c};
} else { usage(); }

if($options{s}) {
	$scoring = $options{s};
} else { usage(); }

if($options{n}) {
	$numHits = $options{n};
} 

if($options{d}) {
	$display = 1;
}

# Assign input arguments to variables
my $querySequencesFile = $query;
my $caseSequencesFile  = $cases;

# Parse input FASTA files and return lists of Sequence objects
my $querySequenceList = V3::Scoring::parseFile($querySequencesFile);
my $caseSequenceList  = V3::Scoring::parseFile($caseSequencesFile);
my $scoringMatrixRef  = V3::Scoring::parseScoringMatrix($scoring);

# print the cutoff and top n for the record
printf("ScoreMatrix: %s\t-> TopN: %3d\n", $scoring, $numHits);

# print header
print "v3id, tropism, CBRCall, X4SimilarityScore, SimilarCase\n";

## loop thru all cases, and run CBR
foreach my $querySequence (@{$querySequenceList}) {
	# Do Case Based Reasoning comparision against query and case library
	V3::Scoring::caseBasedScoring($querySequence, $caseSequenceList, $scoringMatrixRef);
	
	# Sort the cases by their score
	my $sortedCaseList = V3::Scoring::sortByCBRScore($caseSequenceList);
	
	# Show top n results; if -d argument is given
	if ($display) {
	   V3::Scoring::printToDisplay($querySequence, $sortedCaseList)
	}
	
	# Call TROPISM based on any top n is X4
	## for X4: if any of top n is X4 => Call X4
	my $tropismCall = "UNK";
	my $maxSimScore = 0;
	my $maxSimCase = "-";
	
	my $topHit  = shift(@{$sortedCaseList});	
	
	for(2 ... $numHits) {
	   if($topHit->getTropism() eq "DM" || $topHit->getTropism() eq "X4") {

	     if($tropismCall eq "UNK" ) {
		$maxSimScore = $topHit->cbrScore();
		$maxSimCase = $topHit->getAccession();
	     }	     

	     $tropismCall = "X4";
	   }
	   ## go to next case
	   $topHit  = shift(@{$sortedCaseList}) or last;
	} 
		
	## if none of top-n cases X4 => then call is R5
	if($tropismCall eq "UNK" ) {
	   $tropismCall = "R5";
	   
	   # store a negative score for R5
	   $topHit  = shift(@{$sortedCaseList});
	   $maxSimScore = $topHit->cbrScore() * -1;
	   $maxSimCase = $topHit->getAccession();
	}
	
	## print the final call
	if (!$display) {
		print $querySequence->getAccession(),", ",$querySequence->getTropism(), ", ";
		print $tropismCall,", ",$maxSimScore,", ", $maxSimCase,"\n";
	}
	
	## Calculate concordance, sensitivity, specificity
	my $query = $querySequence->getTropism();
	my $tropism="";
	
	if($query eq "X4" || $query eq "DM") {	
		$tropism="X4";
		$totalX4++;

		if($tropismCall eq $tropism ) {
			$numX4Correct++;	
		} 
	}

	$tropism = "R5";
	if($query eq $tropism) {
		$totalR5++;

		if($tropismCall eq $tropism ) {
			$numR5Correct++;
		} 
	}
}

printf("TotalR5:   %3d\nCorrect: %3d %4.1f\n", $totalR5, $numR5Correct, ($numR5Correct/$totalR5)*100);
printf("TotalX4:   %3d\nCorrect: %3d %4.1f\n", $totalX4, $numX4Correct, ($numX4Correct/$totalX4)*100);
printf("Overall Concordance: %4.1f\n", (($numR5Correct+$numX4Correct)/($totalR5+$totalX4))*100);


# Program usage message
sub usage {
	print "usage: v3cbr.pl -q query_sequences -c case_sequences -s scoring_matrix {-d} {-n #}\n";
	print "\t-d 	turn on verbose output\n";
	print "\t-n # 	set the number of top matches to check for X4\n";
	exit;
}

