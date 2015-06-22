#!/usr/bin/perl -w
#
# V3CBR
# HIV V3 Tropism Case Based Reasoning
#
# usage: v3CBREngine.pl -q query_sequence -c case_sequences
#
# Mojgan Haddad
# Dec 2010
# Version 1.0
################################################################

use strict;
use Getopt::Std;

# Include V3 specific modules
use Scoring;

# Parse arguments 
my %options = ();
getopts("s:q:c:n:d",\%options);

my $query = "";
my $cases = "CBR_CaseLibrary.txt";
my $display = 0;
my $numHits = 3; # any of top 3 X4 => X4
my $scoring = "CBR_Score_Matrix.txt";
my $totalX4 = 0;
my $totalR5 = 0;
my $numX4Correct = 0;
my $numR5Correct = 0;

if($options{q}) 
{
   $query = $options{q};
} 
else { usage(); }


# Assign input arguments to variables
my $querySequencesFile = $query;
my $caseSequencesFile  = $cases;

# Parse input FASTA files and return lists of Sequence objects
my $querySequenceList = V3::Scoring::parseFile($querySequencesFile);
my $caseSequenceList  = V3::Scoring::parseFile($caseSequencesFile);
my $scoringMatrixRef  = V3::Scoring::parseScoringMatrix($scoring);

# print header
print "v3id\tCBRCall\tX4SimilarityScore\tSimilarCaseID\n";

## loop thru all cases, and run CBR
foreach my $querySequence (@{$querySequenceList}) 
{
   # Do Case Based Reasoning comparision against query and case library
   V3::Scoring::caseBasedScoring($querySequence, $caseSequenceList, $scoringMatrixRef);
   
   # Sort the cases by their score
   my $sortedCaseList = V3::Scoring::sortByCBRScore($caseSequenceList);
   
   # Show top n results; if -d argument is given
   if ($display) 
   {
      V3::Scoring::printToDisplay($querySequence, $sortedCaseList)
   }
   
   # Call TROPISM based on whether any of top n seq is X4
   # Criteria for X4: if any of top n is X4 => Call as X4
   my $tropismCall = "UNK";
   my $maxSimScore = 0;
   my $maxSimCase = "-";
   
   my $topHit  = shift(@{$sortedCaseList});
   my $nextHit = $topHit;

   for(2 ... $numHits) {
      if($nextHit->getTropism() eq "DM" || $nextHit->getTropism() eq "X4") {

         if($tropismCall eq "UNK" ) {
	   $maxSimScore = $nextHit->cbrScore();
	   $maxSimCase = $nextHit->getAccession();
         }	     

         $tropismCall = "X4";
      }
      ## go to next case
      $nextHit  = shift(@{$sortedCaseList}) or last;
   } 

   ## if none of top-n cases X4 -> then call is R5
   if($tropismCall eq "UNK" ) {
     $tropismCall = "R5";

     # store a negative score for R5
     $maxSimScore = $topHit->cbrScore() * -1;
     $maxSimCase = $topHit->getAccession();
   }
   
   # print final call, score, and most similar case ID
   print $querySequence->getAccession(),"\t", $tropismCall,"\t",$maxSimScore,"\t", $maxSimCase,"\n";
}



# Program usage message
sub usage 
{
	print "usage: v3CBREngine.pl -q query_sequences \n";
	print "\t query with following data columns:\n";
	print "\t mgrm_acc, tropism, seq_len35, nt_r, nt_y, nt_w, nt_k, aa_x, charge_flg, iep_flg, small_flg, basic_flg, tnr5_hmm\n";
	exit;
}

