#!/usr/local/bin/perl

# Use this to remove carriage returns midrow in Oracle patient export data

$infile = $ARGV[0];
$outfile = $ARGV[1];

open(SESAME, "$infile")||die"Can't open  file $infile)";
open(SESAME2, "+>$outfile") || die "Can't create new file $outfile";

$count=0;
$keep_going='no';
$frag='';
while($line=<SESAME>)
{
     $c = $line =~ tr/"//;
     if($c == 1)
     {
        print "++++ found $line\n";
        if($keep_going eq 'no')
        {
           $DB::single=2; # perl debugger
           $line =~ s/\n//;
           $frag = $line;
           $keep_going = 'yes';
           print "this is bad...fixing...\n";
        }
        else
        {
           $DB::single=2; # perl debugger
           print $frag.$line."\n\n";
           print SESAME2 $frag.$line;
           $frag='';
           $keep_going='no';
        }
        
     }
     else 
     {
        print SESAME2 $line;
     }
}
     
close(SESAME);
close(SESAME2);
