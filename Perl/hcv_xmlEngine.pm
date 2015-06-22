#!/usr/bin/perl -w
# A module that dumps db data.
=head1 NAME
   hcv_xmlEngine
   Version 1.0.1
   Created 05.25.2011
   by Mark Evans
   Revised 08.04.2011
=head1 SYNOPSIS
   require hcv_xmlEngine;
=head1 DESCRIPTION
   Module to process xml file.
=head1 COPYRIGHT
   Copyright Monogram BioSciences 2011.
=head1 NOTES
=cut

###############################################################
use Data::Dumper;
use strict;
use DBI;
use XML::Simple;
use Mutations;
my $ENGINE_VERSION = 1.0;

my $GS_ORA_DB          = $ENV{GS_ORA_DB};
my $GS_ORA_DB_USER     = $ENV{GS_ORA_DB_USER};
my $GS_ORA_DB_PASSWD   = $ENV{GS_ORA_DB_PASSWD};

my $GS_MYSQL_DB        = $ENV{GS_MYSQL_DB};
my $GS_MYSQL_DB_USER   = $ENV{GS_MYSQL_DB_USER};
my $GS_MYSQL_DB_PASSWD = $ENV{GS_MYSQL_DB_PASSWD};

################################################################ Begin package Sample
package Sample;

###############
# Sample::new 
###############
sub new 
{
   my $class = shift;
   my $self = {};
   
   bless($self, $class);
   $self->{samples} = ();
   return $self;
}


####################
# Sample::AUTOLOAD 
####################
sub AUTOLOAD 
{
   no strict;
   if ($AUTOLOAD=~/(\w+)$/) 
   {
      my $field = $1;
      *{$field} = sub 
                  {
                     my $self = shift;
                     @_ ? $self->{$field} = shift : $self->{$field};
                  };
      &$field(@_);
   }
   else 
   {
      Carp::confess("Cannot figure out field name from '$AUTOLOAD'");
   }
}


###################
# Sample::destroy 
###################
sub DESTROY { }


################################################################ Begin package Drug
package Drug;

#############
# Drug::new 
#############
sub new 
{
   my $class = shift;
   my $self = {};
   
   bless($self, $class);
   $self->{xml} = $_[0];
   return $self;
}


##################
# Drug::AUTOLOAD 
##################
sub AUTOLOAD 
{
   no strict;
   if ($AUTOLOAD=~/(\w+)$/) 
   {
      my $field = $1;
      *{$field} = sub 
                  {
                     my $self = shift;
                     @_ ? $self->{$field} = shift : $self->{$field};
                  };
      &$field(@_);
   }
   else 
   {
      Carp::confess("Cannot figure out field name from '$AUTOLOAD'");
   }
}

#################
# Drug::DESTROY 
#################
sub DESTROY { }


################################################################ Begin package Function
package Function;

#################
# Function::new
#################
sub new 
{
   my $class = shift;
   my $self = {};
   
   bless($self, $class);
   $self->{xml} = $_[0];
   return $self;
}

######################
# Function::AUTOLOAD
######################
sub AUTOLOAD 
{
   no strict;
   if ($AUTOLOAD=~/(\w+)$/) 
   {
      my $field = $1;
      *{$field} = sub 
                  {
                     my $self = shift;
                     @_ ? $self->{$field} = shift  : $self->{$field};
                  };
      &$field(@_);
   }
   else 
   {
      Carp::confess("Cannot figure out field name from '$AUTOLOAD'");
   }
}

#####################
# Function::Destroy
#####################
sub DESTROY { }


################################################################ Begin package xmlEngine
package xmlEngine;

=head2 new()
=item I<Description:>
   This procedure creates a new xmlEngine object.
=item I<Requires:>
   Parsed XML::Simple object
=item I<Returns:>
   Reference to an xmlEngine object.
=item I<Side Effects:>
   None
=item I<Example:>
   $sample = new xmlEngine();
=cut

##################
# xmlEngine::new
##################
sub new 
{
   my $class = shift;
   my $self = {};
   
   bless($self, $class);
   $self->{xml} = $_[0];
   $self->version($ENGINE_VERSION);
   return $self;
}


##############################
# xmlEngine::getTotalVersion
##############################
sub getTotalVersion 
{
   my ($self,$version);
   $self = shift; ##ref to my object
   $version=$self->xml->{totalVersion};
   return($version);
}


###############################
# xmlEngine::doResistanceFlag
###############################
sub doResistanceFlag 
{
   my ($result)=$_[0];
   if ($result eq 'Sensitive')
   {
      return 0;
   } 
   elsif ($result eq 'Resistant')
   {
      return 1;
   }  
   elsif ($result eq 'Resistance Possible') 
   {
      return 2;
   }
}


#########################
# xmlEngine::codeResult 
#########################
sub codeResult 
{
   my($result,$rams) = @_;
   my ($c,@ramlist,$dr,@r,$r,@v,$v);
   if ($result eq 'Sensitive')
   {
      return 0;
   } 
   elsif ($result eq 'Resistant' || $result eq 'Resistance Possible')
   {
#$DB::single=2; # perl debugger
      $c = 0;
      $rams =~ s/\s+//g;                                     # get rid of spaces to prep for split
      @ramlist = split(/;/,$rams);                       # split into array of region[rams]
      foreach $dr (@ramlist)                             # Loop through regions
      {
         @r = split(/:/,$dr);
         if ($r[1] =~ /.*None.+/)
         {
            $c += 0;
         }
         else 
         {
            @v = split(/,/,$r[1]);
            $c += $#v + 1;
         }
      }
      return $c;
   }  
   else
   {
      return 1;
   }
}


############################
# xmlEngine::getScoreIndex 
############################
sub getScoreIndex 
{
   my ($functions) = @_;   
   my @function_array = @{$functions};
   my ($index, $j, $name);
   
   $index = $#function_array;
   
   for($j=0;$j<=$#function_array;$j++) 
   {
      $name = $function_array[$j]->name;
      
      if ($name =~ /Score/) 
      {
         $index = $j;
      }
   }
   return $index;
}


######################
# xmlEngine::doMySQL 
######################
sub doMySQL 
{
   my ($self, @samples,@drugs,$i,$j,$k,@functions,$valueString, $pr_muts, $rt_muts, $int_muts, $index);
   $self=shift;
   @samples = @{$self->samples};                                                                    
   
   my $mysqlDB = DBI->connect($GS_MYSQL_DB,$GS_MYSQL_DB_USER,$GS_MYSQL_DB_PASSWD) || (print "Can't connect to MySQL database: $!\n<RETURNSTATUS=0>\n" & die);
   
   for($i=0;$i<=$#samples;$i++)
   {
      @drugs=@{$samples[$i]->drugs};
      
      $valueString = "DELETE FROM GT_RULE_RESULT_ROWS WHERE SAMPLE_ID='" . $samples[$i]->id . "'";
      $mysqlDB->do($valueString);
      
      for($j=0;$j<=$#drugs;$j++)
      {
         
         $index = getScoreIndex($drugs[$j]->functions);
         
         $valueString =  'Insert into GT_RULE_RESULT_ROWS
                          (SAMPLE_ID, DRUG_ABRV_NAME, DRUG_CLASS, COMMENTS, RAMs, RESULT, RESISTANT_FLG, MUTATION_SCORE, 
                          CREATED_DATE, RULE_VERSION, ENGINE_VERSION)
                          Values
                          (';
         
         $valueString .= "'".$samples[$i]->id . "'," . "'" . $drugs[$j]->drug_abrv . "',";
         $valueString .= "'".$drugs[$j]->drug_class . "'," . "'" . $drugs[$j]->comments . "',";
         $valueString .= "'".$drugs[$j]->finalrams. "'," . "'" . codeResult($drugs[$j]->result,$drugs[$j]->finalrams) . "',";
         $valueString .= "'".doResistanceFlag($drugs[$j]->result) . "'," . $drugs[$j]->drug_mut_score . ", " . "NOW(),";
         $valueString .= "'".$self->xml->{totalVersion}  . "'," . "'" . $self->version . "' )";
         
         $mysqlDB->do($valueString);
      }
      # print STDERR $samples[$i]->id . "\t" . $drugs[$j]->rams . "n";
   }
   
   $mysqlDB->disconnect();
}


##########################
# xmlEngine::printReport
##########################
sub printReport 
{
   my ($self, @samples,@drugs,$i,$j,$k,@functions,$valueString, $index,$hash,%gene_hash,@genes,$final_rams,$gene,$ram_ref,%ram_hash,@ram_pos,$ram_list);
   $self=shift;
   
   @samples = @{$self->samples};        
   my $header = "sample_id\tdrug_abrv\tdrug_class\tdrug_comments\tRAM_list\tdrug_result\tresistance_flg\txml_version\tversion\tdrug_score\n";
   print STDERR $header;                                       
   for($i=0;$i<=$#samples;$i++)
   {
      @drugs=@{$samples[$i]->drugs};

      for($j=0;$j<=$#drugs;$j++) 
      {
         $index = getScoreIndex($drugs[$j]->functions);

#       $DB::single=2; # perl debugger
         $valueString='';
         $valueString = "'".$samples[$i]->id . "'\t" . "'" . $drugs[$j]->drug_abrv . "'\t";
         $valueString .= "'".$drugs[$j]->drug_class . "'\t" . "'" . $drugs[$j]->comments . "'\t";
         $valueString .= "'".$drugs[$j]->finalrams. "'\t" . "'" . $drugs[$j]->result . "'\t";
         $valueString .= "'".doResistanceFlag($drugs[$j]->result) . "'\t";
         $valueString .= "'".$self->xml->{totalVersion}  . "'\t" . "'" . $self->version . "')\t";
         $valueString .= $drugs[$j]->drug_mut_score."\n";
         
         print STDERR $valueString;
      }
   }
}

####################
# xmlEngine::doOracle
####################
sub doOracle
{

   my ($self, @samples,@drugs,$i,$j,$k,@functions,$valueString, $index);
   $self=shift;
   @samples = @{$self->samples};

   my $dsn = $self->db_dsn;      #'dbi:Oracle:dev';
   my $db_user = $self->db_user; #'rule_engine';
   my $db_pass = $self->db_pass; #'rule_engine';
   
   my $dbh;
   eval
   {
      $dbh = DBI->connect($dsn, $db_user, $db_pass) or die "DB Error: $DBI::errstr";
   };
   if($@)
   {
      print "Can't connect to Oracle database: $!\n<RETURNSTATUS=0>\n" & die;
   }
   
   $dbh->{AutoCommit}      = 1;
   $dbh->{RaiseError}      = 1;
   $dbh->{ora_check_sql}   = 0;
   $dbh->{RowCacheSize}    = 64;
   
   for($i=0;$i<=$#samples;$i++) 
   {
      @drugs=@{$samples[$i]->drugs};
      doDeleteResult($dbh, $samples[$i]->id);
      for($j=0;$j<=$#drugs;$j++) 
      {
         $index = getScoreIndex($drugs[$j]->functions);
#         my($hash,%gene_hash,@genes,$gene);
#         $hash = @{$drugs[$j]->functions}[$index]->rams;    # get compound RAM hash reference
#         %gene_hash = %{$hash};                             # level one dereference
#         @genes = sort keys %gene_hash;                     # Extract gene names that are used as keys
         
         my $rams = $drugs[$j]->finalrams;                      # get compound RAM string
         $rams =~ s/\s+//;                                     # get rid of spaces to prep for split
         my @ramlist = split(/;/,$rams);                       # split into array of region[rams]
         pop(@ramlist);                                        # remove empty last index (split artifact)
#          foreach $gene (@genes)
         foreach my $dr (@ramlist)                             # Loop through regions
         {
            my @r = split(/:/,$dr);                            # Seperate region from ram list
            doInsertResult(
                              dbh   => $dbh,
                              rule_engine_id => $self->version,
                              created_by     => 'rule_engine',
                              drug_abrv_name => $drugs[$j]->drug_abrv,
                              drug_class     => $drugs[$j]->drug_class,
                              comments       => $drugs[$j]->comments,                                 # <== and check this one too please
                              resistance_mutations => $r[1],
                              result         => codeResult($drugs[$j]->result,$drugs[$j]->finalrams), # <== Duc, check that this is correct
                              resistant_flg  => doResistanceFlag($drugs[$j]->result),
                              rule_version   => $self->xml->{totalVersion},
                              aliquotID      => $samples[$i]->id,
                              region         => $r[0]
                           )
         }
      }
      # print STDERR $samples[$i]->id . "\t" . $drugs[$j]->rams . "n";
   }
   $dbh->disconnect();
}



############################
# xmlEngine::doDeleteResult
############################
sub doDeleteResult
{
   my $dbh       = shift;
   my $aliquotID = shift;

   my $retVal;
   my $error_message;
   my $sql = qq{
      BEGIN
         :retVal := ENGINE_RESULTS.deleteHCVRuleResult(:p_aliquotID);
      END;
   };
   my $sth = $dbh->prepare($sql) or die "DB Error: $DBI::errstr";
   $sth->bind_param_inout(":retVal", \$retVal, 2);
   $sth->bind_param(":p_aliquotID", $aliquotID);
   eval
   {
      $sth->execute() or die "DB Error: $DBI::errstr";
   };
   if ($@)
   {
      die "DB Error: $@";
   }
   if ($retVal ne 'OK') {
      print "Can't delete previous result for aliquot $aliquotID: $!\n<RETURNSTATUS=0>\n" & die;
   }
   $sth->finish();
}


###########################
# xmlEngine::doInsertResult
###########################
sub doInsertResult
{
   my %params = @_;
   my $dbh              = $params{'dbh'};
   my $rule_engine_id   = $params{'rule_engine_id'};
   my $created_by       = $params{'created_by'};
   my $drug_abrv_name   = $params{'drug_abrv_name'};
   my $drug_class       = $params{'drug_class'};
   my $comments         = $params{'comments'};
   my $resistance_mutations   = $params{'resistance_mutations'};
   my $result           = $params{'result'};
   my $resistant_flg    = $params{'resistant_flg'};
   my $rule_version     = $params{'rule_version'};
   my $aliquotID        = $params{'aliquotID'};
   my $region           = $params{'region'};
   
   my $retVal;
   my $error_message;
   my $sql = qq{
      BEGIN
         :retVal := ENGINE_RESULTS.insertHCVRuleResult(
            :p_rule_engine_id,
            :p_created_by,
            :p_drug_abrv_name,
            :p_drug_class,
            :p_comments,
            :p_resistance_mutations,
            :p_result,
            :p_resistant_flg,
            :p_rule_version,
            :p_aliquotID,
            :p_region
            );
      END;
   };
   my $sth = $dbh->prepare($sql) or die "DB Error: $DBI::errstr";
   $sth->bind_param_inout(":retVal", \$retVal, 2);
   $sth->bind_param(":p_rule_engine_id", $rule_engine_id);
   $sth->bind_param(":p_created_by", $created_by);
   $sth->bind_param(":p_drug_abrv_name", $drug_abrv_name);
   $sth->bind_param(":p_drug_class", $drug_class);
   $sth->bind_param(":p_comments", $comments);
   $sth->bind_param(":p_resistance_mutations", $resistance_mutations);
   $sth->bind_param(":p_result", $result);
   $sth->bind_param(":p_resistant_flg", $resistant_flg);
   $sth->bind_param(":p_rule_version", $rule_version);
   $sth->bind_param(":p_aliquotID", $aliquotID);
   $sth->bind_param(":p_region", $region);
   eval
   {
      $sth->execute() or die "DB Error: $DBI::errstr";
   };
   if ($@)
   {
      die "DB Error: $@";
   }
   if ($retVal ne 'OK') {
      print "Cannot insert result for aliquot $aliquotID: $! retVal=$retVal\n<RETURNSTATUS=0>\n";
   }
   $sth->finish();
}


########################
# xmlEngine::doSamples 
########################
sub doSamples 
{
   my ($self,%mut_lists,$aliquotID,%samples,$sample,$key,$tag,$sample_obj,@drug_objs);
   $self=shift;
   %samples = %{$self->xml->{sample}};
   
   # Loop through each sample in xml file
   foreach $key (keys(%samples)) 
   {
      # skip dummy sample
      next if ($key eq "dummy");
      
      # Get all gene mutations lists dynamically and add to hash with gene keys
      # keys(%{$samples{$key}})  <= gives inner hash keys NS3Mutations, etc
      
      foreach $tag (keys(%{$samples{$key}}))
      {
         if ($tag =~ m/(.+)Mutations/)
         {
            my $region = $1;
            my $smpl = $samples{$key}{$tag};
            $smpl =~ s/\^/X/g;   ## remove ^ sign for deletion -> put back in after parsing string
            my $obj = new Mutations();
            my $obj2 = $obj->ProcessGeno($smpl);
 #           $DB::single=2; # perl debugger
            $obj2->mutString($smpl);
            $mut_lists{$region} = $obj2;
         }
      }
      
      $sample_obj = new Sample();
      $sample_obj->id($key);
      $sample_obj->subtype(uc $samples{$key}->{Subtype});
      
      # Eval drugs for each mutation set
      @drug_objs = $self->doDrugs(\%mut_lists,$sample_obj);
      
      if ($#drug_objs == -1) 
      { 
         return(0); 
      }      
      
      push(@{$sample_obj->{drugs}},@drug_objs);
      push(@{$self->{samples}},$sample_obj);
   }
   # success
   return(1);
}


##########################
# xmlEngine::doOperation  
##########################
sub doOperation 
{
   my ($var1,$var2,$op) =@_;
   
   if($op eq "GE") 
   {
      return($var1>=$var2);
   } 
   elsif ($op eq "LE") 
   {
      return($var1<=$var2);
   } 
   elsif ($op eq "EQ") 
   {
      return($var1==$var2);
   } 
   elsif ($op eq "GT") 
   {
      return($var1>$var2);
   }
   elsif ($op eq "LT") 
   {
      return($var1<$var2);
   }
   else 
   {
      print "Operation=$op not supported!\n";
      print "<RETURNSTATUS=0>\n";
      return(0);
   }
}


###########################
# xmlEngine::doConditions 
###########################
sub doConditions 
{
   my ($self,$drug,$function_objs) = @_;
   my (@conditions,@sortedConditions,$key,@functions,$priority,$tmp,$i, $operation,$term,$term1,$term2,$tval1,$tval2,$val1,$val2,$cType,$action,$j);
   
   @functions=@{$function_objs};
   @conditions = @{$drug->{condition}};
   @sortedConditions = sort {$a->{"cPriority"} <=> $b->{"cPriority"}} @conditions;
   
   for($i=0;$i<=$#sortedConditions;$i++)
   {
      $operation = $sortedConditions[$i]->{operation};
      $term1=$sortedConditions[$i]->{term}->{tOrder1}->{tValue};
      $term2=$sortedConditions[$i]->{term}->{tOrder2}->{tValue};
      $tval1=$tval2='';
      
      if($term1 =~ /^[\d.]+$/) 
      {
         $tval1 = $term1;
      }
      if($term2 =~ /^[\d.]+$/) 
      {
         $tval2 = $term2;
      }
      
      for($j=0;$j<=$#functions;$j++)
      {
         # print STDERR "Name=". $functions[$j]->name."\n";
         if($functions[$j]->name eq $term1) 
         {
            $tval1 = $functions[$j]->result;
         }
         if($functions[$j]->name eq $term2) 
         {
            $tval2 = $functions[$j]->result;
         }
      }
      if($tval1 eq '' || $tval2 eq '' || $operation eq '') 
      {
         print "Problem with condition tval1=$tval1, tval2=$tval2, operation=$operation\n";
         print "<RETURNSTATUS=0>\n";
         return(0); 
      } 
      elsif (doOperation($tval1,$tval2,$operation)) 
      {
         # Concatenate the winning fxn name with the outcome for later use
         return($term1.":".$sortedConditions[$i]->{action});
      }
   }
   return("Sensitive");
}


#########################
# xmlEngine::doComments
#########################
sub doComments 
{
   my ($self,$drug,$function_objs, $myResult)=@_;
   my (@comments, @functions, $h,$j, $gene, $result, @pos,@com_keys,$mut_string,$ck,@pos_keys,%comment_muts, %comment_hash, %mutations, %display, $mut, $comm, $text);
   
   @functions = @{$function_objs};
   @comments = @{$drug->{comment}};
   $text =""; 
   
   # Create hash of comments, results by comment ID
   foreach $h (@comments)
   {
      if ( %{$h}->{commentID} != 99 )
      {
         $comment_hash{%{$h}->{commentID}} = [%{$h}->{commentText},%{$h}->{result}];
      }
   }
   
   # Begin assessing comment mutations
   %comment_muts = %{$functions[0]{comment_muts}};
   foreach $gene (sort keys(%comment_muts))                                    # Go through each gene looking for comment mutations
   {
      if (scalar keys %{$comment_muts{$gene}})                                 # Check if any mutations recorded for that gene, otherwise skip
      {
         %mutations = %{ $comment_muts{ $gene }};                              # Dereference the hash and get the mutations for that gene
         foreach $mut (sort keys %mutations)                                   # Iterate through each mutation in the list
         {
            @pos = split(/\|/,$mut);
            $result = $comment_hash{ $mutations{ $mut }}[1];                   
            $comm = $comment_hash{ $mutations{ $mut }}[0];
            
            # If the result conditions match
            if ($myResult eq $result || $result eq 'any')
            {
               if ($display{$comm}{$gene})                                    # Use comments as keys, then genes, then position
               {
                  if (exists $display{$comm}{$gene}{$pos[0]})                 # if more than one mutation at same position
                  {                                                           # add it to same entry. Not alpha sorted.
                     $display{$comm}{$gene}{$pos[0]} = $display{$comm}{$gene}{$pos[0]}.", ".$pos[1];
                  }
                  else
                  {
                     $display{$comm}{$gene}{$pos[0]} = $pos[1];
                  }
               }
               else
               {
                  $display{$comm}{$gene}{$pos[0]} = $pos[1];
               }
            }
         }
      }
   }
   
   # Build final comment line
   @com_keys = keys %display;                                           # get comment keys
   foreach $ck (@com_keys)
   {
      my @ms;                                                           # initialize mut_string array 
      my @gene_keys = sort keys %{$display{$ck}};                       # get gene keys
      foreach my $gk (@gene_keys)
      {
         $mut_string='';                                                # init mut_string for this gene
         @pos_keys = sort {$a <=> $b} keys %{$display{$ck}{$gk}};       # get position keys numerically sorted
         $mut_string .= "$gk: ";
         foreach my $pk (@pos_keys)
         {
            $mut_string .= "$display{$ck}{$gk}{$pk}, ";
         }
         $mut_string .= ";";                                            # add ; at end of string build for this gene
         $mut_string =~ s/, ;$/; /;                                     # remove last , before the ;
         $mut_string =~ s/,(?!.*,)/ and/;                               # replace the last , with the word 'and'
         push(@ms,$mut_string);                                         # store mut_string in @ms so it doesn't get overwritten by next gene
      }
      $text .= $ck.join(" ",@ms);                                       # Attach comment to list of mutations
   }
   $text =~ s/; $/./;                                                   # remove last ; and replace with .
   return $text;
}


######################
# xmlEngine::doDrugs 
######################
sub doDrugs 
{
   my ($self,%drugs,$drug_class,$key,$sample_obj,$drug_obj,@function_objs,@drug_objs,$result, $comments,@abrv,%mut_list);
   $self=shift;
   %mut_list = %{$_[0]};
   $sample_obj=$_[1];
   %drugs = %{$self->xml->{drug}};
   
   # Create loop to filter drug rules to contain only subtype-specific rules
   # or generic rules that don't match a subtype rule
   # e.g. if subtype = 1A, then keep BOC-1A and TVR, but eliminate BOC and BOC-1B
   
   
   # Since the drug_abrv is one level in, need to loop through this first to make
   # a loopup hash of drug names indexed by drug_abrv
   my %drugabrv;
   while (my($key,$value) = each %drugs)
   {
      $drugabrv{$value->{'drug_abrv'}} = $key;
   }
   
   # Now loop through drug rules and remove those that are not needed
   while (my ($key, $value) = each %drugs)
   {
      my @ids = split(/-/,$drugs{$key}{'drug_abrv'});
      if ($#ids > 0 and $sample_obj->{subtype} ne $ids[1])
      {
         delete $drugs{$key};
      }
      elsif ($#ids > 0 and $sample_obj->{subtype} eq $ids[1])
      {
         if ( exists $drugabrv{$ids[0]} )
         {
            delete $drugs{$drugabrv{$ids[0]}};
         }
      }
   }
   
   # Loop through each drug rule and process mutations
   foreach $key (keys(%drugs)) 
   {
      $drug_obj = new Drug();
      $drug_obj->name($key);
      @function_objs = $self->doFunctions($drugs{$key}, \%mut_list);
      
      my ($k,$hash,%gene_hash,@genes,$final_rams,$gene,$ram_ref,%ram_hash,@ram_pos,$ram_list);
      
      # Need to extract ram list from the embedded object
      ###################################################
      $hash = $function_objs[0]->rams;                    # get compound RAM hash reference
      %gene_hash = %{$hash};                             # level one dereference
      @genes = sort keys %gene_hash;                     # Extract gene names that are used as keys
      $final_rams='';
      
      # Process mutations for each gene
      foreach $gene (@genes)
      {
         $ram_ref = $gene_hash{$gene};                   # Get gene RAM hash reference
         %ram_hash = %{$ram_ref};                        # level two dereference
         @ram_pos = sort {$a <=> $b} keys %ram_hash;     # Get ram position keys and sort numerically
         $ram_list ='';
#         $DB::single=2; # perl debugger
         # Assemble final_ram display string
         if ($#ram_pos > -1)                 
         {
            foreach $k (@ram_pos)
            {
               if ($ram_list eq '') 
               { 
                 $ram_list = printRAMS($ram_hash{$k},$mut_list{$gene}->{mutString},$gene);
               }
               else 
               { 
                  $ram_list = $ram_list.", ".printRAMS($ram_hash{$k},$mut_list{$gene}->{mutString},$gene); 
               }
            }
         }
         else 
         { 
            $ram_list = "None"; 
         } 
         $final_rams = $final_rams."$gene: $ram_list; ";
      }
      $final_rams =~ s/X/\^/g;                           # Replace the X with the original ^ character.
      
      $drug_obj->finalrams($final_rams);
      
      
     
      # Want to get which rule fxn won and what the mut_score for it was
      my $cmpd_result = $self->doConditions($drugs{$key},\@function_objs);  
      
      # $cmpd_result is in fxn:result format unless = 'Sensitive', so need to break apart
      if ($cmpd_result eq "Sensitive")
      {
         $result = $cmpd_result;
         $drug_obj->rule_fxn('None');
         $drug_obj->drug_mut_score('0');
      }
      else
      {
         my @condition_result = split(/:/,$cmpd_result);
         $result = $condition_result[1];
         $drug_obj->rule_fxn($condition_result[0]);  # Which rule function won, resulting in final call
         my %k;
         # Find the match by fxn name in the function_objs, then get the corresponding result to store in drug_obj
         foreach my $k (@function_objs)
         {
            if (%{$k}->{'name'} eq $condition_result[0])
            {
               $drug_obj->drug_mut_score(%{$k}->{'result'});
            }
         } 
      }
      
      
      $comments = $self->doComments($drugs{$key},\@function_objs, $result);
      
      if ($#function_objs == -1) 
      { 
         # print error
         print "<RETURNSTATUS=0>\n";
         return(0); 
      }

      $drug_obj->result($result);
      $drug_obj->comments($comments);
      $drug_obj->drug_class($drugs{$key}->{drug_class});
      
      # Change those drug abbreviations that are subtype specific to be generic for reporting
      my @ids = split(/-/,$drugs{$key}{'drug_abrv'});
      if ($#ids > 0 )
      {
         $drug_obj->drug_abrv($ids[0]);
      }
      else
      { 
         $drug_obj->drug_abrv($drugs{$key}->{drug_abrv}); 
      }
      
      push(@{$drug_obj->{functions}},@function_objs);
      push(@drug_objs,$drug_obj);
   }
   return(@drug_objs);
}


#######################
# xmlEngine::doSumMax 
#######################
sub doSumMax 
{
   my ($mut_score,$vWeight,$mut_val,$vPosition)=@_;
   my ($pos);
   
   $pos = "Pos".$vPosition;
   if(defined($mut_score->$pos))
   {
      if($mut_score->$pos < $mut_val*$vWeight)
      {
         $mut_score->$pos($mut_val*$vWeight);
      }
   }
   else 
   {
      $mut_score->$pos($mut_val*$vWeight);
   }
   return($mut_score);
}


####################
# xmlEngine::doSum 
####################
sub doSum 
{
   my ($mut_score,$vWeight,$mut_val,$vPosition)=@_;
   my ($pos,$newScore,$oldScore);
   
   $pos = "Pos".$vPosition;
   
   if (defined($mut_score->$pos)) 
   {
      $oldScore = $mut_score->$pos;
   }
   $newScore = $mut_val*$vWeight;
   $mut_score->$pos($oldScore + $newScore);
   
   return($mut_score);
}


#############################
# xmlEngine::applyOperation 
#############################
sub applyOperation 
{
   my ($mut_score,$operation)=@_;
   my ($mykey,$total);
   
   $total = 0; 
   
   foreach $mykey (keys(%{$mut_score})) 
   {
      if ($operation eq "+" || $operation eq "ispositive") 
      {
         $total = $total + $mut_score->$mykey;
      } 
      elsif ($operation eq "-") 
      {
         $total = $total - $mut_score->$mykey;
      }
      elsif ($operation eq "*") 
      {
         $total = $total * $mut_score->$mykey;
      }
      elsif ($operation eq "/")
      {
         $total = $total / $mut_score->$mykey;
      }
      elsif ($operation eq "^") 
      {
         $total = $total ^ $mut_score->$mykey;
      }
      elsif ($operation eq "log") 
      {
         $total = log($mut_score->$mykey);
      }
      elsif ($operation eq "append") 
      {
         $total = $total . ",$mykey";
      }
   }
   if ($operation eq "ispositive") 
   {
      if ($total > 0) 
      {
         return(1);
      }
      else 
      { 
         return(0);
      }
   }
   else 
   {
      return($total);
   }
}


#####################
# xmlEngine::doTerm 
#####################
sub doTerm 
{
   my ($mut_score,$vWeight,$vName, $function_objs)=@_;
   my (@function_objs, $i, $value, $key);
   
   $value = 0;
   
   # lookup function name in function hash and get the result
   @function_objs = @{$function_objs};
   
   for($i=0;$i<=$#function_objs;$i++) 
   {
      if ($function_objs[$i]->name eq $vName) 
      {
         $value = $function_objs[$i]->result;         
      }
   }
   
   $mut_score->$vName($value * $vWeight);
   return($mut_score);
} 


#########################
# xmlEngine::doConstant
#########################
sub doConstant 
{
   my ($mut_score,$vWeight,$vName)=@_;
   my ($pos);
   $pos = $vName;
   $mut_score->$pos($vWeight);
   
   return($mut_score);
}


#####################
# xmlEngine::doRAMS
#####################
sub doRAMS 
{
   my ($rams,$vName,$vRAM)=@_;
   if ($vRAM) 
   {
      $rams->$vName($vName);
   }
   return($rams);
}


########################
# xmlEngine::printRAMS 
########################
sub printRAMS 
{
   my ($rams, $mutString, $gene)=@_;
   my ($tmpMutStr);

   $rams =~ s/\D//g;
   $tmpMutStr = $mutString; 
   $tmpMutStr =~ s/.*(\D$rams[A-Z\/]+).*/$1/g;
   $tmpMutStr =~ s/X/\^/g;
   return $tmpMutStr;
}


#######################
# xmlEngine::checkXML 
#######################
sub checkXML 
{
   my ($self);
}


##########################
# xmlEngine::doFunctions  
##########################
sub doFunctions 
{
   my ($self,$functions, $key, %functionHash, $priority, %results, @vars,$var,$i,$operation,%muts,$drug, $resultType,
   $mut_score,$weight,$vName,$mut_val,$vWeight,$vPosition,$finalScore,$vType,$rams,$vRAM,$function_obj,
   @function_objs,$nbFunctions,%sortedFunctionHash, $vCommentID,$gene,%seen_muts,$g,%comment_muts);
   $self = shift; ##ref to my object
   $drug=$_[0];
   %muts = %{$_[1]}; # this is now a hash ( gene => {mutations})
   $functions = $drug->{function};
   %functionHash=%{$functions};
   $nbFunctions=2;
         
   foreach $key (keys (%functionHash)) 
   {
      #counting $key;
      if ($key eq "name") 
      {
         $nbFunctions--;
      }
   }
   
   ## If more than 1 function in XML, we should never get $key=name
   #  sort {$functionHash{$a}->{"priority"} <=> $functionHash{$b}->{"priority"}} (keys(%functionHash))) {
   
   if($nbFunctions==2) 
   {
      %seen_muts=();
      %comment_muts=();

      foreach my $gene (keys(%muts))
      {
         if (not exists($seen_muts{$gene}))
         {
            $seen_muts{$gene} = new Mutations();
            $comment_muts{$gene} = new Mutations();
         }
      }
      
      foreach $key (sort {$functionHash{$a}->{"priority"} <=> $functionHash{$b}->{"priority"}} (keys(%functionHash))) 
      {
         $finalScore = 0;
         $function_obj = new Function();
         $resultType = $functionHash{$key}->{resultType};
         $priority = $functionHash{$key}->{priority};
         @vars = @{$functions->{$key}->{var}};
         $operation=$functions->{$key}->{operation};
         $mut_score = new Mutations();
         $weight = $functions->{$key}->{weighting};
         
         next if ($function_obj->name($key) eq "dummy");
         
         for($i=0;$i<=$#vars;$i++) 
         {
            #initialize in every loop
            $mut_val = $vRAM = $vCommentID = 0;
            $vName=$vars[$i]->{vName};
            $vWeight=$vars[$i]->{vValue};
            $vPosition=$vars[$i]->{vPosition};
            $vType = $vars[$i]->{varType};
            $gene = $vars[$i]->{vAttribute};
                     
            if($vType eq "constant") 
            {
               $mut_score = doConstant($mut_score,$vWeight,$vName);
            }
            elsif ($vType eq "mutation") 
            {
               if ($weight eq "max") 
               {
                  if($muts{$gene}->{$vName} > 0) 
                  {
                     $mut_val = 1;
                     $mut_score = doSumMax($mut_score,$vWeight,$mut_val,$vPosition);
                  }
               }
               elsif ($weight eq "fraction") 
               {
                  $mut_val = $muts{$gene}->{$vName};
                  $mut_score=doSum($mut_score,$vWeight,$mut_val,$vPosition);
               }
               
               # If vRAMs, record it in RAMs and if a comment mutation, record it
               if($muts{$gene}->{$vName} > 0) 
               {
                  $vRAM = $vars[$i]->{vRAM};
                  if ($vRAM == 1)
                  {
                     $seen_muts{$gene}{$vars[$i]->{vPosition}}=$vName;
                  }

                  $vCommentID = $vars[$i]->{vCommentID};
                  if ($vCommentID)
                  {
 #                    $DB::single=2; # perl debugger
                     my $combo_key = $vPosition."|".$vName;
                     $comment_muts{$gene}{$combo_key} = $vCommentID;                     
                  }
               }                     
            }
            elsif ($vType eq "term") 
            {
               $mut_score = doTerm($mut_score,$vWeight,$vName, \@function_objs);
            } ## end of term
         }
         
         $finalScore = applyOperation($mut_score,$operation);
         $function_obj->rams(\%seen_muts);
         $function_obj->result($finalScore);

 #        $DB::single=2; # perl debugger
         $function_obj->comment_muts(\%comment_muts);
         push(@function_objs,$function_obj);
      }
   }
   else
   {
      ## if there is only 1 drug, the XML structure is different (1 level down)
      
      print "ERROR: needs 2 functions! Number passed:". $nbFunctions . "\n";
      print "<RETURNSTATUS=0>\n";
      return(0);
   }
   return(@function_objs);
}


#######################
# xmlEngine::AUTOLOAD
#######################
sub AUTOLOAD 
{
   no strict;
   if ($AUTOLOAD=~/(\w+)$/) 
   {
      my $field = $1;
      *{$field} = sub 
                  {
                     my $self = shift;
                     # Test_expr ? If_True : Else_If_False
                     @_ ? $self->{$field} = shift : $self->{$field};
                  };
      # Call function $field and pass arguments $_      
      &$field(@_);
   }
   else 
   {
      Carp::confess("Cannot figure out field name from '$AUTOLOAD'");
   }
}


#####################
# xmlEngine:DESTROY
#####################
sub DESTROY { }
1;

