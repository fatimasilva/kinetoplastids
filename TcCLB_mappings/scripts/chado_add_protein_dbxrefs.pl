#! /usr/bin/perl -w 

use strict;
use warnings;
use Getopt::Long;
use DBI;
use IO::Prompt;

# Loads dbxrefs like UniProtKB or RefSeq ids into chado using a tab file with the dbxref accessions 
# and GeneDB gene ids (one pair GeneID:dbxref per row)
#
# ID mapping should use current geneIDs, but the script can check for previous IDs !00 Confirm

my ($f_in, $dbxref_db, $help, $name_out, $dbname, $dbhost, $dbport, $dbuser);

GetOptions(
    'a|id_file=s'          => \$f_in,
    'b|dbxref_db=s'        => \$dbxref_db,
    'h|help'               => \$help,
    'o|name_out=s'         => \$name_out,
    'd|database_name=s'    => \$dbname,
    's|database_host=s'    => \$dbhost,
    'u|database_user=s'    => \$dbuser,
    'p|database_port=s'    => \$dbport,

);

(($f_in && $dbxref_db && $name_out && $dbname && $dbhost && $dbport && $dbuser) && !$help) || die <<USAGE;

Usage: $0
    -a|id_files  <tab file with gene ID (first column) and dbxref to load (second column)>
    -b|dbxref_db     <the name of the database the dbxrefs belong to: 'UniProtKB', 'RefSeq'>
    -o|name_out <prefix for output files>
    -h|help    <This help message>
    -d|database_name
    -u|database_user
    -s|database_host
    -p|database_port
USAGE

# log and out files
my $f_drop = $name_out . ".log";
my $f_out = $name_out . "_changes_chado.tsv";
my $f_test = $name_out . "_test.tsv";
my $f_his = $name_out . "_pep_feature_ids_history.tsv";

open DROP, ">$f_drop";

#### Database connection details
## Prompt the user for the db password
my $dbpass = prompt('Password:', -echo => '*');

$dbpass = "$dbpass";

## Connect to the db
my $dbi_connect = "DBI:Pg:dbname=$dbname;host=$dbhost;port=$dbport";

my $dbh = DBI->connect($dbi_connect, $dbuser, $dbpass,{RaiseError => 0, AutoCommit => 0}) or die "Can't connect to the database $dbname!\n";
print "Success connecting to database $dbname!\n";
my $errflag;

#### Prepare SQL statements
## Check if a dbxref exist in a db (by accession) and get the dbxref_id
my $s_sql_get_dbxref_id=$dbh->prepare('SELECT dbxref_id, db_id FROM dbxref WHERE accession =? and db_id =    ?');

## Insert new dbxref
my $s_sql_insert_dbxref=$dbh->prepare('INSERT INTO dbxref (db_id, accession) VALUES (?,?)');

## Get feature_id by uniquename
my $s_sql_get_feature_id=$dbh->prepare('SELECT feature_id FROM feature WHERE uniquename=?');

## Get feature_id by Previous ID (synonym.name) 
my $s_sql_with_previous_id=$dbh->prepare('SELECT feature_id FROM feature join feature_synonym using (feature_id) JOIN synonym using (synonym_id) where synonym.name=? AND synonym.type_id=26803');

## Get feature_id of mRNA (by feature_id of parent gene)
my $s_sql_mRNA_id=$dbh->prepare('SELECT subject_id from feature_relationship where object_id=?');

## Get feature_id of polypeptide (by feature_id of parent mRNA)
my $s_sql_pep_id=$dbh->prepare('SELECT subject_id from feature_relationship where object_id=? and type_id=69');

## Check if a dbxref already exist in a feature (by feature_id and dbxref_id)
my $s_sql_get_feature_dbxref_id=$dbh->prepare('SELECT feature_dbxref_id FROM feature_dbxref WHERE feature_id=? AND dbxref_id=?');

## Insert a feature_dbxref
my $s_sql_insert_feature_dbxref_id=$dbh->prepare('INSERT INTO feature_dbxref (feature_id, dbxref_id) VALUES (?,?)');


#### Get dbxrefs into hashes key:gene_id, value=dbxref_id
my $dbxrefs_ref = get_dbxrefs($f_in);
my %dbxrefs = %{$dbxrefs_ref};

## Check extraction
open TEST, ">>$f_test";
foreach my $key (sort (keys %dbxrefs)){

    print TEST "$key\t$dbxrefs{$key}\n";

}

#### Add dbxrefs to the polypeptide feaature in Chado, add the dbxrefs to the dbxref table too if needed
## Some db values, add more if needed (add keys also to the USAGE)

my %dbxref_db = ( 'UniProtKB' => 67,'RefSeq' => 406);
my $db_id = $dbxref_db{$dbxref_db};

open OUT, ">> $f_out";
open HIS, ">> $f_his";

Geneid: foreach my $gene_id (sort (keys %dbxrefs)){

    if (defined($errflag)){
        warn "Error inserting feature_dbxref: " . $DBI::errstr; 
        $dbh->rollback();
        # stop processing IDs
        last Geneid;
    }

    #### Get polypeptide feature_id
    my $feature_id = get_feature_id_relation($gene_id);

    #### Unfold multiple uniprot ids
    my @dbxrefs = split (/\|/, $dbxrefs{$gene_id});

    foreach my $dbxref (@dbxrefs){

        ## Get dbxref_id (insert new one if necessary)
        my $dbxref_id = get_dbxref_id($dbxref, $db_id);

        if(defined $feature_id && defined $dbxref_id){
        
            print "$gene_id\t$dbxrefs{$gene_id}\tCHECKING_FEATURE_DBXREF\n";
        
            ## Insert dbxref into feature_dbxref
            insert_feature_dbxref($gene_id, $dbxref, $feature_id, $dbxref_id);

        } elsif (defined $feature_id){ # Failed getting feature_id or $dbxref_id

            print DROP "$gene_id\t$dbxref\tDBXREF_NOT_FOUND\n";

        } else{

            print DROP "$gene_id\t$dbxref\tGENEID_NOT_FOUND\n";
        }

    }

} 

close TEST;
close HIS;
close OUT;
close DROP;

#### Commit changes and close db connection
## uncomment when ready
$dbh->commit() unless (defined($errflag)); 

## close connection
$dbh->disconnect();

########
####
sub get_dbxrefs{ # extract ids from tab file
    my ($s_file) = @_;
    my %s_dbxrefs = ();
    my @s_gene_ids = ();

    open INFILE, "< $s_file";

    while (<INFILE>){
        chomp;
        my $s_line = $_;

        unless ($s_line =~ /^#/){
            my ($s_gene_id, $s_dbxref) = split (/\t/, $s_line);
        
            ## Drop gene ids with no dbxrefs associated
            if ($s_dbxref ne ''){

                if (defined $s_dbxrefs{$s_gene_id}){

                        print DROP "$s_gene_id\t$s_dbxrefs{$s_gene_id}\t$s_dbxref\tMULTIPLE_DBXREF\n";

                        $s_dbxrefs{$s_gene_id} = $s_dbxrefs{$s_gene_id} . '|' . $s_dbxref; 

                }else{

                        $s_dbxrefs{$s_gene_id} = $s_dbxref;
                }

            }else{ # log gene_ids with no dbxref associated

                print DROP "$s_gene_id\tNO_DBXREF\n";

           }
        }
    }

    close INFILE;

    ## Return hash with key:gene_id, value=uniprot_id
    return \%s_dbxrefs;
}


####
sub get_dbxref_id { # check if the ID is present in chado, insert it if it is not, and then return the dbxref_id
    my ($s_dbxref, $s_db_id) = @_;
    my ($s_dbxref_id);

    #### Check if dbxref exists already in the db in chado (by db_id) 
#    my $s_sql_get_dbxref_id=$dbh->prepare('SELECT dbxref_id, db_id FROM dbxref WHERE accession =? and db_id =?');

    my @s_row =$dbh->selectrow_array($s_sql_get_dbxref_id,undef,$s_dbxref,$s_db_id);

    if (@s_row){ # if the dbxref exist in the dbs, get dbxref_id

        ($s_dbxref_id, $s_db_id) = @s_row;

    }else{ # if dbxref does not exist in the db, add dbxref

#        my $s_sql_insert_dbxref=$dbh->prepare('INSERT INTO dbxref (db_id, accession) VALUES (?,?)');

        unless($s_sql_insert_dbxref->execute($s_db_id, $s_dbxref)){

            ## error inserting
            print DROP "$s_dbxref\tERROR_INSERT_DBXREF\n"; # log

            $errflag=1;
            next Geneid; # Go to the beginning of the loop (outside subroutine)
        }

        print OUT "$s_dbxref\tADDED_DBXREF\n";

        #### Get dbxref_id
        $s_dbxref_id = $dbh->selectrow_array($s_sql_get_dbxref_id, undef, $s_dbxref, $s_db_id);

    }

    #### return dbxref_id 
    return $s_dbxref_id;
}

#### 
sub get_feature_id_relation {
    my ($s_gene_id) = @_;
    my ($s_gene_feature_id, $s_mRNA_feature_id, $s_pep_feature_id); 

    #### Get gene feature id, added checking for previous IDs because there were IDs being missed
#   my $s_sql_get_feature_id=$dbh->prepare('SELECT feature_id FROM feature WHERE uniquename=?');

    $s_gene_feature_id = $dbh->selectrow_array($s_sql_get_feature_id, undef, $s_gene_id);

    ## Add check for previous systematic IDs cvterm_id for previous_systematic_id=26803
    unless (defined $s_gene_feature_id){
#       my $s_sql_with_previous_id=$dbh->prepare('SELECT feature_id FROM feature join feature_synonym using (feature_id) JOIN synonym using (synonym_id) where synonym.name=? AND synonym.type_id=26803');

        $s_gene_feature_id=$dbh->selectrow_array($s_sql_with_previous_id, undef, $s_gene_id);

    }

    if (defined $s_gene_feature_id){

        #### Get polypeptide feature id 
#       my $s_sql_mRNA_id=$dbh->prepare('SELECT subject_id from feature_relationship where object_id=?');

        $s_mRNA_feature_id=$dbh->selectrow_array($s_sql_mRNA_id, undef, $s_gene_feature_id);

#       my $s_sql_pep_id=$dbh->prepare('SELECT subject_id from feature_relationship where object_id=? and type_id=69');

        $s_pep_feature_id=$dbh->selectrow_array($s_sql_pep_id, undef, $s_mRNA_feature_id);

        ## Catch the ones coming from the mRNAs directly 
         # (some previous systematic ID are on the mRNA features)

        unless (defined $s_pep_feature_id){

            $s_pep_feature_id=$dbh->selectrow_array($s_sql_pep_id, undef, $s_gene_feature_id);

        }

        print "$s_gene_id\t$s_pep_feature_id\tPEP_FEATURE_ID\n";

    } else{

        print TEST "$s_gene_id\tID_NOT_FOUND\n";

    }

    #### Return feature_id (polypeptide)
    return $s_pep_feature_id;

}

####
sub insert_feature_dbxref {
    my ($s_gene_id, $s_dbxref, $s_feature_id, $s_dbxref_id) = @_;

    #### Check if the feature_dbxref already exists
#   my $s_sql_get_feature_dbxref_id=$dbh->prepare('SELECT feature_dbxref_id FROM feature_dbxref WHERE feature_id=? AND dbxref_id=?');

    my $s_feature_dbxref_id=$dbh->selectrow_array($s_sql_get_feature_dbxref_id, undef, $s_feature_id, $s_dbxref_id);    

    unless (defined $s_feature_dbxref_id){

#       my $s_sql_insert_feature_dbxref_id=$dbh->prepare('INSERT INTO feature_dbxref (feature_id, dbxref_id) VALUES (?,?)');

        unless($s_sql_insert_feature_dbxref_id->execute($s_feature_id, $s_dbxref_id)){

            ## error inserting
            $errflag = 1;
   
            print DROP "$s_gene_id\t$s_feature_id\t$s_dbxref\t$s_dbxref_id\tERROR_INSERT_FEATURE_DBXREF\n";


            next Geneid; # Go to the beginning of the loop (outside subroutine)
        }

        ## Log insertion
        print OUT "$s_gene_id\t$s_feature_id\t$s_dbxref\t$s_dbxref_id\tADDED_FEATURE_DBXREF\n";

        ## Get list of pep feature_ids to add history tags
        print HIS "$s_feature_id\n";
    }

}
