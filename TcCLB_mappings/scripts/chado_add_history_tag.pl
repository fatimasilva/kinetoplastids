#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use DBI;
use IO::Prompt;

#### Adds history tags to a list of peptide feature IDs
#### Needs the list of peptide feature IDs (one column), the  history tag to add ($annot)
#### And type of annotation (other, structural, annotation)
#### the date (input) and user (input, same as dbuser)
#### Looks for the lowest rank available to insert the annotation (feature_cvterm)

my ($f_ids, $name_out, $help, $dbname, $dbuser, $dbport, $dbhost, $date);

GetOptions(
    '-a|ids_file=s'        => \$f_ids,
    '-o|output_prefix=s'   => \$name_out,
    '-h|help=s'            => \$help,
    '-d|dbname=s'          => \$dbname,
    '-s|dbhost=s'          => \$dbhost,
    '-p|dbport=s'          => \$dbport,
    '-u|dbuser=s'          => \$dbuser,
    '-y|date=s'            => \$date,
);

( ($f_ids && $name_out && $dbname && $dbuser && $dbport && $dbhost && $date) && !$help) || die <<USAGE;

Usage: $0
    -a|ids_file
    -o|name_out <prefix for output files>
    -h|help    <This help message>
    -d|database_name
    -u|database_user
    -s|database_host
    -p|database_port
    -y|date YYYYMMDD
USAGE

# log and out files
my $f_drop = $name_out . '.log';
my $f_out = $name_out . '_added_history.tsv';
my $f_extracted_ids = $name_out . '_extracted_ids.txt';

open DROP, ">>$f_drop";

#### Get list of IDs
my $pep_fids_ref = get_id_list($f_ids);
my @pep_fids = @{$pep_fids_ref};

### print list of ids (sanity check)
#open IDS, ">>$f_extracted_ids";
#foreach my $id (@ids){
#    print IDS "$id\n";
#}
#close IDS;

#### Connect to database
## Prompt user for db password
my $dbpass = prompt('Password:', -echo => '*');

$dbpass = "$dbpass";

## Connection details
my $dbi_connect = "DBI:Pg:dbname=$dbname;host=$dbhost;port=$dbport";
my $dbh = DBI->connect($dbi_connect, $dbuser, $dbpass,{RaiseError => 0, AutoCommit =>0}) or die "Can't connect to the database!\n";
print "Success connecting to database $dbname!\n";

## Some db values
##!00 Change annotation if needed
my $annot = 'added_dbxref';
my $errflag;

#### Start processing id list
open OUT, ">>$f_out";
Pepid: foreach my $pep_fid (@pep_fids){

    ## Check if errors inserting
    if (defined $errflag) {
        warn "Error inserting: " . $DBI::errstr;
        $dbh->rollback();
        # stop processing IDs
        last;
    };

    #### Get feature_id
    my $uniquename = get_uniquename($pep_fid);    

    if (defined $uniquename){

        print "$uniquename\t$pep_fid\n";

        #### Insert annotation 
        my $feature_cvterm_id = insert_annotation($uniquename, $pep_fid);

        ## Error inserting?
        if (defined $errflag){ warn "Error inserting annotation: " . $DBI::errstr; next Pepid};

        #### Insert qualifier, date and curator name
        insert_annotation_values($feature_cvterm_id, $uniquename, $annot, $date, $dbuser);

        ## Error inserting?
        if (defined $errflag){ warn "Error inserting annotation values: " . $DBI::errstr; next Pepid};

    } else{ # failed getting uniquename

        print DROP "$uniquename\tID_NOT_FOUND\n";
    }
}

close OUT;
close DROP;

#### Commit changes and close db connection
## uncomment when ready
$dbh->commit() unless (defined($errflag));

## close connection
$dbh->disconnect();

############
sub get_id_list {
    my $s_file = shift @_;
    my @s_ids = ();

    #### Extract IDs from file
    open FIN, "< $s_file";
    while (<FIN>){
        chomp;
        my $s_id = $_;
        push (@s_ids, $s_id);
    }

    close FIN;
    #### return list of ids
    return \@s_ids;
}

####
sub get_uniquename{
    my $s_pep_fid = shift @_;

    ## Prepare statement
    my $s_sql_get_uniquename=$dbh->prepare('SELECT uniquename FROM feature where feature_id=?');

    my $s_uniquename=$dbh->selectrow_array($s_sql_get_uniquename,undef,$s_pep_fid);

    #### Return feature_id
    return $s_uniquename;

}

####
sub insert_annotation{
    my ($s_uniquename, $s_pep_fid) = @_;
    my $s_cvterm_id = 86838; # annotation
    #my $s_cvterm_id = 86839; # structural
    #my $s_cvterm_id = 86840; # other

    #### get the lowest available rank. They do have ranks... need to take that into account to do the inserts
    my $s_rank = get_fcvt_rank($s_pep_fid);

    #### prepare insert statement
    my $s_sql_insert_annotation=$dbh->prepare('INSERT INTO feature_cvterm (feature_id, cvterm_id, pub_id, rank) VALUES (?,?,1,?)');

    unless($s_sql_insert_annotation->execute($s_pep_fid, $s_cvterm_id, $s_rank)){

        # error inserting
        $errflag = 1;
    }
    
    #### success inserting? get new feature_cvterm_id
    my $s_fcvt_id;

    if (defined $errflag) {

        warn "Error loading new annotation into feature_cvterm: " . $DBI::errstr;
        $dbh->rollback();
        # stop processing file
        last;

    }else{ # get new feature_cvterm_id and log

        my $s_sql_get_fcvt_id=$dbh->prepare('SELECT feature_cvterm_id FROM feature_cvterm WHERE feature_id=? AND cvterm_id=? AND rank=?');

        my $s_fcvt_id_ref = $dbh->selectall_arrayref($s_sql_get_fcvt_id,undef,$s_pep_fid,$s_cvterm_id,$s_rank);

        my $s_row_ref = @{$s_fcvt_id_ref}[0]; # get only first value (most recent)
        $s_fcvt_id = @{$s_row_ref}[0]; # get fcvt_id of the most recent annotation change

        print "added_annotation\t$s_pep_fid\t$s_uniquename\t$s_fcvt_id\t$s_rank\n";
        print OUT "$s_pep_fid\t$s_uniquename\t$s_fcvt_id\t$s_rank\tADDED_ANNOTATION\n";

    }
    #### return new feature_cvterm_id
    return $s_fcvt_id;
}

####
sub get_fcvt_rank {
    my $s_pep_fid = shift @_;
    my $s_cvterm_id = 86838; # annotation
    #my $s_cvterm_id = 86839; # structural
    #my $s_cvterm_id = 86840; # other

    #### Get annotation rank (needed for insert but they are not ordered in db)
    my $s_sql_rank=$dbh->prepare('SELECT rank FROM feature_cvterm WHERE feature_id=? AND cvterm_id=? ORDER BY rank DESC');

    my $s_higher_rank=$dbh->selectrow_array($s_sql_rank,undef,$s_pep_fid, $s_cvterm_id);

    my $s_rank = 0;

    if (defined $s_higher_rank){
        $s_rank = $s_higher_rank + 1;

    } else { # (rank =0 ) No annotation history in that gene yet

    }

    #### return next rank available to insert
    return $s_rank;

}

####
sub insert_annotation_values {
    my ($s_fcvt_id, $s_uniquename, $s_annot, $s_date, $s_curator) = @_;
    my $s_t_id_qual = 26760;
    my $s_t_id_date = 1697;
    my $s_t_id_curator = 26774;

    #### insert qualifier
    my $s_sql_avalue=$dbh->prepare('INSERT INTO feature_cvtermprop (feature_cvterm_id, type_id, value) VALUES (?,?,?)');
    unless ($s_sql_avalue->execute($s_fcvt_id, $s_t_id_qual, $s_annot)){

        # error insertng
        $errflag = 1;
     
    }
    #### insert date
    unless ($s_sql_avalue->execute($s_fcvt_id, $s_t_id_date, $s_date)){

        # error inserting
        $errflag = 1;
    }
    #### insert curator email
    unless ($s_sql_avalue->execute($s_fcvt_id, $s_t_id_curator, $s_curator)){

        # error inserting
        $errflag = 1;
    }

    #### log if all suscessful
    if (defined $errflag) {

        warn "Error loading new values into feature_cvtermprop: " . $DBI::errstr;
        $dbh->rollback();
        # stop processing file
        last;

    }else{ # log

        print "added_avalues\t$s_uniquename\t$s_fcvt_id\t$s_annot\t$s_date\t$s_curator\n";
        print OUT "$s_uniquename\t$s_fcvt_id\t$s_annot\t$s_date\t$s_curator\tADDED_AVALUES\n";

    }


}
