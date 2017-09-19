#!/usr/bin/perl -w

############ DECLARE MODULES ##################

use strict;
use Data::Dumper;
use lib "REST/lib/perl5/site_perl/REST";
use lib "JSON/lib/perl5/site_perl";
use Client;
use JSON;

########### DECLARE MAIN VARIABLES ############

my ($line, $gene_details, $type, $depth, $refcount, 
    $altcount, $refpercentage, $altpercentage, 
    $exac_return);

my @vcf_fields;

my $print_header_yes = 1;

########### MAIN ########################################

open(FH, "Challenge_data.vcf") or die "cannot open Challenge_data file for reading\n";

while($line=<FH>){
	
	#skip comments in VCF
	next if $line=~/^#/;
	chomp $line;

	#split fields by tab
	@vcf_fields = split /\t/, $line;
	
	$gene_details = construct_url_rest(($vcf_fields[0], $vcf_fields[1], $vcf_fields[3], $vcf_fields[4]));
	$type = get_most_deleterious_annotation($vcf_fields[7]);
	$depth = get_depth($vcf_fields[7]);
	$refcount = get_refcount($vcf_fields[7]);
	$altcount = $depth - $refcount;
	$refpercentage = $refcount/$depth;
	$altpercentage = 1-$refpercentage;
	$exac_return = call_client($gene_details);
	print_output($type, $depth, $altcount, $altpercentage, $refpercentage, $exac_return, $print_header_yes);
	$print_header_yes+=1;	

}
close FH;	

########## CONSTRUCT QUERY FOR EXAC DB #####################

sub construct_url_rest{
	
	my ($chromosome, $position, $ref, $alt) = @_;
        return($chromosome."-". $position."-". $ref. "-".$alt);
}

######### GET MOST DELETERIOUS ANNOTATION #####################

sub get_most_deleterious_annotation{
	
	my $annotype = (split /=/,$_[0])[-1];
	my @anno = split /,/, $annotype;
	if(scalar(@anno)==1){
		return $anno[0];
	}
	else{
		if (grep(/del/, @anno) && (grep /ins/, @anno)) {
			return "indel";
		}
		elsif(grep(/del/, @anno)){
			return "del";
		}
		elsif(grep(/ins/, @anno)){
			return "ins";
		}
		elsif(grep(/complex/, @anno)){
			return "complex";
		}
		elsif(grep(/mnp/, @anno)){
			return "mnp";
		}
		else{
			return "snp";
		}
	}

}

########### GET READ DEPTH ###################################

sub get_depth{

	my $getdepth = (split /;/, $_[0])[7];
        return ((split /=/, $getdepth)[1]);

}

########### GET REFERENCE READ COUNT ############################

sub get_refcount{

	my $getref = (split /;/, $_[0])[28];
	return ((split /=/, $getref)[1]);

}

########### EXAC CONNECT ########################################

sub call_client{
	my @client_return = ();
	my $gene = shift;
	my $client = REST::Client->new();
	$client->setHost("http://exac.hms.harvard.edu");
	$client->GET(
    	"/rest/variant/$gene"
	);
	
	my $response = from_json($client->responseContent());
	my $allele_freq = $response->{'variant'}->{'allele_freq'};
	
	# Set 0 for no allele frequency info

	if(!($allele_freq)){
		$allele_freq = 0;
	}

	my $variant_id = $response->{'variant'}->{'variant_id'};

	#Set NA for no variant id info from EXAC
	
	if(!($variant_id)){
                $variant_id = "No Exac Info";
        }
	
	my $variant_in_gene = $response->{'variant'}->{'genes'};

	if(!($variant_in_gene)){
		$variant_in_gene = ["No Exac Info"] ;
	}
	
	push(@client_return, $allele_freq);
	push(@client_return, $variant_id);
	push(@client_return, $variant_in_gene);

	return(\@client_return);
}

############## PRINT TOP 10 GENE INFO ###########################

sub print_gene_info{

	my $db_ref = shift;
        my $i; 

	for($i=0; $i<10; $i++){
                if($db_ref->[2]->[$i]){
                        print $db_ref->[2]->[$i];
                	print "\t";
		}
	}
	print "\n";
	return;
}

############# PRINT OUTPUT ######################################

sub print_output{
 	my ($type, $depth, $altcount,$altpercentage,$refpercentage, $exac_return, $print_header_yes) = @_;
	
	if($print_header_yes == 1){
		
		printf("%.8s\t|\t%.8s\t|\t%.8s\t|\t%.8s\t|\t%.8s\t|\t%.8s\t|\t%.50s\t|\t%.100s\n", "----", "-----", "-----", "----", "----", "-------", "-----------", "---------------");
		printf("%.8s\t|\t%.8s\t|\t%.8s\t|\t%.8s\t|\t%.8s\t|\t%.8s\t|\t%.50s\t|\t%.100s\n", "TYPE", "DEPTH", "\#VARS", "\%VAR", "\%REF", "EXAC_AF", "EXAC_VAR_ID", "GENE ANNOTATION");
		printf("%.8s\t|\t%.8s\t|\t%.8s\t|\t%.8s\t|\t%.8s\t|\t%.8s\t|\t%.50s\t|\t%.100s\n", "----", "-----", "-----", "----", "----", "-------", "-----------", "---------------");
		
	}	

	printf("%.8s\t|\t%.8s\t|\t%.8s\t|\t%.3f\t|\t%.3f\t|\t%.3f\t|\t%.50s\t|\t", $type, $depth, $altcount, $altpercentage, $refpercentage, $exac_return->[0], $exac_return->[1]); 
	print_gene_info($exac_return);
 
}
