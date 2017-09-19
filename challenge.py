############## IMPORT REQUIRED MODULES #########################

import requests
import re

###############  MAIN ##########################################

def main():
	
	print_header_yes = 1
	f = open("Challenge_data.vcf", "r");
	for line in f:
		#skip comments in VCF
		if line.startswith('#'):
			continue
		line = line.rstrip('\n');
		#split VCF fields by tab
		vcf_fields = line.split("\t")

		gene_details = construct_url_rest(vcf_fields[0], vcf_fields[1], vcf_fields[3], vcf_fields[4])
		type = get_most_deleterious_annotation(vcf_fields[7])
		depth = get_depth(vcf_fields[7])
		refcount = get_refcount(vcf_fields[7])	
		altcount = (depth - refcount)	
		refpercentage = refcount/depth
		altpercentage = 1-refpercentage
		exac_return = call_client(gene_details)
		print_output(type, depth, altcount, altpercentage, refpercentage, exac_return['allele_freq'], exac_return['variant_id'], exac_return['variant_in_gene'], print_header_yes);
		print_header_yes += 1


############### CONSTRUCT QUERY FOR EXAC DB #####################


def construct_url_rest(chromosome, position, ref, alt):
        return("-". join([chromosome, position, ref,alt]))


############### GET MOST DELETERIOUS ANNOTATION ##################

def get_most_deleterious_annotation(annotype):

        annotype = re.split("=",annotype)[-1];
        anno = []
	anno = re.split(",", annotype)
        if len(anno)==1:
		return anno[0]
        else:
                if filter(lambda x:'del' in x, anno) and filter(lambda x:'ins' in x, anno):
			return "indel"
                elif filter(lambda x:'del' in x, anno):
                        return "del"
                elif filter(lambda x:'ins' in x, anno):
                        return "ins"
                elif filter(lambda x:'complex' in x, anno):
                        return "complex"
                elif filter(lambda x:'mnp' in x, anno):
                        return "mnp";
                else:
			return "snp"


################## GET READ DEPTH ##############################
                
def get_depth(getdepth):

        getdepth = re.split(";", getdepth)[7]
        return int(re.split("=", getdepth)[1])


################ GET REFERENCE READ COUNTS #######################


def get_refcount(getrefcount):
	getrefcount = re.split(";", getrefcount)[28]
        return int(re.split("=", getrefcount)[1])


################# EXAC CONNECT ##################################

def call_client(gene_desc):
	result = requests.get("http://exac.hms.harvard.edu/rest/variant/%s" %  gene_desc)
        if result.status_code!=200:
                print 'No response for gene_desc';
        json_data = result.json()
	if 'allele_freq' in json_data['variant']:
		allele_freq = json_data['variant']['allele_freq']
	else:
		allele_freq = 0
	if 'variant_id' in json_data['variant']:
		variant_id = json_data['variant']['variant_id']
	else:
		variant_id = "No Exac Info"
	if 'genes' in json_data['variant']:
		variant_in_gene = json_data['variant']['genes']
	else:
		variant_in_gene = ["No Exac Info"]
	return {'allele_freq': allele_freq, 'variant_id': variant_id, 'variant_in_gene': variant_in_gene}



############### PRINT OUTPUT #####################################

def print_output(type, depth, altcount, altpercentage, refpercentage, allele_freq, variant_id, genes, print_header_yes):
	if print_header_yes==1:
		fmt = "{type:s}\t|\t{depth:s}\t|\t{altcount:s}\t|\t{altpercentage:s}\t|\t{refpercentage:s}\t|\t{allele_freq:s}\t|\t{variant_id:s}\t|\t{genes:s}"
		print fmt.format(type="----", depth="-----", altcount="-----", altpercentage="-----", refpercentage="-----", allele_freq="-------", variant_id="-----------", genes="---------------")
		print fmt.format(type="TYPE", depth="DEPTH", altcount="#VARS", altpercentage="%VARS", refpercentage="%REFS", allele_freq="EXAC_AF", variant_id="EXAC_VAR_ID", genes="GENE ANNOTATION")
		print fmt.format(type="----", depth="-----", altcount="-----", altpercentage="-----", refpercentage="-----", allele_freq="-------", variant_id="-----------", genes="---------------")
	else:
		fmt = "{type:s}\t|\t{depth:d}\t|\t{altcount:d}\t|\t{altpercentage:0.3f}\t|\t{refpercentage:0.3f}\t|\t{allele_freq:0.3f}\t|\t{variant_id:s}\t|\t{genes:s}"
		newgenes = "\t". join(map(str,genes))	
		print fmt.format(type=type, depth=depth, altcount=altcount, altpercentage=altpercentage, refpercentage=refpercentage, allele_freq=allele_freq, variant_id=variant_id, genes=newgenes)

if __name__ == '__main__':
    main()
