//import processes
include { option_check } from '../modules/option_check.nf'
include { prefilter } from '../modules/prefilter.nf'
include { generate_pcs } from '../modules/flashpca.nf'
include { make_covars } from '../modules/make_covars.nf'
include { gen_r2_list } from '../modules/bcftools.nf'

include { prune;
		  merge_plink;
          make_plink } from '../modules/plink2.nf'

include { merge_plink_results;
          merge_r2
        } from '../modules/merge_processes.nf'
		
include {regenie_step1;
		 regenie_step2;
		 phenofile_from_fam;
		 split_input_phenofile
		} from '../modules/regenie.nf'

include { prune_python_helper } from '../modules/python2.nf'

//function definitions
def get_file_details(filename) {
    def m = filename =~ /\/([^\/]+)(\.vcf\.gz|\.bgen)$/
    return [ m[0][1], m[0][2] ]
	}

def get_chromosome_code(filename) {
    def m = filename =~ /\/([^\/]+).ap_prf.vcf.gz$/
    return m[0][1]
	}


//PIPELINE WORKFLOW
workflow assoc{

	main:

	option_check()
	
	//if(!params.fam){exit 1, "Cannot find fam file"}

	//params.fam_length=check_fam_for_saige( file(params.fam) )

	for_saige_imp = Channel.fromPath(params.input_imputed_glob, checkIfExists: true).map { it ->
    	def match = get_file_details(it)
    	[it, match[1]] }

	prefilter( for_saige_imp,
				Channel.fromPath(params.fam) )

	ch_mapped_prefilter = prefilter.out.map { it -> [it[0], it[1], get_chromosome_code(it[0]), it[2]] }

	//extract_dosage( ch_mapped_prefilter )

	gen_r2_list( ch_mapped_prefilter )

	merge_r2( gen_r2_list.out.collect() )

	make_plink ( ch_mapped_prefilter,
				 Channel.fromPath(params.fam) )

	merge_plink ( make_plink.out.collect() )

	//removed python part from the original prune process
	prune_python_helper(merge_plink.out,
						merge_r2.out )
	prune( merge_plink.out,
		   merge_r2.out,
		   prune_python_helper.out )

//FLASHPCA2
	generate_pcs( prune.out )

	make_covars( generate_pcs.out,
					   Channel.fromPath(params.fam, checkIfExists: true ).ifEmpty { exit 1, "Cannot find fam file"} )

//REGENIE
	if(!params.disable_regenie){
		if(params.phenofile){
			//TODO: check if phenofile exist
			//ch_pheno = Channel.fromPath(params.phenofile, checkIfExists: true ).ifEmpty { exit 1, "Cannot find phenofile"}
			split_input_phenofile( Channel.fromPath(params.phenofile, checkIfExists: true ).ifEmpty { exit 1, "Cannot find phenofile"} )
			ch_pheno = split_input_phenofile.out.flatten().view()
		}else{
			phenofile_from_fam( Channel.fromPath(params.fam, checkIfExists: true ).ifEmpty { exit 1, "Cannot find fam file"} )
			ch_pheno = phenofile_from_fam.out
		}
		//Regenie step1 should be run with less than 1mio SNPs, therefor we use the pruned plink-files
		regenie_step1( prune.out,
					   make_covars.out.covars,
					   ch_pheno )

		regenie_step2( merge_plink.out,
					   make_covars.out.covars,
					   regenie_step1.out,
					   ch_pheno )
	}
}
