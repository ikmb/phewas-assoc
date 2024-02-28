//import processes
include { option_check } from '../modules/option_check.nf'
include { prefilter } from '../modules/prefilter.nf'
include { generate_pcs } from '../modules/flashpca.nf'
include { gen_r2_list } from '../modules/bcftools.nf'

include {	plink2_assoc;
			plink2_assoc_merge} from '../modules/plink2_assoc.nf'

include {	make_covars;
			make_covars_nopca } from '../modules/make_covars.nf'

include {	prune;
			merge_plink;
			make_plink } from '../modules/plink2.nf'

include {	merge_plink_results;
			merge_r2
		} from '../modules/merge_processes.nf'
		
include {	regenie_step1;
			regenie_step2;
			phenofile_from_fam;
			split_input_phenofile;
			awk_regenie
		} from '../modules/regenie.nf'

include {	saige_step1;
			saige_step2} from '../modules/saige.nf'

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

	//Covars channel, content are all covars columns (PC and params.more_covars_cols) as a single comma separated string 
	ch_covars_cols = Channel.of( covars_columns() )

	for_saige_imp = Channel.fromPath(params.input_imputed_glob, checkIfExists: true).map { it ->
		def match = get_file_details(it)
		[it, match[1]] }


	//TODO: Check for sanity of phenofile or fam file regarding phenotypes
	if(params.fam){
		ch_fam_pheno = Channel.fromPath(params.fam, checkIfExists: true ).ifEmpty { exit 1, "Cannot find fam file"}
	}else{
		ch_fam_pheno = Channel.fromPath(params.phenofile, checkIfExists: true ).ifEmpty { exit 1, "Cannot find phenofile file, please specify at least one of '--phenofile' or '--fam' in your pipeline call"}
	}

	prefilter( option_check.out.readystate,
				for_saige_imp,
				ch_fam_pheno )

	ch_mapped_prefilter = prefilter.out.map { it -> [it[0], it[1], get_chromosome_code(it[0]), it[2]] }
    //.filter { it[2].toInteger() <= 22 }

	//extract_dosage( ch_mapped_prefilter )

	gen_r2_list( ch_mapped_prefilter )

	merge_r2( gen_r2_list.out.collect() )

	make_plink ( ch_mapped_prefilter,
					ch_fam_pheno )

	merge_plink ( make_plink.out.collect() )

	//removed python part from the original prune process
	prune_python_helper(merge_plink.out,
						merge_r2.out )
	prune(	merge_plink.out,
			merge_r2.out,
			prune_python_helper.out )

//FLASHPCA2
//TODO: Make PC generation optional if one wants to supply them via covariates file
	if(params.pca_dims != 0){
		generate_pcs( prune.out )

		make_covars( generate_pcs.out,
					ch_fam_pheno )

		ch_covars = make_covars.out.covars
	}else{ //NOT WORKING YET
		make_covars_nopca( option_check.out.readystate,
							ch_fam_pheno )
		
		ch_covars = make_covars_nopca.out.covars
	}

//Phenofile creation
	if(params.phenofile){
		//TODO: completely remove fam file requirement if phenofile is given
		//TODO: check if phenofile exist
		//ch_pheno = Channel.fromPath(params.phenofile, checkIfExists: true ).ifEmpty { exit 1, "Cannot find phenofile"}
		split_input_phenofile( option_check.out.readystate,
								Channel.fromPath(params.phenofile, checkIfExists: true ).ifEmpty { exit 1, "Cannot find phenofile"} )
		
		if (params.trait == "binary" ){
			ch_pheno = split_input_phenofile.out
				.flatten()
				.map{it ->
				def meta = [:]
				meta.valid = check_pheno_for_assoc(it)
				if(!meta.valid){
					println "Phenotype ${it.simpleName} contains fewer than 10 cases or controls!"
				}
			return [meta, it]
			}
		}else{
			ch_pheno = split_input_phenofile.out.flatten().map{ it -> 
				def meta = [:]
				meta.valid = true
				return [meta, it]
				}
		}
	}else{
		phenofile_from_fam( option_check.out.readystate,
							Channel.fromPath(params.fam, checkIfExists: true ).ifEmpty { exit 1, "Cannot find fam file"} )
		if (params.trait == "binary" ){
			ch_pheno = phenofile_from_fam.out
				.flatten()
				.map{it ->
				def meta = [:]
				meta.valid = check_pheno_for_assoc(it)
				if(!meta.valid){
					println "The phenotype contains fewer than 10 cases or controls!"
				}
			return [meta, it]
			}
		}else{
			ch_pheno = phenofile_from_fam.out
				.flatten()
				.map{ it -> 
				def meta = [:]
				meta.valid = true
				return [meta, it]
				}
		}
	}

//REGENIE
//TODO: Implement to remove missing phenotypes for each run separately:
	if(!params.remove_missing_phenotypes){
		ch_regenie_plink = merge_plink.out
	}else{
		ch_regenie_plink = merge_plink.out
	}
	if(!params.disable_regenie){
		ch_regenie1_input = prune.out.combine(ch_covars).combine(ch_pheno)
		//Regenie step1 should be run with less than 1mio SNPs, therefor we use the pruned plink-files
		regenie_step1( ch_regenie1_input )

		ch_regenie2_input = ch_regenie_plink.combine(regenie_step1.out)
		regenie_step2( ch_regenie2_input )
		awk_regenie( regenie_step2.out.sumstat )
	}

	if(params.saige){
		ch_saige1_input = prune.out.combine(ch_covars).combine(ch_pheno)
		//Regenie step1 should be run with less than 1mio SNPs, therefor we use the pruned plink-files
		saige_step1( ch_saige1_input )

		ch_saige2_input = saige_step1.out
		saige_step2( ch_saige2_input )
	}

	if(params.plink_assoc){
		//extract_dosage( prefilter.out.map { it -> [it[0], it[1], get_chromosome_code(it[0]), it[2]] } )
		//TODO: aktuell nur für fam, nicht für multiple phenotypes, benötigt noch einen prozess, der automatisiert neue fam files für jeden phenotype erstellt und besser zu plink2 wechseln
		plink2_assoc( 	ch_pheno.combine( ch_mapped_prefilter ).combine( ch_covars ),
						ch_fam_pheno )

		plink2_assoc_merge( plink2_assoc.out.plinksumstats.groupTuple() )
	}
}

def check_pheno_for_assoc(file) {
    def lines = file.readLines()
        int x = 0
        int y = 0	
    lines.each { String line ->
        if(line.split(" |\t")[2] == "2") x++
        if(line.split(" |\t")[2] == "1") y++
    }
    return(x > 10 && y > 10)
}

def covars_columns() {
    def pcaDims = params.pca_dims ? (1..params.pca_dims) : 0
    def moreCovars = params.more_covars_cols ? params.more_covars_cols : ''
    
    if (pcaDims == 0) {
        return moreCovars
    } else if (pcaDims instanceof List && pcaDims.every { it =~ /^\d+$/ }) {
        def pcaString = pcaDims.collect { "PC$it" }.join(',')
        if (moreCovars) {
            return "${pcaString},${moreCovars}"
        } else {
            return pcaString
        }
    } else {
        return ''
    }
}