class WorkflowMain {
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        
        log.info header(workflow)

        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

    }

    public static String header(workflow) {
        def headr = ''
        def info_line = "IKMB PheWAS association testing pipeline | version ${workflow.manifest.version}"
        headr = """
    ===============================================================================
    ${info_line}
    ===============================================================================
    """
        return headr
    }

    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} --input_imputed_glob /path/to/QCed/'*.noATCG.vcf.gz' --fam /path/to/SNPQCII/Example_QCed.fam -profile standard"
        def help_string = ''
        // Help message
        help_string = """
Usage: nextflow run ikmb/phewas-assoc --input_imputed_glob "/path/to/QCed/"'*.noATCG.vcf.gz' --fam /path/to/SNPQCII/Example_QCed.fam --collection_name Example --output putput/directory

This example will run an association test with regenie with the phenotype supplied in the fam file on the chromosome-wise splitted vcf files that are the output of ikmb/gwas-qc.

Required parameters:
    --input                        A glob of chromosome-wise split vcf.gz files
    --fam                          The corresponding fam file to the input. 
                                   Can contain fewer samples than the vcf, these samples will be removed from analysis

Optional parameters:
    --phenofile                    Phenotype file for multiple phenotype/traits-testing with regenie. 
                                   Tab separated file with columnsheader "FID IID Phenotype1 Phenotype2" 
                                   Entries must be "0" for FID, "FID_IID" for IID and all phenotypes must be 
                                   either binary or quantitaive, don't mix! Missing Samples will be ignored. 
                                   Binary traits should be specified as control=1,case=2,missing=NA.

    --test                         Test algorithm. May be 'firth' (default) or 'spa'.
    --more_covars                  Whitespace-separated list of covariates with columnsheader "FID IID Covar1 Covar2" 
                                   Entries must be "0" for FID, "FID_IID" for IID

    --more_covars_cols             Comma-separated list of covar column header names. Example: "Covar1,Covar2"
    --null_filter [filter]         Bcftools-style formatted INFO filter for generation of the null model. Default: "R2>0.8"
    --additional_regenie_parameter Add additional parameters to step2 of regenie e.g. annotation and mask parameters for gene-based testing.
    --build                        37 or 38 are valid options. No use currently.
Output:
    --output                       Local directory to which all output is written
    --collection_name              Output filename prefix
        """
        return help_string
    }
}