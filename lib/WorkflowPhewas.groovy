class WorkflowPhewas {
    public static void initialise(params, log) {

        //genomeExistsError(params, log)
/*
        if (!params.hasProperty('run_name')) {
            log.info  "No run name specified, creating automatically one"
        }
*/
        if (!params.fam) {
            log.info "No fam file specified! Please specify with --fam ! Exiting now."
            System.exit(1)
        }

        if (!FileValidator.validateFile(params.fam)) {
            log.info "Fam File '${params.fam}' does not exist or is empty, exiting now."
            System.exit(1)
        }

        if (!(params.test == "spa" || params.test == "firth")) {
            log.info  "Parameter 'test' can only be set to 'spa' or 'firth', exiting now."
            System.exit(1)
        }

        if (!(params.trait == "binary" || params.trait == "quantitative")) {
            log.info  "Parameter 'trait' can only be set to 'binary' or 'quantitative', exiting now."
            System.exit(1)
        }

        if (!(params.build == 37 || params.build == 38)) {
            log.info  "Parameter 'build' can only be set to '37' or '38', exiting now."
            System.exit(1)
        }

        if(params.phenofile){
            if (!FileValidator.validateFile(params.phenofile)) {
                log.info "Phenofile '${params.phenofile}' does not exist or is empty, exiting now."
                System.exit(1)
            }
            def rows = FileValidator.count_file_lines(params.phenofile)
            log.info "Phenofile contains ${rows} rows"
        }

        if(params.regenie_step1_input){
            def bed = params.regenie_step1_input + '.bed'
            def bim = params.regenie_step1_input + '.bim'
            def fam = params.regenie_step1_input + '.fam'
            if (!FileValidator.validateFile(${bed}) || !FileValidator.validateFile(${bim}) || !FileValidator.validateFile(${fam}) ) {
                log.info "The plink fileset with prefix '${params.phenofile}' .bim .bed and .fam  for regenie step 1 do not exist, are not complete or are empty, exiting now."
                System.exit(1)
            }
        }

//TODO: more_covars is "." as default
        if(params.more_covars && params.more_covars != "."){
            if (!FileValidator.validateFile(params.more_covars)) {
                log.info"Covariate File '${params.more_covars}' does not exist or is empty, exiting now."
                System.exit(1)
            }

            //Covariate file should contain the same amount of samples
            if(params.phenofile){
                def phenorows = FileValidator.count_file_lines(params.phenofile)
                def covarsrows = FileValidator.count_file_lines(params.more_covars)
                log.info "Covariate file contains ${covarsrows} rows"
                if(phenorows != covarsrows){
                    log.info "Phenofile and Covariate file do not contain the same amount of rows, please check these files. Exiting now."
                    System.exit(1)
                }
            }else{
                def famrows = FileValidator.count_file_lines(params.fam)
                def covarsrows = FileValidator.count_file_lines(params.more_covars)
                log.info "Fam file contains ${famrows} rows"
                log.info "Covariate file contains ${covarsrows} rows"
                if((famrows + 1) != covarsrows){
                    log.info "Fam file and Covariate file do not contain the same amount of samples, please check these files. Exiting now."
                    System.exit(1)
                }
            }
        }    

        if (!params.pca_dims == 0) {
            log.info "No PCA step will be performed."
        }
/*
        if (params.assembly == "CHM13v2") {
            log.info "WARNING!!! Use of the T2T reference is highly experimental and not supported by some of the available tools..."
        }
*/
    }
/*
    private static void genomeExistsError(params, log) {
        if (params.assembly && !params.genomes.containsKey(params.assembly)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome assembly '${params.assembly}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available assemblies  are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }
*/
}

class FileValidator {

    def static validateFile(def filePath) {
        def file = new File(filePath)
        
        if (file.exists() && file.isFile() && file.size() > 0) {
            return true  // File exists and is not empty
        } else {
            return false // File does not exist or is empty
        }
    }

    def static count_file_lines(def filePath) {
        def file = new File(filePath)

        def lines = file.readLines()
        int x = 0	
        lines.each { String line -> x++ }
        return x 
    }

}