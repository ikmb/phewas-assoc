#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/**
===============================
IKMB PheWAS Assoc Pipeline
===============================

This Pipeline performs Association Testing for GWAS, mGWAS and PheWAS on QC'ed (and imputed) chromosome-wise genotyping data in vcf.gz format.

### Homepage / git
git@github.com:ikmb/phewas-assoc.git

**/

// Pipeline version
params.version = workflow.manifest.version

//Help page & Pipeline Header
WorkflowMain.initialise(workflow, params, log)

//Input Parameter Checks
WorkflowPhewas.initialise( params, log)

//If run_name wasn't set by the user, set it to workflow runName.
if(!params.hasProperty('run_name')){
  params.run_name = workflow.runName
}

log.info  "Run_name: '${params.run_name}'\n"

//Import workflow
include { assoc } from './workflows/Assoc' params(params)

//execute workflow
workflow {

	assoc()

}

//completion logs
workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}

