# Usage information

## Quick Start
1. Install and make sure, Singularity (now Apptainer) and Nextflow are working. For example via [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html):
    ```bash
    # Create a new conda environment for Singularity and Nextflow
    conda create --name nf_env -c conda-forge -c bioconda singularity nextflow
    # Activate environment
    conda activate nf_env
    # Try a simple Nextflow demo
    nextflow run hello
    # Check Singularity
    singularity --version
    ```
   
2. Run the [Quality Control Pipeline](https://github.com/ikmb/gwas-qc/blob/master/Readme.md#quick-start) on the example first.
    - All files necessary for the association testing pipeline are automatically generated.

3. Run the phewas-assoc pipeline with like so:
   ```
   nextflow run ikmb/phewas-assoc \
    --input_imputed_glob "/pathto/gwas-qc/example/output/Example/QCed/"'*.noATCG.vcf.gz'" \
    --fam "/pathto/gwas-qc/example/output/Example/SNPQCII/Example_QCed.fam" \
    --collection_name "EXAMPLE" \
    --output "output/Example/Assoc"
    ```

## How to Start

You will need to:
- A set of `.vcf.gz` files with the following specifics:
    - at most one chromosome per `.vcf.gz` file. Multiple chromosomes per files are not supported. If a file happens to have multiple chromosomes, only the first will be analyzed.
    - any chromosome codes are supported (i.e. `chrX`, `X`, `23`, `chr23`, `chromosomeX` are just fine)
    - **your VCF files can contain genotype probabilities (GP tag) for imputed data or genotype calls (GT tag) for genotyped data for the analysis with Plink2. If both tags are found, GP is chosen.**
    - the `INFO` column in the VCF files should contain an imputation score. This is used to filter the input variants to create a good null model for SAIGE. For topmed imputations, we found `R2>0.8` to yield good results. If the given tag is not found in the VCF file, the default value 1.0 is assumed, making all variants pass the filter. This filter is not used for Plink-based association testing. 
- A FAM file to update sex and phenotype from. Only those files in the FAM file will be used from the VCFs. Note that *all* samples from the FAM file must be present in the VCF files. **Note: if the FAM family ID is 0, the sample name in the VCF should be the individual ID. If the FID is not 0, the sample name should be the family ID and the individual ID with an underscore (i.e. FID 1234 IID 9876 should be 1234_9876 in the VCF file)**
- The genome build of the input data, which will result in Regenie handling the sex chromosomes according to the coordinates of the pseudoautosomal regions. Possible values are 37 and 38.
- Optionally, you can specify additional covariates to be used in association testing. By default, 10 principal components are automatically generated and used. If you want additional covariates, have a tab-separated file with a header at hand. The first two columns should be the family and individual ID in `FID\tIID` format, any futher column is treated as a covariate. Specify the covars file with `--more_covars $FILE` and the columns to be used with `--more_covars_cols AGE,SEX`, where `AGE` and `SEX` are the respective column headers from the covar file that you wish to be included.
- The pipeline output and reports will be written to the `--output` directory.

### Parameters

The following list covers all parameters that may be specified for the Association Testing Pipeline:


| Category | Command | Type  | Description |
| --- | --- | --- |  --- |
| [REQUIRED] | `--input_imputed_glob` | [glob]  | A glob expression to specify the .vcf.gz files that should be used for association analysis |
| [REQUIRED] | `--fam` | [filepath]  | A Plink-style FAM file that will be used to select a subset of samples from the provided VCFs |
| [OPTIONAL] | `--phenofile` | [filepath]  | Phenotype file for multiple phenotype/traits-testing with regenie. Tab separated file with columnsheader "FID IID Phenotype1 Phenotype2". Entries must be "0" for FID, "FID_IID" for IID and all phenotypes must be either binary or quantitaive, don't mix! Missing Samples will be ignored. Binary traits should be specified as control=1,case=2,missing=NA. |
| [ADVISED] | `--trait` | [type]  | Trait type to analyze. May be 'binary' (default) or 'quantitative'. For a binary trait use "1" as control and "2" as case in the phenofile/fam. |
| [ADVISED] | `--test` | [type]  | Test algorithm. May be 'firth' (default) or 'spa'. |
| [ADVISED] | `--build` | [integer] | Define the human genome build code. Valid numbers are 37 and 38. |
| [ADVISED] | `--collection_name` | [name] | Output filename prefix |
| [ADVISED] | `--output` | [filepath]  | Output directory. Default: "output/assoc" |
| [OPTIONAL] | `--more_covars` | [filepath] | Whitespace-separated list of covariates. Columnsheader "FID IID Covariate1 Covariate2". |
| [OPTIONAL] | `--more_covars_cols` | [string] | A comma-separated list of covar column header names that should be used from the file that is used with `--more_covars`. Required when `--more_covars` is being used |
| [OPTIONAL] | `--null_filter` | [string] | A bcftools-style formatted INFO filter for generation of the null model. Default: "R2>0.8" |
| [OPTIONAL] | `--additional_bcftools_arg` | [string] | A bcftools-style formatted view filter for removing variants from the association test, e.g. MAF filter with "-q 0.01:minor". Default: "" |
| [OPTIONAL] | `--additional_regenie_parameter` | [string] | Add additional parameters to step2 of regenie e.g. annotation and mask parameters for gene-based testing. |
| [OPTIONAL] | `--pca_dims` | [integer] | Define the limit of how many PCs should be calculated and included in association testing. Expects integer values. 0 would mean, no PCs will be calculated. Default is 10. |
| [OPTIONAL] | `--plink_assoc` |  | Activation-switch to also perform association tests with plink2 --glm. |
| [OPTIONAL] | `--plink2_glm_options` | [string] | When performing plink2 association testing, adjust the --glm parameter within plink2 with modifiers. Default: "omit-ref hide-covar". |
| [OPTIONAL] | `--disable_regenie` |  | Deactivation-switch to deactivate association test calculation with regenie. |
| [EXPERIMENTAL] | `--regenie_step1_input` |  | Use a different plink fileset (bim/bed/fam) as input for regenie in step1. Default: false. Input must be the path to the fileset prefix (e.g. for prefix.bed + prefix.bim + prefix.fam it would be "--regenie_step1_input prefix") |
| [EXPERIMENTAL] | `--saige` |  | Enable Saige for association testing. |
| [OPTIONAL] | `-resume` | | Activation-switch to restart where the pipeline was when cancelled or aborted. May or may not work, depending on your filesystem specifics. |


## Quick start on Kiel medcluster
On medcluster you only need to load the dependencies Singularity and Nextflow. Then the pipeline can directly be executed.

   ```
   module load singularity nextflow

   nextflow run ikmb/phewas-assoc \
    --input_imputed_glob "gwas-qc/example/output/Example/QCed/"'*.noATCG.vcf.gz'" \
    --fam "gwas-qc/example/output/Example/SNPQCII/Example_QCed.fam" \
    --collection_name "EXAMPLE" \
    --output "output/Example/Assoc"   
   ```
