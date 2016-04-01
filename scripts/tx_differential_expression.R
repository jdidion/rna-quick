#' Perform differential expression analysis using sleuth.
#' requires a data table with one row per sample. The only
#' required column is `sample`, which gives the sample name.
#' There should also be at least one additional column that
#' specifies a condition that divides the samples into at
#' at least two groups. There can also be a `path` column to
#' specify that directory that contains the kallisto output 
#' for the sample, otherwise the path is assumed to be
#' <data.dir>/<sample>.
#'
#' The default values for the full model and beta value
#' (which you should almost certainly change) assume 
#' that the model is testing a single response variable 
#' (condition) assumed to have two possible values 
#' ('case' and 'control').
#'
#' @param samples data.frame, data.table, or path to a csv,
#' tsv, or Rdata file containing the data table.
#' @param models model(s), specified as formula. If a single
#' value, taken to be the full model, otherwise a list of
#' model names and formulas. If a model is named 'full' that
#' will be used as the full model, otherwise the first model
#' in the list.
#' @param betas response variable(s) to test for differential
#' expression. Each is a concatanation of the variable name
#' and the factor relative to which beta values should be
#' calculated. If a vector, each beta will be tested for each
#' model; otherwise, a list mapping model name to beta.
#' @param data.dir parent directory containing kallisto 
#' output directories.
#' @param outfile RData file in which to save final sleuth
#' object. Defaults to "sleuth_results.RData" in `data.dir`.
#' @param db data.frame with a target_id column containing
#' the transcript IDs used in the Kallisto index and mapping
#' them to other features (e.g. gene IDs/names). Uses the
#' bioMart Ensembl database by default. Use NA to avoid using
#' transript->gene mappings.
run_sleuth <- function(samples, models=~condition, betas="conditioncontrol", 
                       tx2gene=NULL, max_bootstrap=100) {
    if (length(models) == 1) {
        names(models) <- "full"
        partial <- NULL
    }
    else {
        n <- names(models)
        if (is.null(n)) {
            n <- c("full", paste0("partial", 1:(length(models)-1)))
        }
        else if (!("full" %in% names(models))) {
            n[1] <- "full"
        }
        names(models) <- n
        partial <- setdiff(n, "full")
    }
    
    if (is.null(tx2gene)) {
        target_id <- read.table(samples[1, "path"], sep="\t", header=T, stringsAsFactors=FALSE)$target_id
        tx2gene <- data.frame(
            tx_id=stringr::str_extract(target_id, 'ENSTR?[\\d\\.]+'),
            gene_id=stringr::str_extract(target_id, 'ENSGR?[\\d\\.]+')
        )
    }
    
    message("Preparing data")
    data <- sleuth::sleuth_prep(samples, models$full target_mapping=tx2gene, max_bootstrap=max_bootstrap)
    
    for (mod in partial) {
        message(paste("Fitting model", mod))
        data <- sleuth::sleuth_fit(data, formula=models[[mod]], fit_name=mod)
    }
    
    if (is.list(betas)) {
        for (mod in names(betas)) {
            for (b in betas[[mod]]) {
                message(paste("Testing", b, "for model", mod))
                data <- sleuth::sleuth_wt(data, which_beta=beta, which_mod=mod)
            }
        }
    }
    else {
        for (mod in names(models)) {
            for (b in betas) {
                message(paste("Testing", b, "for model", mod))
                data <- sleuth::sleuth_wt(data, which_beta=beta, which_mod=mod)
            }
        }
    }

    data
}

'usage: tx_DE.R [options]

options:
 -d <dir>, --kallisto-dir <dir>     Diretory with Kallisto results [default: .]
 -m <file>, --metadata-file <file>  Metadata file path [default: libraries.txt]
 -o <file>, --output <file>         RData output file path [default: sleuth_results.RData]
 -x <file>, --tx2gene-file <file>   TSV file with mappings from transcript to gene IDs' -> doc

opts <- docopt::docopt(doc)
metadata <- read.table(opts$metadata_file, sep="\t", header=T, stringsAsFactors=F)
metadata$path <- file.path(opts$data.dir, metadata$UniqueID)

# loop through models and convert using as.formula

sleuth_result <- run_sleuth(metadata, models, betas, opts$tx2gene_file)
save(sleuth_result, file=opts$output)
