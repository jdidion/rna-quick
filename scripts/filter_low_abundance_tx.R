
library(data.table)
library(docopt)

# Return a list of length C, where C is the number of 
# experimental conditions. Each list item is a vector
# with names as sample names and values as abundance.tsv
# file paths.
get_quant_files <- function(metadata_file="libraries.txt", quant_dir=".", condition="Condition") {
    md <- read.table(metadata_file, sep="\t", header=T, stringsAsFactors=F)
    if (!(condition %in% colnames(md))) {
        stop(paste(condition, "is not a column in", metadata_file))
    }
    tapply(file.path(quant_dir, md$UniqueID, "abundance.tsv"), md[,condition], identity, simplify=FALSE)
}

# Load and sum all quant files in quant.files vectors.
# Returns a single vector of TPM values with names as
# transcript names.
sum_abundances <- function(abundance_files) {
    abundances <- NULL
    for (i in 1:length(abundance_files)) {
        sample <- names(abundance_files)[i]
        a <- read.table(quant.files[i], sep="\t", header=T, stringsAsFactors=FALSE)
        if (i == 1) {
            abundances <- a$tpm
            names(abundances) <- a$target_id
        }
        else {
            abundances += a$tpm
        }
    }
    abundances
}

filter_transcripts <- function(abundance_files, cutoff=0.1, tx2gene=NULL) {
    abundances <- data.table(do.call(cbind, lapply(abundance_files, sum_abundances)))
    if (is.null(tx2gene)) {
        abundances$tx_id <- stringr::str_extract(rownames(abundances), 'ENSTR?[\\d\\.]+')
        abundances$gene_id <- stringr::str_extract(rownames(abundances), 'ENSGR?[\\d\\.]+')
    }
    else {
        if (typeof(tx2gene) == "character") {
            tx2gene <- read.table(tx2gene, sep="\t", header=TRUE, stringsAsFactors=FALSE)
        }
        colnames(tx2gene) <- c("tx_id", "gene_id")
        m <- match(rownames(abundances), tx2gene$tx_id)
        if (any(is.na(m))) {
            stop("One or more transcript IDs missing from tx2gene mapping")
        }
        abundances <- cbind(abundances, as.data.table(tx2gene[m,]))
    }
    abundances[, rel_abundance := tpm / sum(tpm), gene_id]
    abundances[rel_abundance < cutoff, tx_id]
}

'usage: filter_low_abundance_tx.R [options]

options:
 -c <name>, --condition <name>      Name of condition used to partition samples [default: Condition]
 -d <dir>, --kallisto-dir <dir>     Diretory with Kallisto results [default: .]
 -m <file>, --metadata-file <file>  Metadata file path [default: libraries.txt]
 -o <file>, --output <file>         Output file path [default: exclude_tx_ids.txt]
 -t <float>, --cutoff <float>       Cutoff value between 0-1 [default: 0.1]
 -x <file>, --tx2gene-file <file>   TSV file with mappings from transcript to gene IDs' -> doc

opts <- docopt(doc)
quant_files <- get_quant_files(opts$metadata_file, opts$kallisto_dir, opts$condition)
exclude_tx_ids <- filter_transcripts(quant_files, as.numeric(opts$cutoff), opts$tx2gene_file)
writeLines(exclude_tx_ids, opts$output)
