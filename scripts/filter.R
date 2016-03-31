library(stringr)
library(data.table)
library(docopt)

# Return a list of length C, where C is the number of 
# experimental conditions. Each list item is a vector
# with names as sample names and values as abundance.tsv
# file paths.
get_quant_files <- function(condition, metadata.file, quant.dir, config) {
	md <- read.table(metadata.file, sep="\t", header=T, stringsAsFactors=F)
	if (!(condition %in% colnames(md))) {
		stop(paste(condition, "is not a column in", metadata.file))
	}
	if ("FileNames" %in% colnames(md)) {
		tapply(md$FileNames, md[,condition], identity, simplify=FALSE)
	}
	else if ("LibraryType" %in% colnames(md)) {

	}
	else {
		if ("NumFiles" %in% colnames(md)) {

		}
		else {

		}
	}
}

# Load and sum all quant files in quant.files vectors.
# Returns a single vector of TPM values with names as
# transcript names.
sum_abundances <- function(quant.files) {
	abundances <- NULL
	for (i in 1:length(quant.files)) {
		sample <- names(quant.files)[i]
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

# TODO: this is specific to Ensembl/GENCODE annotations.
filter_transcripts <- function(quant_files, cutoff=0.1) {
	abundances <- data.table(do.call(cbind, lapply(quant_files, sum_abundances)))
	abundances$tx_id <- str_extract(rownames(abundances), 'ENSTR?[\\d\\.]+')
	abundances$gene_id <- str_extract(rownames(abundances), 'ENSGR?[\\d\\.]+')
	abundances[, rel_abundance := tpm / sum(tpm), gene_id]
	abundances[rel_abundance < cutoff, tx_id]
}

'usage: filter_low_abundance_tx.R [options]

options:
 -c <name>, --condition <name>      Name of condition used to partition samples [default: Condition]
 -d <dir>, --kallisto-dir <dir>     Diretory with Kallisto results [default: .]
 -f <file>, --config-file <file>    Config file path [default: ]
 -m <file>, --metadata-file <file>  Metadata file path [default: metadata.txt]
 -o <file>, --output <file>         Output file path [default: exclude_tx_ids.txt]
 -t <float>, --cutoff <float>       Cutoff value between 0-1 [default: 0.1]' -> doc

opts <- docopt(doc)
quant_files <- get_quant_files(opts$condition, opts$metadata_file, opts$kallisto_dir, opts$config_file)
exclude_tx_ids <- filter_transcripts(quant_files, as.numeric(opts$cutoff))
writeLines(exclude_tx_ids, opts$output)
