tximport <- function(ids, abundance_files, tx2gene=NULL, ignoreTxVersion=FALSE) {
    matrices <- do.call(abind::abind, lapply(abundance_files, function(f) {
        a <- read.table(f, sep="\t", header=T, stringsAsFactors=FALSE)
        abundances <- a[,c("tpm", "est_counts", "eff_length")]
        rownames(abundances) <- a$target_id
        abundances
    }, along=3)
    dimnames(matrices)[3] <- ids
    
    if (is.null(tx2gene)) {
        tx2gene <- data.frame(
            tx_id=stringr::str_extract(rownames(abundance_matrix), 'ENSTR?[\\d\\.]+'),
            gene_id=stringr::str_extract(rownames(abundance_matrix), 'ENSGR?[\\d\\.]+')
        )
    }
    else {
        if (typeof(tx2gene) == "character") {
            tx2gene <- read.table(tx2gene, sep="\t", header=TRUE, stringsAsFactors=FALSE)
        }
        colnames(tx2gene) <- c("tx_id", "gene_id")
        m <- match(rownames(abundance_matrix), tx2gene$tx_id)
        if (any(is.na(m))) {
            stop("One or more transcript IDs missing from tx2gene mapping")
        }
        tx2gene <- tx2gene[m,]
    }

    txi <- list(abundance=matrices[,1,], counts=matrices[,2,], length=matrices[,3,],
                countsFromAbundance=NULL)

    tximportData::summarizeToGene(txi, tx2gene, ignoreTxVersion, "no")
}

'usage: aggregate_tx.R [options]

options:
 -d <dir>, --kallisto-dir <dir>     Diretory with Kallisto results [default: .]
 -m <file>, --metadata-file <file>  Metadata file path [default: libraries.txt]
 -o <file>, --output <file>         RData output file path [default: gene_data.RData]
 -x <file>, --tx2gene-file <file>   TSV file with mappings from transcript to gene IDs' -> doc

opts <- docopt::docopt(doc)
metadata <- read.table(opts$metadata_file, sep="\t", header=T, stringsAsFactors=F)
abundance_files <- file.path(opts$kallisto_dir, metadata$UniqueID, "abundance.tsv")
gene_data <- tximport(metadata$UniqueID, abundance_files)
save(gene_data, file=opts$output)
