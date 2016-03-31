'usage: gene_DE.R [options]

options:
 -d <dir>, --data-file <dir>           Path to RData file with gene-level data [default: gene_data.RData]
 -i <file>, --save-intermediate <file> Path to RData file for intermediate data
 -m <file>, --metadata-file <file>     Metadata file path [default: libraries.txt]
 -o <file>, --output <file>            RData output file path [default: sleuth_results.RData]
 -t <N>, --threads <N>                 Number of threads to use [default: 1]' -> doc
 
opts <- docopt::docopt(doc)
threads <- as.integer(opts$threads)

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds <- DESeq(dds, parallel=threads > 1, BPPARAM=MulticoreParam(workers=threads))
res <- results(dds)
