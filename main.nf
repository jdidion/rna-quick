#!/usr/bin/env nextflow

/* Try to create a URL; return false if url is not a value URL. */
def parseUrl(url) {
    try {
        new URL(url)
    } 
    catch (MalformedURLException e) {
        false
    }
}

/* Get a URL for a genome shortcut, or throw an exception if the shortcut is invalid. */
def parseGenomeShortcut(shortcut) {
    def parts = shortcut.split('.')
    if (parts[0] == "ensembl") {
        release = int(parts[1])
        if (parts[2][0:3] == "GRCh") {
            organism = "homo_sapiens"
            build = parts[2]
        }
        else {
            organism = parts[2]
            build = parts[3]
        }
        URL('ftp://ftp.ensembl.org/pub/release-${release}/fasta/${organism}/dna/${organism}.${build}.dna.toplevel.fa.gz')
    }
    else if (parts[0] == "ucsc") {
        build = parts[1]
        URL('http://hgdownload.cse.ucsc.edu/goldenPath/${build}/bigZips/chromFa.tar.gz')
    }
    else {
        throw Exception('Invalid genome shortcut: $shortcut')
    }
}

/* Get a URL for an annotation shortcut, or throw an exception if the shortcut is invalid. */
def parseAnnotationShortcut(shortcut) {
    def parts = shortcut.split('.')
    if (parts[0] == "gencode") {
        release = parts[1]
        content = "annotation" if parts[2] == "comprehensive" else "long_noncoding_RNAs"
        regions = "" if parts[3] == "all" else "chr_patch_hapl_scaff"
        URL('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${release}/gencode.v${release}.${regions}${content}.gtf.gz')
    }
    else {
        throw Exception('Invalid genome shortcut: $shortcut')
    }
}

/* Download the contents of a URL to a file. */
def download(url, outfile) {
    java.nio.ReadableByteChannel input = java.nio.Channels.newChannel(url.openStream())
    java.nio.FileOutputStream output = new java.nio.FileOutputStream(outfile)
    output.getChannel().transferFrom(input, 0, Long.MAX_VALUE)
}

/* Download genome and annotation files, if necessary. */
process init {
    exec:

    def genome = params.genome
    if (!file(genome).exists) {
        url = parseUrl(genome)
        if (!url) {
            url = parseGenomeShortcut(genome)
        }
        outfile = path.join(
            params.reference_dir,
            url.substring(url.lastIndexOf('/')+1, url.length())
        )
        if (!file(outfile).exists) {
            download(url, outfile)
        }
        params.genome = outfile
    }
    
    def annotation = params.annotation
    if (!file(annotation).exists) {
        url = parseUrl(annotation)
        if (!url) {
            url = parseAnnotationShortcut(annotation)
        }
        outfile = path.join(
            params.reference_dir,
            url.substring(url.lastIndexOf('/')+1, url.length())
        )
        if (!file(outfile).exists) {
            download(url, outfile)
        }
        params.annotation = outfile
    }
}

