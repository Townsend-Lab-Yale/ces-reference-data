library(devtools)
load_all()

ucsc_info = getFromNamespace(".UCSC_cached_chrom_info", "GenomeInfoDb")[["hg19"]]
ucsc_info = GenomeInfoDb:::.add_ensembl_column(ucsc_info, "hg19")
ncbi_info = getFromNamespace(".NCBI_cached_chrom_info", "GenomeInfoDb")[["GCF_000001405.25"]]

cached_chromInfo = list()
cached_chromInfo[['UCSC']] = list(name = "hg19", value = ucsc_info)
cached_chromInfo[['NCBI']] = list(name = "GCF_000001405.25", value = ncbi_info)
saveRDS(cached_chromInfo, 'inst/refset/cached_chromInfo.rds')