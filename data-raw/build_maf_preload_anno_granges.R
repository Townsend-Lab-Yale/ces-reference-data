library(data.table)
library(GenomicRanges)

## Step 1: Population variant frequency data source and generation
# Annovar has publicly available gnomAD data in an easy-to-use table format (large files)
# Need these columns from exome file: #Chr, Start, End, non_cancer_AF_popmax (columns 1, 2, 3, 21)
# awk -F "\t" 'BEGIN {OFS = "\t"} $21 > .01 {print $1, $2, $3, $21}' hg19_gnomad211_exome.txt > common_germline_from_annovar_gnomad211_exome.txt
# 
# For genome file, there is no "non-cancer" cohort, so use AF_popmax (field 7) instead
# awk -F "\t" 'BEGIN {OFS = "\t"} $21 > .01 {print $1, $2, $3, $21}' hg19_gnomad211_exome.txt > common_germline_from_annovar_gnomad211_exome.txt

# Read in data files
common_exome = fread("common_germline_from_annovar_gnomad211_exome.txt")

# verify all numeric, and already confirmed that they are MAF-like intervals (1-based, closed)
common_exome[, all(is.numeric(non_cancer_AF_popmax) & ! is.na(non_cancer_AF_popmax))]
common_exome = common_exome[, .(chr = `#Chr`, start = Start, end = End)]
exome_gr = makeGRangesFromDataFrame(common_exome, starts.in.df.are.0based = F)


common_genome = fread("common_germline_from_annovar_gnomad211_genome.txt")
common_genome[, all(is.numeric(AF_popmax) & ! is.na(AF_popmax))]
common_genome = common_genome[, .(chr = `#Chr`, start = Start, end = End)]
genome_gr = makeGRangesFromDataFrame(common_genome, starts.in.df.are.0based = F)
combined = sort(reduce(c(exome_gr, genome_gr)))

# .xz compression seems to be most efficient on GRanges (worth it for large GRanges)
setattr(combined, "anno_col_name", "germline_variant_site")
saveRDS(combined, 'gnomad_common_variation_granges.rds', compress = 'xz')



## Step 2: RepeatMasker
# Download RepeatMasker hg19 annotations from UCSC genome browser
rmsk = fread("~/Documents/reference/rmsk.hg19.bed") # will respect BED coordinate format
rmsk = rmsk[, .(genoName, genoStart, genoEnd)]
rmsk = makeGRangesFromDataFrame(rmsk, seqnames.field = "genoName", start.field = "genoStart", end.field = "genoEnd",
                                keep.extra.columns = F, ignore.strand = T, starts.in.df.are.0based = T)
seqlevelsStyle(rmsk) = "NCBI" # remove chr prefixes to match curr_maf data
rmsk = rmsk[seqnames(rmsk) %in% c(1:22, 'X', 'Y')] # restrict to primary contigs
rmsk = reduce(rmsk)
setattr(rmsk, "anno_col_name", "repetitive_region")
saveRDS(rmsk, "repeat_masked_hg19_granges.rds", compress = 'xz')


## Step 3: COSMIC mutations
# Using COSMIC census (coding) mutations from release 92
# Note file also has a GRCh38 column, if you want to replicate this for hg38
cmc = fread("/Users/Jeff/Downloads/cmc_export.v92.tsv/cmc_export.tsv", 
            select = c("Mutation genome position GRCh37", "MUTATION_SIGNIFICANCE_TIER"))


cmc[, c("chr", "start", "end") := tstrsplit(`Mutation genome position GRCh37`, split='[-:]')]
cmc = cmc[chr %in% 1:24 & ! is.na(start) & ! is.na(end), .(chr, start, end, cosmic_mut_tier = MUTATION_SIGNIFICANCE_TIER)]
cmc[chr == "23", chr := "X"]
cmc[chr = "24", chr := "Y"]

muts_by_tier = makeGRangesFromDataFrame(cmc, keep.extra.columns = T)
muts_by_tier = as.list(split(cosmic_muts, cosmic_muts$cosmic_mut_tier))
muts_by_tier[[1]] = reduce(muts_by_tier[[1]])
muts_by_tier[[2]] = reduce(setdiff(muts_by_tier[[2]], muts_by_tier[[1]]))
muts_by_tier[[3]] = reduce(setdiff(muts_by_tier[[3]], c(muts_by_tier[[2]], muts_by_tier[[1]])))
muts_by_tier[[4]] = reduce(setdiff(muts_by_tier[[4]], c(muts_by_tier[[3]], muts_by_tier[[2]], muts_by_tier[[1]])))


setattr(muts_by_tier, "anno_col_name", "cosmic_site_tier")
saveRDS(muts_by_tier, "cmc_92_gr_by_tier.rds", compress = 'xz')



