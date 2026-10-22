#==============================================================================#
# Title: Prenatal Exposure to Pesticides and placental DNA methylation in the 
# BiSC and SEPAGES cohorts. 
# Author: Joana Llauradó
# Contact: joana.llaurado@isglobal.org                              
# Date: Stay in Grenoble October - November 2026
#==============================================================================#
# [DESCRIPTION]: Functional annotation of all CpGs of each EWAS model:
# - hg19 annotation: position, relation to CpG island, gene name and group
# - 15 chromatin states in fetal placenta (FP_15_E091)
# - Partially methylated domains in placenta (PMD_placenta)
# - Germline differentially methylated regions / imprinted regions (IR_Placenta)

# [INPUT]:
#           db/annotation_data/FP_15_E091.RData, PMD_placenta.RData, IR_Placenta.RData
# [OUTPUT]: results/5_Annotation/cpg_annotation_all.rds


#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(GenomicRanges)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

output_folder <- "results/5_Annotation/"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
#==============================================================================#----
# ANNOTATION                                                                   #----
#==============================================================================#----

# Overlap CpGs with a database of genomic regions; one row per CpG, multiple
# overlapping regions collapsed with ";"
# - gr_cpg: GRanges of CpGs with mcols "name" = CpG ID
# - gr_db:  GRanges of the database; value_col = column with the region label
overlap_collapse <- function(gr_cpg, gr_db, value_col = "name") {
  hits <- findOverlaps(gr_cpg, gr_db)
  if (length(hits) == 0) return(data.frame(CpG = character(0), value = character(0)))
  df_hits <- data.frame(CpG   = gr_cpg$name[queryHits(hits)],
                        value = as.character(mcols(gr_db)[[value_col]][subjectHits(hits)]))
  aggregate(value ~ CpG, data = df_hits, FUN = function(x) paste(unique(x), collapse = ";"))
}

# Unique elements of a ";"-separated string (e.g. "GENE1;GENE1;GENE2" -> "GENE1;GENE2")
unique_split <- function(x) {
  sapply(strsplit(as.character(x), ";"), function(v) paste(unique(v), collapse = ";"))
}


#==============================================================================#----
# LOAD ANNOTATION DATABASES                                                    #----
#==============================================================================#----

# 1. Illumina hg19 annotation (EPIC)
annotation <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
annotation <- annotation[, c("chr", "pos", "Relation_to_Island",
                             "UCSC_RefGene_Name", "UCSC_RefGene_Group")]
annotation$CpG <- rownames(annotation)
annotation$gene_name <- unique_split(annotation$UCSC_RefGene_Name)
annotation$gene_group <- unique_split(annotation$UCSC_RefGene_Group)
annotation <- annotation[, c("CpG", "chr", "pos", "Relation_to_Island", "gene_name", "gene_group")]
dim(annotation)   # ~866k CpGs


# 2. Chromatin states - fetal placenta, 15 states (GRanges)
load("db/annotation_data/FP_15_E091.RData")    # FP_15_E091
FP_15_E091
names(mcols(FP_15_E091))                        

# 3. PMDs placenta (data.frame -> GRanges)
load("db/annotation_data/PMD_placenta.RData")  # PMD_placenta
gr_pmd <- GRanges(seqnames = PMD_placenta$Chr_PMD,
                  ranges   = IRanges(PMD_placenta$Start_PMD, end = PMD_placenta$End_PMD),
                  name     = paste(PMD_placenta$Chr_PMD, PMD_placenta$Start_PMD,
                                   PMD_placenta$End_PMD, sep = ":"))

# 4. gDMRs / imprinted regions placenta (data.frame -> GRanges)
load("db/annotation_data/IR_Placenta.RData")   # IR_Placenta
gr_gdmr <- GRanges(seqnames = IR_Placenta$Chr_DMR,
                   ranges   = IRanges(IR_Placenta$Start_DMR, end = IR_Placenta$End_DMR),
                   name     = paste(IR_Placenta$Chr_DMR, IR_Placenta$Start_DMR,
                                    IR_Placenta$End_DMR, sep = ":"))

# Check: same chromosome naming in all databases (e.g. "chr1")
head(unique(annotation$chr)); head(seqlevels(FP_15_E091))
head(seqlevels(gr_pmd));      head(seqlevels(gr_gdmr))

#==============================================================================#----
# ANNOTATE ALL CpGs (ONCE)                                                     #----
#==============================================================================#----
# CpGs as genomic ranges
gr_cpg <- GRanges(seqnames = annotation$chr,
                  ranges   = IRanges(start = annotation$pos, end = annotation$pos),
                  name     = annotation$CpG)

# Overlaps with each database (NA = no overlap)
states <- overlap_collapse(gr_cpg, FP_15_E091)
pmd    <- overlap_collapse(gr_cpg, gr_pmd)
gdmr   <- overlap_collapse(gr_cpg, gr_gdmr)

annotation$States15_FP <- states$value[match(annotation$CpG, states$CpG)]
annotation$PMD         <- pmd$value[match(annotation$CpG, pmd$CpG)]
annotation$gDMR        <- gdmr$value[match(annotation$CpG, gdmr$CpG)]

# Check: no duplicated CpGs, % with each annotation
sum(duplicated(annotation$CpG))                                                  # 0
round(100 * colMeans(!is.na(annotation[, c("States15_FP", "PMD", "gDMR")])), 1)
head(annotation)

saveRDS(annotation, paste0(output_folder, "cpg_annotation_all.rds"))


#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----