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
# - EPIraction placenta trophoblast enhancer-gene targets (optional, top 5000)
# The annotation depends only on the CpG -> done once and merged with each model.
# Summary tables of top CpGs are made in a separate script.
#
# [INPUT]:  results/3_Clean_Results/main_t{2,3}/{24d,3pba,tcpy,dmp,dep,dap}_t{2,3}.rds
#           db/annotation_data/FP_15_E091.RData, PMD_placenta.RData, IR_Placenta.RData
#           (from https://github.com/isglobal-brge/easier/tree/main/data)
#           db/annotation_data/hg38_EPIC_annotation.bed (EPIraction only)
# [OUTPUT]: results/4_Annotation/cpg_annotation_all.rds
#           results/4_Annotation/main_t{2,3}/<exposure>_t{2,3}_annotated.rds
#           (all CpGs, ordered by p-value; + EPIraction_placenta if run)

#==============================================================================#----
# SET UP                                                                       #----
#==============================================================================#----
library(GenomicRanges)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
source("script/00_functions.R")

input_folder  <- "results/3_Clean_Results/"
output_folder <- "results/4_Annotation/"

exposures <- c("24d", "3pba", "tcpy", "dmp", "dep", "dap")

run_epiraction <- FALSE   # TRUE to add EPIraction (needs internet: run locally)

dir.create(paste0(output_folder, "main_t2"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(output_folder, "main_t3"), recursive = TRUE, showWarnings = FALSE)

# Check: all input files exist
files_in <- c(paste0(input_folder, "main_t2/", exposures, "_t2.rds"),
              paste0(input_folder, "main_t3/", exposures, "_t3.rds"))
files_in[!file.exists(files_in)]   # should be character(0)

#==============================================================================#----
# LOAD ANNOTATION DATABASES                                                    #----
#==============================================================================#----

# 1. Illumina hg19 annotation (EPIC)
annotation <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
annotation <- annotation[, c("chr", "pos", "Relation_to_Island",
                             "UCSC_RefGene_Name", "UCSC_RefGene_Group")]

# 2. Chromatin states - fetal placenta, 15 states (GRanges)
load("db/annotation_data/FP_15_E091.RData")    # FP_15_E091
FP_15_E091
names(mcols(FP_15_E091))                        # check: state label column = "name"

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

all_cpgs <- unique(unlist(lapply(files_in, function(f) readRDS(f)$CpG)))
length(all_cpgs)

cpg_annot <- annotate_cpgs(data.frame(CpG = all_cpgs), annotation,
                           FP_15_E091, gr_pmd, gr_gdmr)

# Check: CpGs without hg19 position, and % with each annotation
sum(is.na(cpg_annot$chr))
round(100 * colMeans(!is.na(cpg_annot[, c("States15_FP", "PMD", "gDMR")])), 1)

saveRDS(cpg_annot, paste0(output_folder, "cpg_annotation_all.rds"))

#==============================================================================#----
# MERGE ANNOTATION WITH EACH MODEL                                             #----
#==============================================================================#----

for (tp in c("t2", "t3")) {
  for (i in seq_along(exposures)) {
    
    file_in <- paste0(input_folder, "main_", tp, "/", exposures[i], "_", tp, ".rds")
    cat("Exposure:", exposures[i], "| Time:", tp, "\n")
    ewas_results <- readRDS(file_in)
    
    ewas_annot <- merge(ewas_results, cpg_annot, by = "CpG", all.x = TRUE)
    ewas_annot <- ewas_annot[order(ewas_annot$pvalue), ]
    stopifnot(nrow(ewas_annot) == nrow(ewas_results))   # merge kept all CpGs once
    
    saveRDS(ewas_annot, paste0(output_folder, "main_", tp, "/", exposures[i], "_", tp,
                               "_annotated.rds"))
  }
}

#==============================================================================#----
# EPIraction - TOP 5000 CpGs PER MODEL (optional)                              #----
#==============================================================================#----
# Placenta trophoblast enhancer-gene interactions from EPIraction
# (https://epiraction.crg.eu/). Input: hg38 positions of the top 5000 CpGs.
# Needs internet access and the hg38 EPIC annotation made once with UCSC
# liftOver (steps in the diet project annotation script).
# 
# if (run_epiraction) {
#   library(httr)
#   library(rvest)
#   library(R.utils)
#   library(data.table)
#   
#   hg38_annotation <- read.table("db/annotation_data/hg38_EPIC_annotation.bed",
#                                 sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#   colnames(hg38_annotation) <- c("chr_hg38", "start_hg38", "end_hg38", "CpG")
#   
#   for (tp in c("t2", "t3")) {
#     for (i in seq_along(exposures)) {
#       
#       file_base  <- paste0(output_folder, "main_", tp, "/", exposures[i], "_", tp)
#       ewas_annot <- readRDS(paste0(file_base, "_annotated.rds"))   # ordered by p
#       
#       # Top 5000 CpGs with hg38 position (EPIraction format: chr pos CpG)
#       top5000 <- head(ewas_annot[, c("CpG", "pvalue")], 5000)
#       top5000 <- merge(top5000, hg38_annotation, by = "CpG")
#       top5000 <- top5000[order(top5000$pvalue), c("chr_hg38", "start_hg38", "CpG")]
#       
#       file_txt <- paste0(file_base, "_top5000_hg38.txt")
#       file_gz  <- paste0(file_txt, ".gz")
#       write.table(top5000, file_txt, sep = "\t", col.names = FALSE,
#                   row.names = FALSE, quote = FALSE)
#       R.utils::gzip(file_txt, destname = file_gz, remove = FALSE, overwrite = TRUE)
#       
#       # Query EPIraction (all tissues), keep placenta trophoblast
#       epi_file <- get_annot_EPIraction(gzip_file = file_gz,
#                                        out_folder = paste0(output_folder, "main_", tp))
#       if (is.null(epi_file)) next
#       epi <- fread(epi_file, data.table = FALSE)
#       epi_pla <- epi[epi$tissue == "Placenta.Trophoblast", ]
#       
#       # One row per CpG: target genes collapsed with ";"
#       epi_by_cpg <- aggregate(symbol ~ annotation, data = epi_pla,
#                               FUN = function(x) paste(unique(x), collapse = ";"))
#       colnames(epi_by_cpg) <- c("CpG", "EPIraction_placenta")
#       
#       ewas_annot$EPIraction_placenta <- epi_by_cpg$EPIraction_placenta[match(ewas_annot$CpG, epi_by_cpg$CpG)]
#       saveRDS(ewas_annot, paste0(file_base, "_annotated.rds"))
#       cat(exposures[i], tp, "- CpGs with EPIraction target:",
#           sum(!is.na(ewas_annot$EPIraction_placenta)), "\n")
#     }
#   }
# }

#==============================================================================#----
# END OF SCRIPT                                                                #----
#==============================================================================#----