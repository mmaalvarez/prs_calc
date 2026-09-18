#!/usr/bin/env Rscript

library(tidyverse)
library(vcfR)
library(GenomicRanges)
build="38"
library(paste0("BSgenome.Hsapiens.UCSC.hg", build), character.only = T)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("first", "dplyr")

if(interactive()){setwd("/g/strcombio/fsupek_data/users/malvarez/projects/lucia/prs/external/TCGA_WGS/test_manual_calculation/")}

args = commandArgs(trailingOnly=TRUE)

## which vcf to calculate PRS
vcf_name = ifelse(interactive(),
                  yes = "/g/strcombio/fsupek_franklin/malvarez/test_pgsc_calc/vcfs_for_manual_prs_calc/all_chr.withheader.vcf",
                  #yes = "/g/strcombio/fsupek_franklin/malvarez/test_pgsc_calc/vcfs_for_manual_prs_calc/f792fc7b-6cbf-4b4d-ac5e-8d82322f63d9_final.vcf",
                  no = args[1])

## get models, to keep only those variants from the VCF
prs_model = ifelse(interactive(),
                   yes = paste0("/g/strcombio/fsupek_data/users/malvarez/projects/lucia/data/prs_models/published_PRS_models/hung_etal_2021/PGS000740_hmPOS_GRCh", build, ".txt"),
                   no = args[2]) %>% 
  read_tsv(., comment = "#") %>% 
  (\(x) {
    if ("other_allele" %in% names(x) && !"hm_inferOtherAllele" %in% names(x))
      rename(x, any_of(c(other_allele = "hm_inferOtherAllele")))
    else
      x
  })() %>%
  # warning: make sure chr_position are in the required build, sometimes chr_position is in hg37 and hm_position is hg38
  select(effect_allele, other_allele, effect_weight, chr_name, chr_position) %>% 
  unite("coord", chr_name, chr_position, sep = "_")


## get VCF info
# read VCF
vcf = read.vcfR(vcf_name)
gc()

# get which are the REF (0) and ALT (1) alleles
ref_alt_alleles = vcf@fix %>% 
  as_tibble %>% 
  select(CHROM, POS, REF, ALT) %>% 
  mutate(CHROM = gsub("chr", "", CHROM),
         ALT = gsub("<NON_REF>|,<NON_REF>", "", ALT),
         ALT = ifelse(ALT == "", NA, ALT)) %>% 
  unite("coord", CHROM, POS, sep = "_")
gc()

## run predictions
predictions = vcf %>% 
  # get genotypes
  extract.gt
gc()

rownames(predictions) = paste(gsub("chr","",getCHROM(vcf)), getPOS(vcf), sep = "_")

predictions = predictions %>%  
  data.frame %>% 
  rownames_to_column("coord") %>% 
  as_tibble() %>%
  pivot_longer(cols = !matches("coord"), names_to = "sampleId", values_to = "genotype") %>%
  separate(genotype, into = c("A1", "A2"), sep = "/|\\|") %>% 
  # keep only the SNPs in prs model
  right_join(prs_model)
gc()
rm(vcf)
gc()

predictions = predictions %>%  
  # get which are the REF (0) and ALT (1) alleles
  left_join(ref_alt_alleles) %>% 
  # WARNING: assuming REF/REF (i.e. 0/0) for missing variants
  rowwise %>% 
  mutate(A1 = ifelse(is.na(A1), "0", A1),
         A2 = ifelse(is.na(A2), "0", A2)) %>% 
  ungroup
gc()
rm(ref_alt_alleles)
gc()

## effect_allele in PRS model may not refer to the ALT allele, so retrieve REF/ALT alleles info
# 1. Create a sample dataframe with your coordinates
missing_REF_info = predictions %>% 
  filter(is.na(REF)) %>% 
  separate(coord, into = c("chrom", "pos"), sep = "_") %>% 
  mutate(chrom = paste0("chr", chrom),
         chrom = gsub("chrchr", "chr", chrom)) %>% 
  select(chrom, pos)

# 2. Convert dataframe to a GRanges object
# For a single nucleotide (SNP), the start and end positions are the same
gr = GRanges(seqnames = missing_REF_info$chrom,
             ranges = IRanges(start = as.numeric(missing_REF_info$pos), end = as.numeric(missing_REF_info$pos)))

# 3. Retrieve the REF sequences
ref_seqs = getSeq(get(paste0("BSgenome.Hsapiens.UCSC.hg", build)), gr)

# 4. Convert the output (DNAStringSet) to character and add to your dataframe
missing_REF_info$REF = as.character(ref_seqs)
missing_REF_info = missing_REF_info %>% 
  unite("coord", chrom, pos, sep = "_") %>% 
  mutate(coord = gsub("chr", "", coord))


### 
## continue with predictions...
predictions_final = predictions %>% 
  ## add the REF info
  rows_patch(missing_REF_info, by = "coord") %>% 
  # A1 and A2 should have a "0" for the non-effect allele, and a "1" for the effect allele
  # NOTE: this only matters for homozygous 1/1 (i.e. ALT/ALT) when the effect allele in the model refers to the REF allele --> should then become 0/0 (i.e. no effect)
  mutate(A1_harmonised = case_when(effect_allele == REF ~ 1 - as.numeric(A1),           # Effect is REF
                                   !is.na(ALT) & effect_allele == ALT ~ as.numeric(A1), # Effect is ALT
                                   TRUE ~ 0),
         A2_harmonised = case_when(effect_allele == REF ~ 1 - as.numeric(A2),
                                   !is.na(ALT) & effect_allele == ALT ~ as.numeric(A2),
                                   TRUE ~ 0),
         allele_dosage = A1_harmonised + A2_harmonised,
         effect_by_dosage = effect_weight * allele_dosage) %>% 
  # fill NAs in sampleId column with the sample ID
  mutate(sampleId = replace_na(sampleId, first(na.omit(sampleId)))) %>% 
  group_by(sampleId) %>% 
  summarise(prs = sum(effect_by_dosage)) %>% 
  select(sampleId, prs)

write_tsv(predictions_final,
          paste0(gsub("vcfs\\/|\\.vcf","",vcf_name), "_predictions.tsv"))
