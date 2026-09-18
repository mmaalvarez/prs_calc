#!/usr/bin/env Rscript

library(tidyverse)
library(vcfR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")


if(interactive()){setwd("/g/strcombio/fsupek_data/users/malvarez/projects/lucia/prs/external/hartwig/calculate_PRS/bin")}

args = commandArgs(trailingOnly=TRUE)

## which sample to calculate PRS in this channel
sampleID = ifelse(interactive(),
                  yes = "ACTN01020190R",
                  no = args[1])

## info on which samples are lung tumor and which are non-lung tumors, and paths to VCFs
metadata_table = ifelse(interactive(),
                        yes = "/g/strcombio/fsupek_data/users/malvarez/projects/lucia/data/vcf_bam/external/hartwig/hartwig__lung_vs_control_tissues.tsv",
                        no = args[2]) %>% 
  read_tsv

## get models, to keep only those variants from the VCF
prs_model = ifelse(interactive(),
              yes = "/g/strcombio/fsupek_data/users/malvarez/projects/lucia/data/prs_models/published_PRS_models/young_etal_2025/MDPI_SuppTableS3_hmPOS_GRCh37.txt",
              no = args[3]) %>% 
  read_tsv(., comment = "#") %>% 
  (\(x) {
    if ("other_allele" %in% names(x) && !"hm_inferOtherAllele" %in% names(x))
      rename(x, any_of(c(other_allele = "hm_inferOtherAllele")))
    else
      x
  })() %>%
  select(effect_allele, effect_weight, chr_name, chr_position) %>% 
  unite("coord", chr_name, chr_position, sep = "_")


## get VCF info
# read VCF
vcf = read.vcfR(Sys.glob(metadata_table %>% filter(IID == sampleID) %>% pull(path)))
# keep only snps
snp_mask = nchar(vcf@fix[,"REF"]) == 1  &  nchar(vcf@fix[,"ALT"]) == 1
vcf_snps = vcf[snp_mask]
rm(vcf)
gc()
# get which are the REF (0) and ALT (1) alleles
ref_alt_alleles = vcf_snps@fix %>% 
  as_tibble %>% 
  select(CHROM, POS, REF, ALT) %>% 
  unite("coord", CHROM, POS, sep = "_")

## run predictions
predictions = vcf_snps %>% 
  # get genotypes
  extract.gt %>% 
  data.frame %>% 
  rownames_to_column("coord") %>% 
  as_tibble() %>%
  pivot_longer(cols = !matches("coord"), names_to = "sampleId", values_to = "genotype") %>%
  separate(genotype, into = c("A1", "A2"), sep = "/|\\|") %>% 
  # keep only the SNPs in prs model
  right_join(prs_model) %>% 
  # get which are the REF (0) and ALT (1) alleles
  left_join(ref_alt_alleles) %>% 
  # WARNING: assuming REF/REF (i.e. 0/0) for missing variants...
  rowwise %>% 
  mutate(A1 = ifelse(is.na(A1), "0", A1),
         A2 = ifelse(is.na(A2), "0", A2),
         # fill sampleId as well, it's missing when the variant was not in the VCF
         IID = sampleID) %>% 
  ungroup


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
ref_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg19, gr)

# 4. Convert the output (DNAStringSet) to character and add to your dataframe
missing_REF_info$REF = as.character(ref_seqs)
missing_REF_info = missing_REF_info %>% 
  unite("coord", chrom, pos, sep = "_") %>% 
  mutate(coord = gsub("chr", "", coord))


### 
## continue with predictions...
predictions = predictions %>% 
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
  group_by(IID) %>% 
  summarise(prs = sum(effect_by_dosage)) %>% 
  left_join(metadata_table) %>% 
  mutate(y = case_when(type == "Control" ~ 0,
                       type == "Lung_tumor" ~ 1)) %>% 
  select(IID, prs, type, y)

write_tsv(predictions,
          "predictions.tsv")
