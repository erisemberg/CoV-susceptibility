library(readr)
suppressMessages(library(tidyverse))
library(stringi)
source('/Users/ellenrisemberg/Documents/ValdarFerris/scripts/qtl_functions.R')

setwd("/Users/ellenrisemberg/Documents/ValdarFerris/Coronavirus/virus-research")

pc_var_types <- c('missense_variant', 'splice_region_variant', 'stop_gained',
                  'stop_lost', 'splice_donor_variant', 'splice_acceptor_variant',
                  'inframe_deletion', 'initiator_codon_variant',  'coding_sequence_variant', 
                  'incomplete_terminal_codon_variant', 'inframe_insertion', 'frameshift_variant',
                  'splice_', 
                  'stop_lost&NMD__variant')
reg_var_types <- c('downstream_gene_variant', 'intron_variant', 'upstream_gene_variant',
                   'non_coding_transcript_exon_variant', 'non_coding_transcript_variant',
                   '3_prime_UTR_variant', '5_prime_UTR_variant', 'synonymous_variant',
                   'NMD_transcript_variant', 'mature_miRNA_variant', 'stop_retained_variant',
                   'non_coding__exon_variant')


#-------------------------------CC006xCC044------------------------------------#
# read VCFs 
#vcf644_1 <- read.csv(file = "source_data/chr9-variants/sanger.founders.CC006xCC044.chr9.89.54-96Mb.NZO_HlLtJ.NOD_ShiLtJ.genes.csv")
vcf644_2 <- read.csv(file = "source_data/chr9-variants/sanger.founders.CC006xCC044.chr9.98.92-110Mb.129S1_SvImJ.NOD_ShiLtJ.vcf.genes.csv")
vcf644_3 <- read.csv(file = "source_data/chr9-variants/sanger.founders.CC006xCC044.chr9.101-114Mb.NZO_HlLtJ.NOD_ShiLtJ.vcf.genes.csv")
vcf644_4 <- read.csv(file = "source_data/chr9-variants/sanger.founders.CC006xCC044.chr9.110-122Mb.NZO_HlLtJ.WSB_EiJ.vcf.genes.csv")
vcf644_5 <- read.csv(file = "source_data/chr9-variants/sanger.founders.CC006xCC044.chr9.118-123.56Mb.129S1_SvImJ.WSB_EiJ.vcf.genes.csv")

# change column names - col 9-10 indicate compared strains / make CC006 haplotype first 
#colnames(vcf644_1)[9:10] <- c('CC044', 'CC006') # originally NOD, NZO
#vcf644_1 <- vcf644_1[,c(names(vcf644_1)[1:8], 'CC006', 'CC044')] # switch order so CC006 is listed first 
colnames(vcf644_2)[9:10] <- c('CC006', 'CC044') # originally 129S1, NOD
colnames(vcf644_3)[9:10] <- c('CC044', 'CC006') # originally NOD, NZO 
vcf644_3 <- vcf644_3[,c(names(vcf644_3)[1:8], 'CC006', 'CC044')] # switch order so CC006 is listed first
colnames(vcf644_4)[9:10] <- c('CC006', 'CC044') # originally NZO, WSB
colnames(vcf644_5)[9:10] <- c('CC006', 'CC044') # originally 129S1, WSB

vcfs644 <- list(vcf644_2, vcf644_3, vcf644_4, vcf644_5) # list of all dfs (removed vcf644_1)
allvcf644 <- do.call("rbind", vcfs644) # combine into one dataframe
# remove redundant 

allvcf644 <- allvcf644[,c('gene_name', names(allvcf644)[1:7], names(allvcf644)[9:10])] # re-order so gene_name is first
allvcf644 <- allvcf644[order(allvcf644$gene_name),] # order by gene name 
allvcf644 <- allvcf644[which(allvcf644$ANN != "intergenic_variant"),] # filter out intergenic
genes644 <- unique(allvcf644$gene_name) # candidate gene list 


#-------------------------------CC011xCC074------------------------------------#
vcf1174_1 <- read.csv(file = "source_data/chr9-variants/sanger.founders.CC011xCC074.chr9.74.94-96Mb.A_J.vcf.genes.csv")
vcf1174_2 <- read.csv(file = "source_data/chr9-variants/sanger.founders.CC011xCC074.chr9.92-107Mb.PWK_PhJ.A_J.vcf.genes.csv")
vcf1174_3 <- read.csv(file = "source_data/chr9-variants/sanger.founders.CC011xCC074.chr9.103-124Mb.A_J.vcf.genes.csv")
vcf1174_4 <- read.csv(file = "source_data/chr9-variants/sanger.founders.CC011xCC074.chr9.120-124.07Mb.PWK_PhJ.A_J.vcf.genes.csv")

# change column names - col 9-10 indicate compared strains / make CC074 haplotype (always A/J)  first 
colnames(vcf1174_1)[9] <- c('CC074') # originally A/J
vcf1174_1$CC011 = rep(0, nrow(vcf1174_1)) # add B6/CC011 column
colnames(vcf1174_2)[9:10] <- c('CC074', 'CC011') # originally A/J, PWK
colnames(vcf1174_3)[9] <- c('CC074') # originally A/J
vcf1174_3$CC011 = rep(0, nrow(vcf1174_3)) # add B6/CC011 column
colnames(vcf1174_4)[9:10] <- c('CC074', 'CC011') # originally A/J, PWK

vcfs1174 <- list(vcf1174_1, vcf1174_2, vcf1174_3, vcf1174_4) # list of all dfs
allvcf1174 <- do.call("rbind", vcfs1174) # combine into one dataframe

allvcf1174 <- allvcf1174[,c('gene_name', names(allvcf1174)[1:7], names(allvcf1174)[9:10])] # re-order so gene_name is first
allvcf1174 <- allvcf1174[order(allvcf1174$gene_name),] # order by gene name 
allvcf1174 <- allvcf1174[which(allvcf1174$ANN != "intergenic_variant"),] # filter out intergenic
genes1174 <- unique(allvcf1174$gene_name) # candidate gene list 


#-------------------------------intersection-----------------------------------#
# intersect 
intrsct <- intersect(genes644, genes1174)

#------------------------intersect variant analysis----------------------------#
# variant analysis for genes with variants segregating in both crosses 
allvcf644$cross <- 'CC006xCC044'
colnames(allvcf644)[9:10] <- c('susceptible', 'resistant')
allvcf1174$cross <- 'CC011xCC074'
colnames(allvcf1174)[9:10] <- c('susceptible', 'resistant')

allvcf <- rbind(allvcf644, allvcf1174)

# only want to look at genes that are in the intersection of candidate genes (filtering from there)
# allvcf <- completevcf[which(completevcf$gene_name %in% intrsct),]
# sort by gene name then variant ID 
# allvcf <- allvcf[order(allvcf$gene_name, allvcf$ID),] 

genelist <- intrsct
dqdgenes <- vector()

cand_gene_sum <- data.frame(gene_name = character(),
                            category = character(),
                            pc_vars = logical())

ix = 0
for (i in 1:length(intrsct)){ # for each gene 
  gene <- intrsct[i]
  if(is.na(gene)) { next }
  
  df <- allvcf[which(allvcf$gene_name == gene),]
  
  # remove genes from the list if there are only variants segregating in both crosses that are in opposite phase  
  if (sum(!duplicated(df$susceptible)) == 1) { next } # all variants are in the same phase 
  
  dupIDs <- df$ID[duplicated(df$ID)] # duplicated IDs 
  dupdf <- df[which(df$ID %in% dupIDs),]
  
  disqualifying_vars <- 0
  vars_processed <- 0
  category <- NA 
  # go through this dataframe of duplicated IDs, as soon as you find one where the phase is not a disqualifer, move on
  # (we would need all variants to be a disqualifier to disqualify the gene)
  for (j in 1:length(dupIDs)){
    vars_processed <- vars_processed + 1
    cur_dupID <- dupIDs[j]
    dupIDdf <- dupdf[which(dupdf$ID == cur_dupID),]
    # only care about IDs present across crosses (could be duplicated within a cross bc of haplotype overlap)
    if (sum(!duplicated(dupIDdf$cross)) == 1) { break } # all variants are from the same cross
    # variants are from both crosses. If variants are in the same phase, still a candidate gene 
    if (sum(!duplicated(dupIDdf$susceptible)) != 1){ # variants are not all in the same phase 
      category <- "both_shared" 
    } else { break } # all variants are in the same phase
    # all variants that made it to this point are potentially disqualifiers
    disqualifying_vars <- disqualifying_vars + 1 
  }
  
  # If all vars looked at are disqualifiers, gene is disqualified, remove gene from intrsct 
  if (disqualifying_vars/vars_processed == 1){
    dqdgenes <- append(dqdgenes, gene)
    genelist <- genelist[which(genelist != gene)]
    next 
  }
  
  # gene is not disqualified and does not have shared variants 
  if (is.na(category)) {  
    category <- "both_distinct"
  }
  
  ix = ix+1
  cand_gene_sum[ix,'gene_name'] <- gene 
  cand_gene_sum[ix,'category'] <- category 
}

cand_genes <- data.frame(gene = genelist,
                         num_pc_variants = rep(NA, length(genelist)),
                         num_reg_variants = rep(NA, length(genelist)),
                         vartypes = rep(NA, length(genelist)))

cand_gene_sum$pc_vars <- FALSE

# go through list again and document variant types associated with each gene 
for (k in 1:length(genelist)){
  gene = genelist[k]
  df <- allvcf[which(allvcf$gene_name == gene),]
  vartypes <- unique(df$ANN)
  
  cand_genes[k, 'num_pc_variants'] <- length(which(vartypes %in% pc_var_types))
  cand_genes[k, 'num_reg_variants'] <- length(which(vartypes %in% reg_var_types))
  cand_genes[k, 'vartypes'] <- paste(vartypes, collapse = ",")
  
  if (length(which(vartypes %in% pc_var_types)) > 0){
    cand_gene_sum[which(cand_gene_sum$gene_name == gene), 'pc_vars'] <- TRUE 
  } 
}

cand_genes$from_which_cross <- "both"


#-------------------------------CC006xCC044------------------------------------#
# genes from CC006xCC044 cross that are not shared between crosses  
# (this includes genes that were in the intersection but were DQ'd)
genes_only644 <- genes644[which(genes644 %notin% genelist)]

# count number of pc and reg variants and list candidate genes 
cand_genes644 <- data.frame(gene = genes_only644,
                            num_pc_variants = rep(NA, length(genes_only644)),
                            num_reg_variants = rep(NA, length(genes_only644)),
                            vartypes = rep(NA, length(genes_only644)))

# candidate gene summary 
cand_gene_sum644 <- data.frame(gene_name = character(),
                               category = character(),
                               pc_vars = logical())
ix = 0

# go through list again and document variant types associated with each gene 
for (k in 1:length(genes_only644)){
  gene = genes_only644[k]
  df <- allvcf[which(allvcf$gene_name == gene),]
  vartypes <- unique(df$ANN)
  
  cand_genes644[k, 'num_pc_variants'] <- length(which(vartypes %in% pc_var_types))
  cand_genes644[k, 'num_reg_variants'] <- length(which(vartypes %in% reg_var_types))
  cand_genes644[k, 'vartypes'] <- paste(vartypes, collapse = ",")
  
  ix = ix+1
  cand_gene_sum644[ix, 'gene_name'] <- gene
  if (length(which(vartypes %in% pc_var_types)) > 0){ 
    cand_gene_sum644[ix, 'pc_vars'] <- TRUE 
  } else {
    cand_gene_sum644[ix, 'pc_vars'] <- FALSE
  }
}

cand_gene_sum644$category <- "CC006xCC044"
cand_genes644$from_which_cross <- "CC006xCC044"

#-------------------------------CC011xCC074------------------------------------#
# genes from CC011xCC074 cross that are not shared between crosses  
# (this includes genes that were in the intersection but were DQ'd)
genes_only1174 <- genes1174[which(genes1174 %notin% genelist)]

# count number of pc and reg variants and list candidate genes
cand_genes1174 <- data.frame(gene = genes_only1174,
                             num_pc_variants = rep(NA, length(genes_only1174)),
                             num_reg_variants = rep(NA, length(genes_only1174)),
                             vartypes = rep(NA, length(genes_only1174)))

# candidate gene summary 
cand_gene_sum1174 <- data.frame(gene_name = character(),
                                category = character(),
                                pc_vars = logical())
ix = 0

# go through list again and document variant types associated with each gene 
for (k in 1:length(genes_only1174)){
  gene = genes_only1174[k]
  df <- allvcf[which(allvcf$gene_name == gene),]
  vartypes <- unique(df$ANN)
  
  cand_genes1174[k, 'num_pc_variants'] <- length(which(vartypes %in% pc_var_types))
  cand_genes1174[k, 'num_reg_variants'] <- length(which(vartypes %in% reg_var_types))
  cand_genes1174[k, 'vartypes'] <- paste(vartypes, collapse = ",")
  
  ix = ix+1
  cand_gene_sum1174[ix, 'gene_name'] <- gene
  if (length(which(vartypes %in% pc_var_types)) > 0){ 
    cand_gene_sum1174[ix, 'pc_vars'] <- TRUE 
  } else {
    cand_gene_sum1174[ix, 'pc_vars'] <- FALSE
  }
}

cand_gene_sum1174$category <- "CC011xCC074"
cand_genes1174$from_which_cross <- "CC011xCC074"


#------------------------------------save--------------------------------------#
write_csv(cand_genes, file = "results/filtered-cand-genes.csv")
all_cand_genes <- rbind(cand_genes, cand_genes644, cand_genes1174)
write_csv(all_cand_genes, file = "results/candidate-genes.csv")

# candidate gene summary 
cand_gene_sum <- rbind(cand_gene_sum, cand_gene_sum644, cand_gene_sum1174)
cand_gene_sum <- as_tibble(cand_gene_sum)

cand_gene_sum %>% 
  group_by(category, pc_vars) %>% 
  dplyr::summarize(count = n())

tst <- cand_gene_sum$gene_name[which(cand_gene_sum$category == 'CC011xCC074' & cand_gene_sum$pc_vars == TRUE)]
paste(tst, collapse = ",")



