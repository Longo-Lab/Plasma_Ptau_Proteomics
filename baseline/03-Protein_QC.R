#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 03-Protein_QC.R
##
## Version: 0.0.2
##
## Purpose of script: Prepare LnPlasmaPTau data for gene set analyses
##
## Author: Robert B Butler III
##
## Date Created: 2022-12-21
##
## Copyright (c) 2022
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   Post linear regression, filters based on protein qualities, dropping 
##   non-human and non-protein markers. Then matches to Ensembl gene IDs.
##
##   For spots detecting gene complexes it splits the genes to duplicate rows, 
##   then for each gene, it selects the row with the highest effect size by 
##   absolute value of beta. Interesting genes should be checked to see if they
##   were detected across multiple spots.
##
##   Finally, prints some plots, and calculates some lists of genes for 
##   downstream analysis, including a gprofiler enrichment check using whole 
##   SomaLogic set as custom background
##
##   REMEMBER, WHEN COMMENTING IN CODE, DO NOT REFERENCE INDIVIDUAL INFORMATION:
##   SAMPLES/PIDNS, DATES OR ANY INDIVIDUAL LEVEL INFO
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(data.table)
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(gprofiler2)
library(patchwork)
library(qqman)

# directory locations
rootdir = '/mnt/box/plasma_ptau_proteomics'
wk_dir = paste0(rootdir, '/baseline')
setwd(wk_dir)
nameset = "dataProt_SS-205063.pTau.baseline"
ensver = "109"

# Functions -------------------------------------
convertGeneSymbol = function(dt, colname, keep.symbols=F){
  # for a given column of NCBI IDs, replace them with ensembl gene ids 
  # includes a filter for autosomal genes that are not pseudo or small RNAs
  # dt: data table to replace names
  # colname: name of the column containing NCBI gene IDs
  # keep.symbols: boolean, retain the symbols column?
  # returns dt with a GENE column
  require("biomaRt")
  print(paste(length(unique(dt[[colname]])), "genes to look up"))
  ensembl = useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", version=ensver)
  # match with ensembl_gene_id across both
  x = data.table(getBM(attributes=c("ensembl_gene_id", "entrezgene_id",
                                    "gene_biotype", "chromosome_name"),
                       filters="entrezgene_id", 
                       values=unique(dt[[colname]]), 
                       mart=ensembl))
  print(paste(nrow(x), "total genes returned"))
  # ensembl transcript types to keep, excluding pseudogenes and TECs
  keep_categ = grep("pseudo", unique(x$gene_biotype), value=T, invert=T)
  keep_categ = grep("TEC", keep_categ, value=T, invert=T)
  x = x[gene_biotype %in% keep_categ]
  print(paste(nrow(x), "genes correct biotype"))
  # use only canonical chromosomes
  use_chrs = c(1:22, "X", "Y", "MT")
  x = x[chromosome_name %in% use_chrs]
  print(paste(nrow(x), "genes in chr 1:22, X, Y, MT"))
  setnames(x, c("ensembl_gene_id", "entrezgene_id"), c("GENE", "NCBI_ID"))
  # replace gene names and write to file
  dt = merge(x[, .(GENE, NCBI_ID)], dt, all.y=T, by.x="NCBI_ID", by.y=colname)
  dt = dt[!is.na(GENE)]
  print(paste(length(unique(dt[["GENE"]])), "unique genes with ensembl ids"))
  if (keep.symbols == FALSE) dt[, NCBI_ID := NULL]
  return(unique(dt))
}

get_mcols = function(str, dds, col) {
  # for a string of ens ids: explode, look up symbols, return symbols as string
  # string: comma separated ens_ids
  # dds: DESeq Dataset with annotations in mcol
  # col: colname to grab
  require("stringr")
  # x = str_split(str, ",")
  y = dds[[col]][match(unlist(str_split(str, ",")), dds[["GENE"]])]
  return(paste0(y, collapse=","))
}


# Main ------------------------------------------
# Load assocaition table
a = fread(paste0(nameset, '.linreg_glm.csv.gz'))
a = a[ix != ""]
# a$variable = gsub("^X\\.", "\\.", a$variable)

# map to reference key
b = fread("../Somalogic_reference_Key.txt")
b = b[Organism == "Human" & Type == "Protein"]
b[, c("Organism", "Units", "Type", "EntrezGeneSymbol_v2") := NULL]
c = merge(a, b, by.x="ix", by.y="Key_2")

# padjust and clean up (total set of human protein probes)
c[, padj := p.adjust(pval, method = "fdr")]
c[, pbonf := p.adjust(pval, method = "bonferroni")]
c = c[EntrezGeneSymbol != ""]
# missing annotations
c[EntrezGeneSymbol == "SHANK3"]$EntrezGeneID = 85358
c[EntrezGeneSymbol == "PAX4"]$EntrezGeneID = 5078
# fix piped rows of multiple NCBI IDs (explode to multiple rows)
c$EntrezGeneID = gsub("\\s+", "|", c$EntrezGeneID, perl=T)
# c$EntrezGeneID = sub("\\|.*$", "", c$EntrezGeneID) # take only the first NCBI_ID in the row
cols = names(c)[names(c) != "EntrezGeneID"]
c = c[, .(EntrezGeneID = tstrsplit(EntrezGeneID, '\\|')), by=cols]
c$EntrezGeneID = as.integer(c$EntrezGeneID)

# drop lower duplicates, use absolute value of beta
# c[duplicated(EntrezGeneSymbol) | duplicated(EntrezGeneSymbol, fromLast=T)][order(EntrezGeneSymbol, -beta)]
c = c[order(-abs(beta))][, head(.SD, 1), by="EntrezGeneID"]

# lookup NCBI IDs for matching
c = convertGeneSymbol(c, "EntrezGeneID", keep.symbols=T)

# plot histogram
c_cut = c[padj < 0.05]
dens = density(c_cut$beta)
qs = quantile(c_cut$beta, c(0.1, 0.9))
dt = data.frame(x=dens$x, y=dens$y)
dt$quant = factor(findInterval(dt$x, c(qs[[1]], 0, qs[[2]])))
# all betas
p1 = ggplot(c, aes(x=beta)) + geom_density(fill="#999999") + 
  ggtitle("All proteins") + theme_classic() + xlim(-0.5, 0.5) +
  theme(plot.title=element_text(hjust=0.5))
# Bonferroni significant betas
p2 = ggplot(dt, aes(x,y)) + geom_line() + xlab("beta") + ylab("density") + 
  geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) + 
  scale_fill_brewer(type="div", palette="RdYlBu", direction=-1, name="",
                    labels=c("Bottom 10%", "Less than 0", "Greater than 0", "Top 10%")) +
  ggtitle("Proteins with Padj < 0.05") + theme_classic() +
  theme(plot.title=element_text(hjust=0.5, face="bold", size=20),
        legend.position=c(0.9, 0.5))
# calculate genomic inflation
lambda.c = median(qchisq(c$pval, 1, lower.tail=FALSE))/qchisq(0.5,1)
# plot QQ
p3 = ggplot(c, aes(sample=-log10(pval))) + 
  stat_qq(distribution=qexp, dparams=list(rate=log(10))) + 
  geom_abline(col="red", intercept=0, slope=1) + 
  ggtitle(sprintf("QQ (Lambda: %.3f)", lambda.c)) +
  theme_classic() + theme(plot.title=element_text(hjust=0.5))
ggsave(plot=p2 + 
         inset_element(p1, left=0.05, top=1, bottom=0.6, right=0.3) +
         inset_element(p3, left=0.7, top=1, bottom=0.6, right=0.95), 
       filename=paste(nameset, "beta_distributions.pdf", sep="."),
       width=10, height=6, device=cairo_pdf)

# write gene lists
sets = list(all = c[padj < 0.05, .(GENE, beta)][order(-abs(beta))],
            up = c[padj < 0.05 & beta > 0, .(GENE, beta)][order(-beta)],
            down = c[padj < 0.05 & beta < 0, .(GENE, beta)][order(beta)],
            top10 = c[padj < 0.05 & beta > qs[[2]], .(GENE, beta)][order(-beta)],
            bottom10 = c[padj < 0.05 & beta < qs[[1]], .(GENE, beta)][order(beta)],
            top.bonf = c[pbonf < 0.05 & beta > 0, .(GENE, beta)][order(-beta)],
            bottom.bonf = c[pbonf < 0.05 & beta < 0, .(GENE, beta)][order(beta)])
lapply(seq(sets), function(i) {
  fwrite(sets[[i]], file=paste(nameset, names(sets)[i], "txt", sep="."), 
         sep="\t", col.names=F)
})
fwrite(c, file=paste(nameset, "csv", sep="."))

# gprofiler enrichment
g = gost(query=sets, organism="hsapiens", ordered_query=T, evcodes=T, 
         custom_bg=c$GENE)
# write to file
if(is.list(g$result)){
  # get symbols for terms
  g$result$symbols = lapply(g$result$intersection, get_mcols, c, "EntrezGeneSymbol")
  # write to file
  fwrite(g$result, paste0(nameset, ".gprofiler.csv"))
}



