# March 2019 - September 2021
# Differential accessibility analysis of mouse CD8 T cells treated
# with DMSO or UK5099 (mitochondrial pyruvate carrier inhibitor)


# load libraries
library(GenomicRanges)
library(DiffBind)
library(ChIPpeakAnno)
library(calibrate)
library(org.Mm.eg.db)
library(reactome.db)
library(rGREAT)
library(limma)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicFeatures)
library(org.Mm.eg.db)
library(refGenome)
library(ggrepel)
library(clusterProfiler)

# set working directory
# setwd()

# Created a .csv sample sheet to indicate sample names and their associated broadpeak and bam files.
samples_2<-read.csv("samplesheet_2.csv")
# samples_2
#         SampleID  Tissue Factor Condition Treatment Replicate
# 1   Mouse_1_DMSO T_cells   None      DMSO      DMSO         1
# 2   Mouse_2_DMSO T_cells   None      DMSO      DMSO         2
# 3   Mouse_3_DMSO T_cells   None      DMSO      DMSO         3
# 4 Mouse_1_UK5099 T_cells   None    UK5099    UK5099         1
# 5 Mouse_3_UK5099 T_cells   None    UK5099    UK5099         3
#                                                                                                      bamReads
# 1   /Volumes/DOF_TW/09_mathias_ATACseq/alignments/output_Mouse_1_DMSO_S1_R1_001_2_sortedFMp_markdupRG_chr.bam
# 2   /Volumes/DOF_TW/09_mathias_ATACseq/alignments/output_Mouse_2_DMSO_S2_R1_001_2_sortedFMp_markdupRG_chr.bam
# 3   /Volumes/DOF_TW/09_mathias_ATACseq/alignments/output_Mouse_3_DMSO_S3_R1_001_2_sortedFMp_markdupRG_chr.bam
# 4 /Volumes/DOF_TW/09_mathias_ATACseq/alignments/output_Mouse_1_UK5099_S4_R1_001_2_sortedFMp_markdupRG_chr.bam
# 5 /Volumes/DOF_TW/09_mathias_ATACseq/alignments/output_Mouse_3_UK5099_S6_R1_001_2_sortedFMp_markdupRG_chr.bam
#   ControlID bamControl                                          Peaks PeakCaller
# 1        NA         NA macs_callpeak/broad_DMSO_bampe_peaks.broadPeak       narrow
# 2        NA         NA macs_callpeak/broad_DMSO_bampe_peaks.broadPeak       narrow
# 3        NA         NA macs_callpeak/broad_DMSO_bampe_peaks.broadPeak       narrow
# 4        NA         NA macs_callpeak/broad_UK5099_bampe_peaks.broadPeak     narrow
# 5        NA         NA macs_callpeak/broad_UK5099_bampe_peaks.broadPeak     narrow

# Create DBA object
atac_seq_2<-dba(minOverlap=2, sampleSheet = samples_2, peakFormat = "narrow")

atac_seq_2<-dba.count(atac_seq_2,minOverlap = 2,score=DBA_SCORE_TMM_READS_FULL)
# save(atac_seq_2,file="atac_seq_2.RData")
# load("data/atac_seq_2.RData")

# DA analysis:
atac_seq_DE_2<-dba.contrast(atac_seq_2, categories = DBA_CONDITION, minMembers = 2)
atac_seq_DE_2<-dba.analyze(atac_seq_DE_2,method=DBA_EDGER, bSubControl = F)

# Count the number of DA peaks:
# Regions less open in DMSO and more open in UK5099:
length(which(atac_seq.DEpeaks_2$Fold<0))
# [1] 981

# Regions more open in DMSO and less open in UK5099:
length(which(atac_seq.DEpeaks_2$Fold>0))
# [1] 503

# Retrieve DA peaks:
atac_seq.DEpeaks_2<-dba.report(atac_seq_DE_2,contrast =1, method= DBA_EDGER)
# extract all peaks, also the none-significant ones:
atac_seq.DEpeaks_all<-dba.report(atac_seq_DE_2, contrast = 1, method = DBA_EDGER,
                                 th=1)

# Annotate peaks to closest gene:
data("TSS.mouse.GRCm38") # TSS annotation data for Mm obtained from biomaRt
TSS.mouse.GRCm38<-TSS.mouse.GRCm38[which(seqnames(TSS.mouse.GRCm38) %in% 
                                           seqnames(atac_seq.DEpeaks_2))]
macs.anno_2 <- annotatePeakInBatch(atac_seq.DEpeaks_2, AnnotationData=TSS.mouse.GRCm38)
macs.anno_2 <- addGeneIDs(annotatedPeak=macs.anno_2, 
                          orgAnn="org.Mm.eg.db", 
                          IDs2Add="symbol")

# Export DA regions as bed files for Homer motif discovery:
tmp<-cbind(as.character(seqnames(macs.anno_2)),
           start(macs.anno_2),end(macs.anno_2),macs.anno_2$Fold, macs.anno_2$FDR)
tmp<-subset(tmp, tmp[,5]<=0.05)

write.table(subset(tmp, tmp[,4]>0)[,1:3], "MAR.bed",
            row.names = F, col.names = F, quote = F, sep="\t")
write.table(subset(tmp, tmp[,4]<0)[,1:3], "LAR.bed",
            row.names = F, col.names = F, quote = F, sep="\t")

# Install Homer, run on terminal:
# export PATH to homer source scripts
# cd /data/homer_TFBS/bin
# perl findMotifsGenome.pl /data/09_mathias_ATACseq/MAR.bed mm10 /data/09_mathias_ATACseq/output_homer_MAR -mask -size given -mset vertebrates 
# perl findMotifsGenome.pl /data/09_mathias_ATACseq/LAR.bed mm10 /data/09_mathias_ATACseq/output_homer_LAR -mask -size given -mset vertebrates 

# Generated fasta sequences for identification of Runx1 motif locations within DA regions,
# use bedtools:
# bedtools getfasta -fi /mnt/biodata/genomes/Mmusculus/mm10/seq/mm10.fa -bed /data/09_mathias_ATACseq/LAR.bed -fo /data/09_mathias_ATACseq/09_mathias_ATACseq_lar.fa

# Searched for Runx1 motif locations in regions more open in UK5099 (LAR.bed)
# Annotate locations of Runx1 motif to closest gene:
# merge with gene symbol from GTF file
# Export as tables, merge with gtf gene_biotype
mouse_gtf = refGenome::ensemblGenome()
refGenome::read.gtf(mouse_gtf, useBasedir = F, 
                    filename="Mus_musculus.GRCm38.95.gtf")
genes_mm10<-mouse_gtf@ev$genes[ ,c("gene_id","gene_name","gene_biotype")]
head(genes_mm10)
#               gene_id     gene_name         gene_biotype
# 1  ENSMUSG00000102693 4933401J01Rik                  TEC
# 4  ENSMUSG00000064842       Gm26206                snRNA
# 7  ENSMUSG00000051951          Xkr4       protein_coding
# 25 ENSMUSG00000102851       Gm18956 processed_pseudogene
# 28 ENSMUSG00000103377       Gm37180                  TEC
# 31 ENSMUSG00000104017       Gm37363                  TEC

# Import FIMO output of Runx1 motif in LAR (regions less open in DMSO and more open
# in UK5099; 09_mathias_ATACseq_lar.fa) using FIMO:
fimo_lar<-read.table("Runx1_MA0002.2_09_lar_fimo.tsv",
                     header = T)
head(fimo_lar)
#     motif_id motif_alt_id            sequence_name start stop strand   score  p.value q.value matched_sequence
# 1 MA0002.2        RUNX1  chr13:44088999-44089376    56   66      + 15.6034 2.49e-07  0.0433      GTCTGTGGTTT
# 2 MA0002.2        RUNX1  chr10:28489060-28489259    60   70      - 15.6034 2.49e-07  0.0433      GTCTGTGGTTT
# 3 MA0002.2        RUNX1  chr13:30785786-30786850    86   96      - 15.6034 2.49e-07  0.0433      GTCTGTGGTTT
# 4 MA0002.2        RUNX1 chr6:119408355-119409262   484  494      + 15.6034 2.49e-07  0.0433      GTCTGTGGTTT
# 5 MA0002.2        RUNX1  chr13:24778694-24779672   639  649      + 15.6034 2.49e-07  0.0433      GTCTGTGGTTT
# 6 MA0002.2        RUNX1  chr10:94848290-94849600   951  961      + 15.6034 2.49e-07  0.0433      GTCTGTGGTTT

# Convert the exact location of the motif to GenomicRanges before annotation to genes
fimo_lar$chr<-sapply(strsplit(as.character(fimo_lar$sequence_name), ":"), '[', 1)
fimo_lar$sequence_name_2<-sapply(strsplit(as.character(fimo_lar$sequence_name), ":"), '[', 2)
fimo_lar$start_seq<-sapply(strsplit(as.character(fimo_lar$sequence_name_2), "-"), '[', 1)
fimo_lar$start_motif<-as.numeric(as.character(fimo_lar$start_seq))+as.numeric(as.character(fimo_lar$start))
fimo_lar$end_motif<-as.numeric(as.character(fimo_lar$start_seq))+as.numeric(as.character(fimo_lar$stop))
nrow(fimo_lar) # 585
# Export to bed for IGV
write.table(cbind(fimo_lar$chr, fimo_lar$start_motif, fimo_lar$end_motif, as.character(fimo_lar$matched_sequence)),
            "fimo_Runx1_TFBS_location_LAR.bed",
            row.names = F, quote = F, sep="\t", col.names = F)

# Convert to genomic ranges:
fimo_lar_GR<-GRanges(seqnames = fimo_lar$chr, 
                     IRanges(fimo_lar$start_motif, end = fimo_lar$end_motif),
                     mcols=fimo_lar)
# Annotate to gene Ensembl and gene symbol
fimo_lar.anno <- annotatePeakInBatch(fimo_lar_GR, AnnotationData=TSS.mouse.GRCm38)

# Add gene symbol:
fimo_lar.anno <- addGeneIDs(annotatedPeak=fimo_lar.anno, 
                            orgAnn="org.Mm.eg.db", 
                            IDs2Add="symbol")

# Export to Excel to view gene list:
tmp<-as.data.frame(cbind(chr=as.character(seqnames(fimo_lar.anno)),
                         start=as.character(start(fimo_lar.anno)), end=as.character(end(fimo_lar.anno)),
                         gene_id=as.character(fimo_lar.anno$feature),
                         insideFeature=as.character(fimo_lar.anno$insideFeature),
                         symbol=as.character(fimo_lar.anno$symbol),
                         sequence_name=as.character(fimo_lar.anno$mcols.sequence_name),
                         matched_sequence=as.character(fimo_lar.anno$mcols.matched_sequence),
                         strand=as.character(fimo_lar.anno$mcols.strand),
                         score=as.character(fimo_lar.anno$mcols.score),
                         p.value=as.character(fimo_lar.anno$mcols.p.value),
                         q.value=as.character(fimo_lar.anno$mcols.q.value)))
# merge with gene symbol from GTF file.
# Export as tables, merge with gtf gene_biotype
head(genes_mm10)
#               gene_id     gene_name         gene_biotype
# 1  ENSMUSG00000102693 4933401J01Rik                  TEC
# 4  ENSMUSG00000064842       Gm26206                snRNA
# 7  ENSMUSG00000051951          Xkr4       protein_coding
# 25 ENSMUSG00000102851       Gm18956 processed_pseudogene
# 28 ENSMUSG00000103377       Gm37180                  TEC
# 31 ENSMUSG00000104017       Gm37363                  TEC
tmp<-merge(tmp, genes_mm10)
write.table(tmp, "FIMO_annotated_Runx1_LAR.xls", sep="\t",
            row.names = F, quote = F)

# Extract the regions with Runx1 motif in LARs to search for other
# enriched TFs using Homer, to check for TFs co-occuring 
# in LARs with a Runx1 motif. Create bed file to use with Homer:
write.table(cbind(chr=sapply(strsplit(as.character(fimo_lar$sequence_name), ":"), '[', 1),
                  start_seq=sapply(strsplit(as.character(fimo_lar$sequence_name_2), "-"), '[', 1),
                  end_seq=sapply(strsplit(as.character(fimo_lar$sequence_name_2), "-"), '[', 2)),
            "LARs_with_Runx1_motif.bed",
            sep="\t", row.names = F, col.names = F, quote = F)

# HOMER:
# perl findMotifsGenome.pl /data/09_mathias_ATACseq/LARs_with_Runx1_motif.bed mm10 /data/09_mathias_ATACseq/output_homer_LAR_Runx1 -mask -size given -mset vertebrates 


# Create volcano plot including all chromatin regions, use
# use annotation to genes by GREAT
gr_great_all<-as.data.frame(cbind(paste("chr",as.character(seqnames(atac_seq.DEpeaks_all)),sep=""),
                                  as.numeric(start(atac_seq.DEpeaks_all)),
                                  as.numeric(end(atac_seq.DEpeaks_all))))
gr_great_all$V2<-as.numeric(as.character(gr_great_all$V2))
gr_great_all$V3<-as.numeric(as.character(gr_great_all$V2))
gr_great_all<-gr_great_all[-c(grep("GL",gr_great_all$V1)),]
gr_great_all<-gr_great_all[-c(grep("chrJ",gr_great_all$V1)),]
gr_great_all<-gr_great_all[-c(grep("chrMT",gr_great_all$V1)),]

great_query_all<-submitGreatJob(gr=gr_great_all, bg=NULL, # bg=bg_great,
                                species="mm10",includeCuratedRegDoms = T, rule= "basalPlusExt",
                                request_interval = 0,
                                version = "3.0.0")

great_query_all<-submitGreatJob(gr_great_all,species="mm10",version = "3.0", request_interval = 0, max_tries=100)
great_gene_association_all <- plotRegionGeneAssociationGraphs(great_query_all)

# Merge genomic region fold change and p-values with GREAT region annotation:
tmp1<-as.data.frame(cbind(chr_start=paste(paste("chr",as.character(seqnames(atac_seq.DEpeaks_all)),sep=""),
                                          start(atac_seq.DEpeaks_all),sep="_"),
                          chr=paste("chr",as.character(seqnames(atac_seq.DEpeaks_all)),sep=""),
                          range=paste(start(atac_seq.DEpeaks_all), end(atac_seq.DEpeaks_all), sep="-"),
                          fold=atac_seq.DEpeaks_all$Fold,
                          p_val=atac_seq.DEpeaks_all$FDR))
tmp2<-as.data.frame(cbind(chr_start=paste(as.character(seqnames(great_gene_association_all)),
                                          start(great_gene_association_all),sep="_"),
                          gene_GREAT_annot=great_gene_association_all$gene,
                          distance_to_TSS=as.numeric(as.character(great_gene_association_all$distTSS))))
tmp2$distance_to_TSS<-as.numeric(as.character(tmp2$distance_to_TSS))
tmp3<-merge(tmp1,tmp2)
tmp3$distance_to_TSS<-as.numeric(as.character(tmp3$distance_to_TSS))

highlight_genes<-read.csv("genes_to_highlight.csv",
                          header=F)

# Generate volcano plots with ggplot2
tmp3$fold<-as.numeric(as.character(tmp3$fold))
tmp3$p_val<-as.numeric(as.character(tmp3$p_val))

tmp3$color_mem<-ifelse(tmp3$gene_GREAT_annot %in% highlight_genes$V1[grep("^MEM", highlight_genes$V2)]&
                         tmp3$p_val<0.05,
                       T, F)

# saveRDS(tmp3, "data/tmp3_09_atac_DE_peaks_all_df.rds")

# ggplot volcano plot of peaks highlighting memory genes:
tmp3_subset<-subset(tmp3, tmp3$gene_GREAT_annot %in% highlight_genes$V1[grep("^MEM", highlight_genes$V2)]&
                      tmp3$p_val<0.05)


ggplot(tmp3, aes(x = -1*fold, 
                 y = -log10(p_val)),
       color="grey87") +
  geom_point(color="grey87")  +
  ggtitle("") +
  theme_bw() +
  geom_text_repel(data = tmp3_subset, 
                  aes(x=-1*fold,y=-log10(p_val),label=gene_GREAT_annot),
                  box.padding = 0.4,
                  size=5) +
  geom_point(data=tmp3_subset, col="dodgerblue2") +
  theme(axis.title = element_text(size = rel(1.25)),
        axis.text = element_text(size = rel(1.15))) +
  theme(plot.title = element_text(size = rel(1.5))) +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank()) +
  theme(plot.title=element_text(hjust=0.5)) +
  scale_x_continuous(name = "log2(fold change) MPCi vs DMSO") +
  scale_y_continuous(name = "-Log10 p-value") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  ggsave("volcano_plots_MEM_ATACseq_17092021.eps",
         device="eps",
         width = 6, height = 6.5)

# ggplot volcano plot of peaks highlighting metabolism genes:
# create color code for genes:
tmp3_subset_metabo<-subset(tmp3, tmp3$gene_GREAT_annot %in% highlight_genes$V1[grep("^METABO", highlight_genes$V3)]&
                             tmp3$p_val<0.05)
tmp3_subset_metabo<-merge(tmp3_subset_metabo, highlight_genes,
                          by.x="gene_GREAT_annot",
                          by.y="V1")
tmp3_subset_metabo$color<-as.character(tmp3_subset_metabo$V2)
tmp3_subset_metabo$color[which(tmp3_subset_metabo$color=="FAO")]<-"brown1"
tmp3_subset_metabo$color[which(tmp3_subset_metabo$color=="glycolysis")]<-"darkmagenta"
tmp3_subset_metabo$color[which(tmp3_subset_metabo$color=="GO")]<-"darkgoldenrod1"
tmp3_subset_metabo$color[which(tmp3_subset_metabo$color=="oxphos")]<-"forestgreen"
tmp3_subset_metabo$color[which(tmp3_subset_metabo$color=="PPP")]<-"cyan2"


ggplot(tmp3, aes(x = -1*fold, 
                 y = -log10(p_val)),
       color="grey87") +
  geom_point(color="grey87")  +
  ggtitle("") +
  theme_bw() +
  geom_text_repel(data = tmp3_subset_metabo, 
                  aes(x=-1*fold,y=-log10(p_val),label=gene_GREAT_annot),
                  box.padding = 0.4,
                  size=5) +
  geom_point(data=tmp3_subset_metabo, col=tmp3_subset_metabo$color) +
  theme(axis.title = element_text(size = rel(1.25)),
        axis.text = element_text(size = rel(1.15))) +
  theme(plot.title = element_text(size = rel(1.5))) +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank()) +
  theme(plot.title=element_text(hjust=0.5)) +
  scale_x_continuous(name = "log2(fold change) MPCi vs DMSO") +
  scale_y_continuous(name = "-Log10 p-value") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  ggsave("volcano_plots_METABO_ATACseq_17092021.eps",
         device="eps",
         width = 6, height = 8)


#  GSEA and barcode plot of memory genes of Dominguez et al:

# library(clusterProfiler)
# this was performed on R v 4.1.0 with clusterProfiler 

mem_dominguez<-read.xlsx("TE MP gene signature - Dominguez Kaech JEM 2015.xlsx")
mem_dominguez<-mem_dominguez$X3
mem_dominguez<-mem_dominguez[-1]
mem_dominguez<-as.data.frame(cbind(term=rep("Memory_dominguez", length(mem_dominguez)),
                                   gene=mem_dominguez))

# only keep 1 value per gene, the highest absolute fold value:
# tmp3<-readRDS("tmp3_09_atac_DE_peaks_all_df.rds")
tmp3<-tmp3[order(abs(tmp3$fold), decreasing = T),]
tmp3<-tmp3[-c(which(duplicated(tmp3$gene_GREAT_annot))),]

# create geneList:
gl<-(-1)*as.numeric(as.character(tmp3$fold))
names(gl)<-as.character(tmp3$gene_GREAT_annot)
length(gl) # 16880
gl<-sort(gl, decreasing = T)

set.seed(1234)
gseaMem<-GSEA(gl, TERM2GENE = mem_dominguez, 
              eps = 1e-50, pvalueCutoff = 1, seed=T)
gseaMem@result
#               ID      Description setSize enrichmentScore      NES      pvalue    p.adjust
# Memory_dominguez Memory_dominguez      44       0.4242936 1.631258 0.004456283 0.004456283
# qvalues rank                   leading_edge
#      NA 2878 tags=45%, list=17%, signal=38%
#        core_enrichment
#  Tnfsf8/Ccr7/Aff3/Il7r/Polr3e/Dusp4/Smox/Dapl1/Tcf7/Id3/St8sia1/Ctla4/Hpgds/Cyp2s1/Dusp6/Clybl/Sell/Crtam/Tha1/Socs3

pdf("barcodePlot_MemoryP_Dominguez_17092021.pdf",
    width = 8, height = 6)
gseaplot(gseaMem, geneSetID = "Memory_dominguez",
         title="Memory_dominguez")
dev.off()



