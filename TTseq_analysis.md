## Initial Processing of FASTQ files 

FASTQ files were combined into a single file for each direction of read using:

```bash
cat *1.fq.gz > <name>1.fq.gz
cat *2.fq.gz > <name>2.fq.gz
```

FASTQC was run on each file to generate quality control assessments:

```bash
module load FastQC/0.12.1
fastqc *.fq.gz --extract -o analysis/fastqc
```

Reads were trimmed to remove adaptors using TrimGalore:

```bash
module load TrimGalore/0.6.6
trim_galore --paired --cores 16 --adapter $adapter_file --output_dir $output_dir $file1 $file2
```


## Read alignment

SAMtools and Bowtie2 was used to align reads to hg38. A quality score was not used at this stage as many of the genes of interest (HSPs) have many isoforms and quality scores can result in there loss from the data.

```bash
module load samtools/1.6
module load bowtie/2.2.6
bowtie2 -p 16 -x $bowtieindex -1 $file1 -2 $file2 > $shortfile.sam
samtools view -bS $shortfile.sam > $shortfile.bam
rm $shortfile.sam
samtools flagstat $shortfile.bam | tee $shortfile.flagstat.txt
samtools sort $shortfile.bam > $shortfile.sorted.bam
samtools index $shortfile.sorted.bam
```

A drosophila spike in was used, as test the number of reads of drosophila relative to human (generated with fastq screen) were used to produce initial scaling factors:

```bash
module load fastq_screen/0.11.4
fastq_screen --conf fastq_screen.conf $FILE
```

fastq_screen.conf contained locations of genomes of interest.
The drosophila reads were very inconsistant and were not used for further analysis.


## Bigwig generation

Initial bigwigs were generated to view data before size factors were produced, instead RPKM was used to normalise data. A blacklist was used to remove reads from repetitive regions:

```bash
module load BEDTools/2.25.0
module load libpng/1.6.18
module load ucsc/v326
module load python/3.5.2
bamCoverage -bl hg38-blacklist.v2.bed -p 4 -r chr14 --binSize 100 --normalizeUsing RPKM -b $1.sorted.bam -o $1.chr14.100.bw
```

Bigwigs were then viewed and plots were generated on the UCSC browser.


## Read counting

Feature counts was used to count the number of reads of each sample at each gene in gencode v44:

```bash
module load subread/1.5.2
featureCounts -a gencode.v44.annotation.gtf -t gene -o counts.txt *sorted.bam
```

This generates a txt file with the count data for all samples that was then further analysed using R.


## DESEQ2 analysis of count data

Packages loaded:

```R
library(pasilla)
library(tidyverse)
library(tximport)
library(here)
library(janitor)
library(devtools)
library(yaml)
library(SummarizedExperiment)
library(DESeq2)
library(gtools)
```

Count data was imported, the first five columns (chr, start, end, strand, length) were removed, 'sorted.bam' was removed from column names and the data was converted to a matrix:

```R
countdata <- read.table("counts.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)]
colnames(countdata) <- gsub("\\.sorted.bam$", "", colnames(countdata))
countdata <- as.matrix(countdata)
```

Conditions were assigned to the data (Control NHS, Control HS etc.):

```R
(group2 <- factor(c(rep("ctrl_NHS", 2), rep("ctrl_HS", 2), rep("KD_NHS", 2), rep("KD_HS", 2))))
```

A coldata frame was generated for the DESeqDataSet. See ?DESeqDataSetFromMatrix
```
(coldata <- data.frame(row.names=colnames(countdata), group2))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~group2)
```

The end of the ensembl references (ID versions) were removed, Ensembl info was imported, genes with multiple symbols were collapsed and the annotation was added to the count data:

```R
rownames(dds) <- sub("\\..*", "", rownames(dds))

require(biomaRt)
ensembl <- useMart("ensembl")
ensembl = biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)
gene_annotation <-
  getBM(
    attributes = c(
      "gene_biotype",
      "hgnc_symbol",
      "ensembl_gene_id",
      "chromosome_name",
      "start_position",
      "end_position",
      "percentage_gene_gc_content",
      "description"
    ),
    filters = "ensembl_gene_id",
    values = rownames(dds),
    mart = ensembl
  )

gene_annotation_uni<-gene_annotation %>% 
  group_by(ensembl_gene_id) %>%  
  mutate(gene_symbol=paste0(hgnc_symbol,collapse=",")) %>% 
  dplyr::select(-hgnc_symbol) %>%  
  unique()

rd<-data.frame(gene_id=rownames(dds)) %>% 
  left_join(gene_annotation_uni,
            by=c("gene_id"="ensembl_gene_id")) 
rd <-  DataFrame(rd,row.names = rd$gene_id)
rowData(dds)<- rd
```

The rownames were checked and NAs were added where need: 

```R
sum((rownames(rowData(dds))==rownames(counts(dds)==FALSE)))
rowData(dds)$gene_biotype<-factor(rowData(dds)$gene_biotype)

rowData(dds)$gene_biotype <- addNA(rowData(dds)$gene_biotype)
```


### Normalising counts using DESEQ2

Multiple methods of normalisation were tested including using all gene, using a panel of genes that are known to not change upon HS was found to be the best result.

I normalised the count data based on the reads of ~100 genes that are not affected by HS found here:
https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0538-z/tables/4
https://www.gsea-msigdb.org/gsea/msigdb/cards/HSIAO_HOUSEKEEPING_GENES

The list of genes was input gene and used to generate size factors:

```R
dds_orig<-dds
gene_prefix=stringr::str_sub(rownames(counts(dds)),1,4)

hk_genes<-c(
  "ENSG00000075624", #actb
  "ENSG00000111640", #gapdh
  "ENSG00000124942", #ahnak
  "ENSG00000130402", #actn4
  "ENSG00000100083", #gga1
  "ENSG00000167770", #otub1
  "ENSG00000173039", #rela
  "ENSG00000144028", #snrnp200
  "ENSG00000184009", #ACTG1
  "ENSG00000122359", #ANXA11
  "ENSG00000177879", #APS3S1
  "ENSG00000143761", #ARF1 
  "ENSG00000175220", #ARHGAP1
  "ENSG00000163466", #ARPC2
  "ENSG00000152234", #ATP5F1A
  "ENSG00000166710", #B2M
  "ENSG00000185825", #BCAP31
  "ENSG00000126581", #BECN1
  "ENSG00000108561", #C1QBP
  "ENSG00000131236", #CAP1
  "ENSG00000118816", #CCNI
  "ENSG00000135535", #CD164
  "ENSG00000110651", #CD81
  "ENSG00000099622", #CIRBP
  "ENSG00000213719", #CLIC1
  "ENSG00000122705",	#CLTA
  "ENSG00000093010",	#COMT
  "ENSG00000006695",	#COX10
  "ENSG00000126267", #COX6B1
  "ENSG00000112695", #COX7A2
  "ENSG00000160213", #CSTB
  "ENSG00000175203", #DCTN2
  "ENSG00000167986", #DDB1
  "ENSG00000136271", #DDX56
  "ENSG00000125868", #DSTN
  "ENSG00000104529", #EEF1D
  "ENSG00000175390", #EIF3F
  "ENSG00000196924", #FLNA
  "ENSG00000168522", #FNTA
  "ENSG00000169727", #GPS1
  "ENSG00000167468", #GPX4
  "ENSG00000163041", #H3-3A
  "ENSG00000068001", #HYAL2
  "ENSG00000166333", #ILK
  "ENSG00000184216", #IRAK1
  "ENSG00000150093", #ITGB1
  "ENSG00000100605", #ITPK1
  "ENSG00000111144", #LTA4H
  "ENSG00000133030", #MPRIP
  "ENSG00000147065", #MSN
  "ENSG00000196498", #NCOR2
  "ENSG00000147862", #NFIB
  "ENSG00000100906", #NFKBIA
  "ENSG00000143799", #PARP1
  "ENSG00000007372", #PAX6
  "ENSG00000107438", #PDLIM1
  "ENSG00000177700", #POLR2L
  "ENSG00000196262", #PPIA
  "ENSG00000277791", #PSMB3
  "ENSG00000175166", #PSMD2
  "ENSG00000092010", #PSME1
  "ENSG00000156471", #PTDSS1
  "ENSG00000187514", #PTMA
  "ENSG00000184007", #PTP4A2
  "ENSG00000172053", #QARS1
  "ENSG00000157916", #RER1
  "ENSG00000188846", #RPL14
  "ENSG00000063177", #RPL18
  "ENSG00000166441", #RPL27A
  "ENSG00000108107", #RPL28
  "ENSG00000198918", #RPL39
  "ENSG00000118705", #RPN2
  "ENSG00000142534", #RPS11
  "ENSG00000115268", #RPS15
  "ENSG00000134419", #RPS15A
  "ENSG00000233927", #RPS28
  "ENSG00000213741", #RPS29
  "ENSG00000135972", #RPS9
  "ENSG00000168028", #RPSA
  "ENSG00000197747", #S100A10
  "ENSG00000031698", #SARS1
  "ENSG00000087365", #SF3B2
  "ENSG00000075415", #SLC25A3
  "ENSG00000115306", #SPTBN1
  "ENSG00000138385", #SSB
  "ENSG00000148290", #SURF1
  "ENSG00000106052", #TAX1BP1
  "ENSG00000104964", #TLE5
  "ENSG00000139644", #TMBIM6
  "ENSG00000047410", #TPR
  "ENSG00000130726", #TRIM28
  "ENSG00000221983", #UBA52
  "ENSG00000175063", #UBE2C
  "ENSG00000165637", #VDAC2
  "ENSG00000100219", #XBP1
  "ENSG00000128245", #YWHAH
  "ENSG00000167232" #ZNF91
)

dds <- DESeq2::estimateSizeFactors(dds, controlGenes = which(rownames(dds) %in% hk_genes))

sizeFactors(dds)
1.5472687  0.8763560  1.2698893  1.9628579  0.6354695  0.5025626  1.2364784  0.7836489 
```

Deseq normalisation was run using the generated size factors:

```R
dds2<-dds
dds2 <-DESeq(dds2)
plotDispEsts(dds2)
resultsNames(dds2)
```


### Produce comparisons between conditions

A list of comparisons was produced e.g. Ctrl NHS v Ctrl HS etc. and gene annotation and names were added:

```R
contrast_list<-list(c("group2","ctrl_HS","ctrl_NHS"),
                    c("group2","KD_NHS","ctrl_NHS"),
                    c("group2","KD_HS","ctrl_HS"),
                    c("group2","KD_HS","KD_NHS"),
                    c("group2","KD_HS","ctrl_NHS"))

names(contrast_list)<-contrast_list %>% map_chr(function(x){paste0(x[2],"_V_",x[3])})
res_list<-contrast_list %>% map(function(x){results(dds2,contrast=x)})
names(res_list)<-names(contrast_list)

res_list_tidy<-res_list%>% 
  map(function(x){as_tibble(x,rownames="gene_id") %>% 
      left_join(data.frame(rowData(dds2)), by = join_by(gene_id))})

names(res_list_tidy)<-names(contrast_list)

```

Save comparisons and normalised count data

```R
res_dir = here("deseq2_results","DEG",Sys.Date())
if(!file.exists(res_dir)){dir.create(res_dir,recursive = T)}
names(res_list) %>% 
  map(function(x){
    out<-res_list[[x]] %>% 
      as_tibble(rownames="gene_id")  %>% 
      left_join(data.frame(rowData(dds2))[,1:8])
    write_csv(out,file = here(res_dir,paste0(x,"_deg.csv")))
  })

nc <- as_tibble(counts(dds2,T),rownames="gene_id") %>% left_join(data.frame(rowData(dds2))[,1:8])
write_csv(nc,here(res_dir,"deseq2_norm_counts.csv"))
```
