library(magrittr)
library(clusterProfiler)
library(rWikiPathways)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(biomaRt)
library(RCy3)


listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)


filters = listFilters(ensembl)
filters[1:5,]

attributes = listAttributes(ensembl)
attributes[1:5,]


genes = getBM(attributes=c("chromosome_name",
                           "start_position",
                           "end_position",
                           "ensembl_gene_id_version",
                           "external_gene_name",
                           "entrezgene",
                           "strand",
                           "gene_biotype"), 
              filters = c("biotype"), 
              values = list(c("protein_coding")), 
              mart = ensembl)

row.names(genes) = genes$ensembl_gene_id_version
names(genes) = c("chr",'start',"end","ensID","gene_name","entrezID","strand","gene_biotype")



#downloadPathwayArchive(organism="Mus musculus", format="gmt")
#downloadPathwayArchive(organism="Drosophila melanogaster", format="gmt")
#downloadPathwayArchive(organism="Homo sapiens", format="gmt")

d <- read.csv("~/ferrari/PhD_project/reference_datasets/Ferrari_iNPC_DMSOvsEPZ_RNA-Seq/downstream_analysis/DESEq2/output_DESeq2_analysis/DE_genes_noLFCthr_shrinked.tsv",sep="\t")
d$ensID = row.names(d)
## assume that 1st column is ID
## 2nd column is fold change

entrez_de = merge(d,genes)


## feature 1: numeric vector
geneList <- entrez_de$log2FoldChange[(entrez_de$padj < 0.05)  & (entrez_de$log2FoldChange<0)]

## feature 2: named vector
names(geneList) <- as.character(entrez_de$entrezID[(entrez_de$padj < 0.05) & (entrez_de$log2FoldChange<0)] )
geneList = na.omit(geneList)
## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)


#data(geneList, package="DOSE")
gene <- names(geneList)



#wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
#wp2gene <- read.gmt(wpgmtfile)
#wp2gene_1 = read.csv("wikipathways-20190510-gmt-Mus_musculus.gmt",sep="\t", header = F)
#wp2gene_1 = t(wp2gene_1)
#write.table(wp2gene_1, "./transpose_wikipath.gmt",sep="\t", row.names = F, col.names = F,quote = F)

wp2gene = read.csv("wikipath_mouse.gmt",sep="\t")
wp2gene = na.omit(wp2gene)

wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp)


ewp <- setReadable(ewp, org.Mm.eg.db, keyType = "ENTREZID")

head(ewp)

dotplot(ewp)
emapplot(ewp)
heatplot(ewp, foldChange=geneList)
cnetplot(ewp, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
upsetplot(ewp)

### visualize pathways

gpml = getPathway(pathway="WP4")

svg = getColoredPathway(pathway="WP1842", graphId=c("dd68a","a2c17"),
                        color=c("FF0000", "00FF00"));
writeLines(svg, "pathway.svg")


url <- getPathwayInfo("WP179")[2]
browseURL(url)


xrefs = getXrefList(pathway="WP4", systemCode="S")
pathways = findPathwaysByXref("HMDB00001", "Ch")
pathways = findPathwaysByXref(identifier="HMDB00001", systemCode="Ch")
pathways = findPathwaysByXref(
  identifier=c("HMDB00001", "HMDB00002"),
  systemCode=c("Ch", "Ch") 
)





library(msigdbr)
msigdbr_show_species()
m_df <- msigdbr(species = "Mus musculus")
m_t2g <- msigdbr(species = "Mus musculus", category = c("H","C1")) %>% 
  dplyr::select(gs_name, entrez_gene)
em <- enricher(gene, TERM2GENE=m_t2g)
em <- setReadable(em, org.Mm.eg.db, keyType = "ENTREZID")
cnetplot(em, foldChange=geneList, circular = T, colorEdge = TRUE, showCategory=5)
emapplot(em, showCategory=10)
as.data.frame(em)
