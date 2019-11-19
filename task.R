BiocManager::install("biomaRt")
BiocManager::install("RCy3")
install.packages("igraph")
install.packages("plyr")
#library(igraph)
#library(plyr)
library(biomaRt)
library(RCy3)
library(dplyr)
#setting working directory
project <- ('C:/Users/Yousuf/Desktop/project')
setwd(project)
# reading file having variant list as a table
my_data <- read.table("variant_list.txt")
#selecting mart
mart.snp <- useMart("ENSEMBL_MART_SNP", "hsapiens_snp")
listDatasets(mart.snp)
listAttributes(mart.snp)
listFilters(mart.snp)
# extracting information related to variants
getHGNC2ENSG = getBM(attributes=c('refsnp_id','ensembl_gene_stable_id','chr_name'),
                     filters ='snp_filter', values = my_data, mart = mart.snp)
#Connect to cytoscape
cytoscapePing ()

cytoscapeVersionInfo ()
# setting up the default layout for the network

setNodeShapeDefault ('ELLIPSE')
setNodeColorDefault ('#AAFF88')
setNodeSizeDefault  (30)
setNodeFontSizeDefault (15)
# define the entities that will be  nodes in your network
nodes <- data.frame(id=union(getHGNC2ENSG$refsnp_id,getHGNC2ENSG$ensembl_gene_stable_id),
                    stringsAsFactors=FALSE)

edges <- data.frame(source=getHGNC2ENSG$refsnp_id,
                    target=getHGNC2ENSG$ensembl_gene_stable_id,  
                    stringsAsFactors=FALSE)

createNetworkFromDataFrames(nodes,edges, title="my first network", collection="DataFrame Example")
full.path=paste(getwd(),'SNP_GENES',sep='/')
exportImage(full.path, 'PNG',zoom=400) #.png scaled by 200%
exportImage(full.path, 'PDF') #.pdf
setVisualStyle('default')
if("cytargetlinker" %in% commandsHelp("")) print("Success: the CyTargetLinker app is installed") else print("Warning: CyTargetLinker app is not installed. Please install the CyTargetLinker app before proceeding.")

genes <- getHGNC2ENSG$ensembl_gene_stable_id
wp <- file.path(getwd(), "LinkSets", "wikipathways-20190610-hsa.xgmml")
CTLextend.cmd = paste('cytargetlinker extend idAttribute="shared name" linkSetFiles="', wp, '" network=current direction=SOURCES', sep="")
commandsRun(CTLextend.cmd)
layoutNetwork()
filter1.cmd = "network select edgeList=all"
filter2.cmd = "network select extendEdges=true"
filter3.cmd = "network create nodeList=selected edgeList=selected networkName=selection source=current"
commandsRun(filter1.cmd)
commandsRun(filter2.cmd)
commandsRun(filter3.cmd)
full.path=paste(getwd(),'final_image',sep='/')
exportImage(full.path, 'PNG',zoom=400) #.png scaled by 200%
exportImage(full.path, 'PDF') #.pdf
?exportImage
save.image()
