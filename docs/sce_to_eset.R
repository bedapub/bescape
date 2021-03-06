library(patchwork)
library(Seurat)
library(BisqueRNA)


args= commandArgs(trailingOnly=TRUE)

if(length(args) !=3){
    stop("You need to specify sceformat.RDS and dblabel.csv files as input and filename that should be output to")
} else {
    sceinput = args[1]
    annot_csv= args[2]
    output = args[3]
}

sce_to_seurat <- function(sce_path, sc_anno_path, filename='eset.RDS'){
    sce <- readRDS(sce_path)
    sc_anno <- read.csv(sc_anno_path, header=FALSE)

    seurat <- as.Seurat(sce, counts=NULL, data='X')

# TODO add celltype_dream as param
    seurat@meta.data$celltype_dream <- sc_anno$V1
    Idents(seurat)<-"celltype_dream"
    seurat@assays$RNA@counts<-seurat@assays$RNA@data
    if(!is.null(filename)){
        saveRDS(seurat, file=filename)
    }
    return(seurat)

}

seurat_to_eset <- function(seurat, delim='-', idx=2, filename='seurat.RDS'){
    # Idents(liver_h5ad) %>% head() # to check sample names
    out.eset <- BisqueRNA::SeuratToExpressionSet(seurat, delim, idx, version = 'v3')
    # fix GeneID
    out.eset@featureData@data$GeneID <- featureNames(out.eset)
    #assign cluster to just characters
    out.eset@phenoData@data$cluster<-as.character(out.eset@phenoData@data$cellType)
    #update vardata
    varMetadata(out.eset)<-data.frame(labelDescription=colnames(out.eset@phenoData@data),stringsAsFactors=FALSE)

    out.eset@featureData@varMetadata<-data.frame(labelDescription="GeneID",stringsAsFactors=FALSE)

    rownames(out.eset@featureData@varMetadata)<-"GeneID"

    saveRDS(out.eset, file=filename)
    return(out.eset)
}


seurat <- sce_to_seurat(sce_path=sceinput, sc_anno_path=annot_csv, filename=NULL)
eset <- seurat_to_eset(seurat, delim='_', idx=2, filename=output)

