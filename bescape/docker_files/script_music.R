#library(Seurat)
suppressMessages(library(MuSiC))
suppressMessages(library(BisqueRNA))
suppressMessages(library(ganalyse))
suppressMessages(library(xbioc))
suppressMessages(library(MCMCpack))

myArgs <- commandArgs(trailingOnly = TRUE)
path_gep <- myArgs[1]
path_bulk <- myArgs[2]
name_temp_file <- myArgs[3]
path_temp_file <-  paste('./', name_temp_file, sep='')

# prep data for deconvolution
gep.eset<-readRDS(file=path_gep)

# load bulk
bulk <- readRDS(file=path_bulk)

if(class(bulk)[1]=="DGEList"){
    print("Converting DGEList class to eset")
    b.eset <- ganalyse::as.eset(bulk)

    print("Setting Gene Symbols as featureNames")
    featureNames(b.eset) <- b.eset@featureData@data$GeneSymbol
} else if (class(bulk)[1]=="ExpressionSet"){
    b.eset <- bulk
} else {
    return("Error: input bulk dataset not an ExpressionSet or DGEList")
}

tryCatch(
{
    print("Performing deconvolution using MuSiC")
    # perform deconvolution
    out.prop <- music_prop(bulk.eset = b.eset, sc.eset = gep.eset, clusters = 'cellType', samples='SubjectName')

    # write cell fractions into temp file, which is to be loaded by the python wrapper
    out.frac <- out.prop$Est.prop.weighted
    write.table(out.frac, path_temp_file, sep="\t")
},
    error = function(e){
        return(paste("Error:", e))
}
)

#saveRDS(liver.frac, 'output/DGE_454_liverGEP.RDS')
#save(liver.frac, file='output/DGE_454_liverGEP.Rdata')


