suppressMessages(library("EPIC"))
suppressMessages(library("Biobase"))
# python handler checks for proper input of arguments
myArgs <- commandArgs(trailingOnly = TRUE)

#mydata <- read.table(myArgs[1], row.names=1, header=TRUE)
path_gep <- myArgs[1]
path_bulk <- myArgs[2]
name_temp_file <- myArgs[3]
path_temp_file <-  paste('./', name_temp_file, sep='')

# load bulk
bulk <- readRDS(file=path_bulk)

if(class(bulk)[1]=="DGEList") {
    print("Converting DGEList class to eset")
    b.eset <- ganalyse::as.eset(bulk)
    b.eset <- exprs(b.eset)
    print("Setting Gene Symbols as featureNames")
    featureNames(b.eset) <- b.eset@featureData@data$GeneSymbol
} else if (class(bulk)[1]=="ExpressionSet"){
    b.eset <- exprs(bulk)
} else {
    warning("Error: input bulk dataset not an ExpressionSet or DGEList")
}


tryCatch(
    {
        print('Performing deconvolution using EPIC')
        if (path_gep=='epic'){
            out <- EPIC(b.eset)
            write.table(out$cellFractions, path_temp_file, sep="\t")
        } else {
            gep.eset <- readRDS(path_gep)
            gep.epic <- list(refProfiles=exprs(gep.eset), sigGenes='')
            out <- EPIC(b.eset, gep.epic)
            write.table(out$cellFractions, path_temp_file, sep="\t")
        }
    },
    error = function(e) {
        return(paste("Error:",e))
    }
)
