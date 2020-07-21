library(tidyverse)
library(Biobase)
library(SCDC)

## Read in the round and sub-Challenge-specific input file 
## listing each of the datasets
myArgs <- commandArgs(trailingOnly = TRUE)
#Declare parameters needed for scdc(TODO)
scfilestring <- gsub(":", "'", myArgs[1])
scfilelist <- eval(parse(text=scfilestring))

celltypevar = myArgs[5] #variable name containing the cell type annot in @phenoData of the eset
samplevar = myArgs[6] # variable name in @phenoData@data$... identifying sample name

celltypesel_str <- gsub(":", "'", myArgs[4])
celltypesel <- eval(parse(text=celltypesel_str))

path_input <- myArgs[2]
input_df <- readr::read_csv(paste0(path_input, "/input.csv"))

path_output <- paste0(myArgs[3], '/predictions_scdc.csv')

## Extract the names of each dataset
dataset_names <- input_df$dataset.name

## Extract the names of the expression files that use 
## Hugo symbols as gene identifiers
expression_files  <- input_df$hugo.expr.file

## Form the paths of the expression files
expression_paths <- paste0(path_input, "/", expression_files)

## Load sc datasets: martin and brca
scdata = list(0)
nscdata<-length(scfilelist)
for(index in 1:nscdata){
scdata[[index]]<-readRDS(scfilelist[index])
}


do_scdc <- function(expression_path, dataset_name){
    
    # This reads in the input file and converts to a matrix which will be
    # input to SCDC
    expression_matrix <- expression_path %>% 
        readr::read_csv() %>% 
        as.data.frame() %>%
        tibble::column_to_rownames("Gene") %>% 
        as.matrix()

    #creating phenoData with just sample ID (TKT)
    pData<-as.data.frame(colnames(expression_matrix))
    colnames(pData)<-"sample"
    metadata<-data.frame(labelDescription=c("sample"),row.names=c("sample"))
    phenoData<-new("AnnotatedDataFrame",data=pData,varMetadata=metadata)
    row.names(phenoData)<-colnames(expression_matrix)

    bulk<-Biobase::ExpressionSet(assayData=expression_matrix,phenoData=phenoData)#no phenotype data available
    
    ens <- SCDC_ENSEMBLE(bulk.eset = bulk, sc.eset.list = scdata, ct.varname = celltypevar, sample = samplevar, truep = NULL, ct.sub =  celltypesel, search.length = 0.01, grid.search = T)
    
    result_matrix <-  as.data.frame(wt_prop(ens$w_table[1, 1:nscdata], ens$prop.only))
    
    
    
    
    #For SCDC based on current output/predictions.csv(TKT)
    #1) put the sample into column, use as key
    #2) convert wide to long format using all the remaining variables(cell types)   
    result_df_wide<- result_matrix %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("sample.id")
    result_df_wide$sample.id<-as.factor(result_df_wide$sample.id)
    result_df <- tidyr::gather(result_df_wide,cell.type,prediction,-sample.id,factor_key= TRUE) 
    
    
    # Add dataset column
    result_df <- dplyr::mutate(result_df, dataset.name = dataset_name)
}

## Run MCP-Counter on each of the expression files
result_dfs <- purrr::map2(expression_paths, dataset_names, do_scdc) 

## Combine all results into one dataframe
combined_result_df <- dplyr::bind_rows(result_dfs)


## Write result into output directory
readr::write_csv(combined_result_df, path_output)

for(index in 1:nscdata){
    basis = SCDC_basis(scdata[[index]], ct.varname = 'cellType', sample = 'SubjectName')
    saveRDS(basis, file=paste0(myArgs[3],'/scdc_basis_', index, '.RDS'))
    
}
    
