suppressMessages(library(tidyverse))
suppressMessages(library(MuSiC))
suppressMessages(library(BisqueRNA))
#suppressMessages(library(ganalyse))
suppressMessages(library(xbioc))
suppressMessages(library(MCMCpack))

myArgs <- commandArgs(trailingOnly = TRUE)
path_gep <- myArgs[1]
gep.eset<-readRDS(file=path_gep)

path_input <- myArgs[2]
input_df <- readr::read_csv(paste0(path_input, "/input.csv"))
path_output <- paste0(myArgs[3], '/predictions_music.csv')

## Extract the names of each dataset
dataset_names <- input_df$dataset.name

## Extract the names of the expression files that use 
## Hugo symbols as gene identifiers
expression_files  <- input_df$hugo.expr.file

## Form the paths of the expression files
expression_paths <- paste0(path_input, '/', expression_files)

## Load sc datasets: martin and brca

do_music <- function(expression_path, dataset_name){
    
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
    
    out.prop <- music_prop(bulk.eset = bulk, sc.eset = gep.eset, clusters = 'cellType', samples='SubjectName')

    # write cell fractions into temp file, which is to be loaded by the python wrapper
    out.frac <- out.prop$Est.prop.weighted
    result_matrix<-out.frac
    result_matrix<-cbind(result_matrix,bulk@phenoData@data[match(rownames(bulk@phenoData@data),rownames(result_matrix)),])

    

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
result_dfs <- purrr::map2(expression_paths, dataset_names, do_music) 

## Combine all results into one dataframe
combined_result_df <- dplyr::bind_rows(result_dfs)


## Write result into output directory
readr::write_csv(combined_result_df, path_output)

    
## export basis vector
basis = music_basis(gep.eset, clusters = 'cellType', samples = 'SubjectName')
saveRDS(basis, file=paste0(myArgs[3],'/music_basis.RDS'))
