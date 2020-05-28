library(rmarkdown)
myArgs <- commandArgs(trailingOnly = TRUE)

lst = list(
   scfilelist = myArgs[1], #list of sc files specified at rendering
   bulkfile =  myArgs[2],
   name_temp_file = myArgs[3],
   hbulk = './data/MIX3cl_scESET.rds',
   conditionvar = 'hba1c_class2', #variable where the main condition is of interest
   condition_level1 = 'Normal',
   condition_level2 = 'T2D',
   celltypevar =  myArgs[4], #variable name containing the cell type annot in @phenoData of the eset
   celltypesel=  myArgs[5], #cell types of interest to estimate
   samplevar = myArgs[6] #variable name in @phenoData identifying the sample name

)



render('SCDC_tkt_simplified.Rmd', params=lst)
