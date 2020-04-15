library(rmarkdown)
myArgs <- commandArgs(trailingOnly = TRUE)

lst = list(
   scfilelist = myArgs[2], #list of sc files specified at rendering
   bulkfile =  myArgs[1],
   hbulk = './data/MIX3cl_scESET.rds',
   conditionvar = 'hba1c_class2', #variable where the main condition is of interest
   condition_level1 = 'Normal',
   condition_level2 = 'T2D',
   celltypevar =  myArgs[3], #variable name containing the cell type annot in @phenoData of the eset
   celltypesel=  c("alpha","beta","delta","gamma","acinar","ductal"), #cell types of interest to estimate
   samplevar = "sample" #variable name in @phenoData identifying the sample name

)

render('SCDC_tkt_simplified.Rmd', params=lst)
