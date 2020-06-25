#Seger simulated using seger based on MUSIC
render('script_music_run_extractgep.Rmd', params=list(scfilelist="c('../../docs/datasets/scdc/pancreatic/segerstolpe_raw_eset.RDS')", celltypevar='cluster', expression_path="../../docs/datasets/new/simulated_bulk/segerstolpe_pancreas/simulated_blk_segerstolpe_hugo.csv",name_temp_file='MUSIC_justSeger.txt',samplevar='SubjectName'))


#Seger simulated using just seger based on SCDC
render('~/bescape/docs/datasets/scdc/SCDC_tkt_simplified_SimBulk_justSeger_gepextract.Rmd', params=list(scfilelist="c('./pancreatic/segerstolpe_raw_eset.RDS')", celltypevar='cluster', expression_path="../new/simulated_bulk/segerstolpe_pancreas/simulated_blk_segerstolpe_hugo.csv",name_temp_file='Just_Seger_SimBulk.txt'))
