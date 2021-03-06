render('SCDC_tkt_simplified_withPeng.Rmd', params=list(scfilelist="c('./pancreatic/segerstolpe_raw_eset.RDS','./pancreatic/baron_raw_eset.RDS','./pancreatic/peng_raw_eset.RDS')", celltypevar='cluster', celltypesel="c("pancreatic acinar cell","pancreatic ductal cell","type B pancreatic cell")", bulkfile='./input/fadista_77.rds', name_temp_file='All_temp.txt'))


render('SCDC_tkt_simplified_withBaronSeger.Rmd', params=list(scfilelist="c('./pancreatic/segerstolpe_raw_eset.RDS','./pancreatic/baron_raw_eset.RDS')", celltypevar='cluster', celltypesel="c('pancreatic A cell','pancreatic acinar cell','pancreatic D cell','pancreatic ductal cell','type B pancreatic cell','PP cell')",bulkfile='./input/fadista_77.rds', name_temp_file='Just_Barron_Segerstolpe_temp.txt'))

render('SCDC_tkt_simplified_justPeng.Rmd', params=list(scfilelist="c('./pancreatic/peng_raw_eset.RDS')", celltypevar='cluster', bulkfile='./input/fadista_77.rds', name_temp_file='Just_Peng_temp.txt'))


#Seger simulated using just seger
render('SCDC_tkt_simplified_SimBulk_justSeger.Rmd', params=list(scfilelist="c('./pancreatic/segerstolpe_raw_eset.RDS')", celltypevar='cluster', expression_path="../new/simulated_bulk/segerstolpe_pancreas/simulated_blk_segerstolpe_hugo.csv",name_temp_file='Just_Seger_SimBulk.txt'))


#Seger simulated using baron and  seger
render('SCDC_tkt_simplified_SimBulk_BaronSeger.Rmd', params=list(scfilelist="c('./pancreatic/segerstolpe_raw_eset.RDS','./pancreatic/baron_raw_eset.RDS')", celltypevar='cluster', expression_path="../new/simulated_bulk/segerstolpe_pancreas/simulated_blk_segerstolpe_hugo.csv",name_temp_file='BaronSeger_SimBulk.txt'))


#Seger simulated using peng and baron and  seger
render('SCDC_tkt_simplified_SimBulk_PengBaronSeger.Rmd', params=list(scfilelist="c('./pancreatic/peng_raw_eset.RDS','./pancreatic/segerstolpe_raw_eset.RDS','./pancreatic/baron_raw_eset.RDS')", celltypevar='cluster', expression_path="../new/simulated_bulk/segerstolpe_pancreas/simulated_blk_segerstolpe_hugo.csv",name_temp_file='BaronSeger_SimBulk.txt'))
