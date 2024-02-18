library(ggplot2)
library(ggridges)
theme_set(theme_minimal())
library(RColorBrewer)
library(viridis)
library(plyr)
#library(vcd)
#library(ridgeline)
#library(hexbin)
library(stringr)

#BiocManager::install("stringr")
#install("ridgeline")
#install.packages("ggridges")
#CONTAINER=/group/bioinf_biomarkers_rna/softwares/R3_6.3_ws1.sif
#singularity exec -B /group $CONTAINER R
#module load apps/singularity/3.7.3 
library(data.table)
setDTthreads(threads = 120)


df_Dat <- fread("/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/SCOPE_LEV_HH_blast_results.csv",  nThread=120, select = c('scope_Family',
       'is_same_scope_Superfamily', 'SS_score', 'AA_score', 'key_id', 'tm_min',
       'tm_max'))


df_Dat

df1= df_Dat[,c('SS_score','scope_Family')]
# df2= df_Dat[,c('tm_min','is_same_scope_Superfamily')]
df1$type='SS_score'
# df2$type='TM_score'
colnames(df1)[1]='data'
# colnames(df2)[1]='data'
# df3=rbind(df1,df2)

# print('homo-figur4')
# ### plot 3
svg('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/ridgeline_ss_family-scop.svg')
ggplot(df1, aes(x = data, y = interaction(scope_Family, type), fill = factor(stat(quantile))))+
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = 5
  )+
  scale_fill_viridis_d(name = "Quintiles", direction = 1, option = "C") +
  theme_ridges()
dev.off()

