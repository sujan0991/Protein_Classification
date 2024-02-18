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

df_Dat <- fread("/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/lev_TM_results_all_safe_cath_after_ck1-2-3.csv",  nThread=120, select = c("SS_score","AA_score","TM_min","key_id"))

ids <- str_split_fixed(df_Dat$key_id, '-', 2)
df_Dat$id1 = ids[,1]
df_Dat$id2 = ids[,2]

cath = fread("/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/csv_folder/cath-domain-description-v4_3_0.csv", nThread=120, select = c("id", "Class", "Arch","Topol","Homol"))
#df_Dat <- df_Dat[sample(nrow(df_Dat), 5000), ]

df_Dat=merge(df_Dat, cath, by.x = "id1", by.y = "id", all = FALSE)
df_Dat=merge(df_Dat, cath, by.x = "id2", by.y = "id", all = FALSE)




## data  prepair -homo
# index_homo = df_Dat$Homol.x == df_Dat$Homol.y
# index_homo_not = df_Dat$Homol.x != df_Dat$Homol.y
# df_Dat$cath_homo =0
# df_Dat$cath_homo[index_homo]=1

# ## data  prepair -class
# index_class = df_Dat$Class.x == df_Dat$Class.y 
# index_class_not = df_Dat$Class.x != df_Dat$Class.y
# df_Dat$cath_class =0
# df_Dat$cath_class[index_class]=1


## data  prepair -Arch 
index_Arch = df_Dat$Arch.x == df_Dat$Arch.y 
index_Arch_not = df_Dat$Arch.x != df_Dat$Arch.y
df_Dat$cath_Arch =0
df_Dat$cath_Arch[index_Arch]=1

# data  prepair -Topol
index_Topol = df_Dat$Topol.x == df_Dat$Topol.y 
index_Topol_not = df_Dat$Topol.x != df_Dat$Topol.y
df_Dat$cath_Topol =0
df_Dat$cath_Topol[index_Topol]=1

#df_Dat <- df_Dat[sample(nrow(df_Dat), 5000), ]

# print('ss_fig1')
# ### plot 1
# df_Dat$TM_min_level = round(df_Dat$TM_min,1)
# svg('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/ggplot_tmLevels-ss_fig1.svg')
# ggplot(
#   df_Dat, 
#   aes(y = as.factor(TM_min_level), x =SS_score, fill = stat(x))) + geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
#   scale_fill_viridis_c(name = "SS_score", direction = 1,option = "C") +
#   labs(title = 'SS vs TM_min') 
# dev.off()


# print('homo-figur2')
# ############ HOMO START
# ### plot2
# df1= df_Dat[,c('SS_score','cath_homo')]
# df2= df_Dat[,c('TM_min','cath_homo')]
# df1$type='SS_score'
# df2$type='TM_score'
# colnames(df1)[1]='data'
# colnames(df2)[1]='data'
# df3=rbind(df1,df2)


# svg('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/ridgeline_tmLevels-ss_homo-figur2.svg')
# ggplot(df3,
#   aes(
#     x = data, y = interaction(cath_homo, type), fill = 0.5 - abs(0.5-stat(ecdf))
#   )) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE)+ 
#   scale_fill_viridis_c(name = "Tail probability", direction = 1, option = "C") +
#   theme_ridges()
# dev.off()

# print('homo-figur3')
# svg('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/ridgeline_tmLevels-ss_homo-figur3.svg')
# ggplot(df3,
#   aes(
#     x = data, y = cath_homo, group = interaction(cath_homo, type), fill = 0.5 - abs(0.5-stat(ecdf))
#   )) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE)+ 
#   scale_fill_viridis_c(name = "Tail probability", direction = 1, option = "C") +
#   theme_ridges()
# dev.off()


# print('homo-figur4')
# ### plot 3
# svg('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/ridgeline_tmLevels-ss_homo-figur4.svg')
# ggplot(df3, aes(x = data, y = interaction(cath_homo, type), fill = factor(stat(quantile))))+
#   stat_density_ridges(
#     geom = "density_ridges_gradient",
#     calc_ecdf = TRUE,
#     quantiles = 5
#   )+
#   scale_fill_viridis_d(name = "Quintiles", direction = 1, option = "C") +
#   theme_ridges()
# dev.off()

################## HOMO END






############ CLASS START

# print('CLASS-figur2')

# ### plot2
# df1= df_Dat[,c('SS_score','cath_class')]
# df2= df_Dat[,c('TM_min','cath_class')]
# df1$type='SS_score'
# df2$type='TM_score'
# colnames(df1)[1]='data'
# colnames(df2)[1]='data'
# df3=rbind(df1,df2)


# svg('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/ridgeline_tmLevels-ss_CLASS-figur2.svg')
# ggplot(df3,
#   aes(
#     x = data, y = interaction(cath_class, type), fill = 0.5 - abs(0.5-stat(ecdf))
#   )) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE)+ 
#   scale_fill_viridis_c(name = "Tail probability", direction = 1, option = "C") +
#   theme_ridges()
# dev.off()

# print('CLASS-figur3')
# svg('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/ridgeline_tmLevels-ss_CLASS-figur3.svg')
# ggplot(df3,
#   aes(
#     x = data, y = cath_class, group = interaction(cath_class, type), fill = 0.5 - abs(0.5-stat(ecdf))
#   )) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE)+ 
#   scale_fill_viridis_c(name = "Tail probability", direction = 1, option = "C") +
#   theme_ridges()
# dev.off()


# print('CLASS-figur4')
# ### plot 3
# svg('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/ridgeline_tmLevels-ss_CLASS-figur4.svg')
# ggplot(df3, aes(x = data, y = interaction(cath_class, type), fill = factor(stat(quantile))))+
#   stat_density_ridges(
#     geom = "density_ridges_gradient",
#     calc_ecdf = TRUE,
#     quantiles = 5
#   )+
#   scale_fill_viridis_d(name = "Quintiles", direction = 1, option = "C") +
#   theme_ridges()
# dev.off()


########### CLASS END





############ CLASS cath_Arch

print('Arch-figur2')

### plot2
df1= df_Dat[,c('SS_score','cath_Arch')]
df2= df_Dat[,c('TM_min','cath_Arch')]
df1$type='SS_score'
df2$type='TM_score'
colnames(df1)[1]='data'
colnames(df2)[1]='data'
df3=rbind(df1,df2)


# svg('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/ridgeline_tmLevels-ss_Arch-figur2.svg')
# ggplot(df3,
#   aes(
#     x = data, y = interaction(cath_Arch,type), fill = 0.5 - abs(0.5-stat(ecdf))
#   )) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE)+ 
#   scale_fill_viridis_c(name = "Tail probability", direction = 1, option = "C") +
#   theme_ridges()
# dev.off()

# print('Arch-figur3')
# svg('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/ridgeline_tmLevels-ss_Arch-figur3.svg')
# ggplot(df3,
#   aes(
#     x = data, y = cath_Arch, group = interaction(cath_Arch,type), fill = 0.5 - abs(0.5-stat(ecdf))
#   )) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE)+ 
#   scale_fill_viridis_c(name = "Tail probability", direction = 1, option = "C") +
#   theme_ridges()
# dev.off()


print('Arch-figur4')
### plot 3
svg('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/ridgeline_tmLevels-ss_Arch-figur4.svg')
ggplot(df3, aes(x = data, y = interaction(cath_Arch,type), fill = factor(stat(quantile))))+
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = 5
  )+
  scale_fill_viridis_d(name = "Quintiles", direction = 1, option = "C") +
  theme_ridges()
dev.off()


########### cath_Arch END


############ CLASS cath_Topol

print('Topol-figur2')

### plot2
df1= df_Dat[,c('SS_score','cath_Topol')]
df2= df_Dat[,c('TM_min','cath_Topol')]
df1$type='SS_score'
df2$type='TM_score'
colnames(df1)[1]='data'
colnames(df2)[1]='data'
df3=rbind(df1,df2)




print('Topol-figur4')
### plot 3
svg('/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/SS_article_plots/ridgeline_tmLevels-ss_Topol-figur4.svg')
ggplot(df3, aes(x = data, y = interaction(cath_Topol, type), fill = factor(stat(quantile))))+
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = 5
  )+
  scale_fill_viridis_d(name = "Quintiles", direction = 1, option = "C") +
  theme_ridges()
dev.off()


########### cath_Arch END









