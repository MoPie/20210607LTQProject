setwd('~/Desktop/20210423Litianqing/')

methylation_data = read.csv('./20210504GenomeWideMethylation.csv',
                            header = TRUE, stringsAsFactors = FALSE)


plot_order = c('Blastocyst','WGBS_SSEANeg','H9_Naive','h2B43-ES1','h2B43-ES2','h2B43-ES3','hB434',
               'WGBS_UCLA1_Primed','H9_Primed','PH1','PH9P117','PS4')

plot_order = c('Blastocyst','h2B43-ES1','h2B43-ES2','h2B43-ES3','hB434',
               'PH1','PH9P117','PS4')
rownames(methylation_data) = methylation_data$Sample
methylation_data_order = methylation_data[plot_order,]

plot_color = c('cyan3','cadetblue3','cadetblue3','cadetblue3','cadetblue3',
               'lightsalmon3','lightsalmon3','lightsalmon3')
pdf('20210528GenomeWideCGMethylation.pdf',width = 10, height = 5)
barplot(methylation_data_order$GenomeWideMethylation * 100,
        names.arg = methylation_data_order$Sample,
        las = 2, cex.names = 0.4,border = FALSE,
        col = plot_color,main = 'Methylation(GenomeWide)',
        ylab = 'CG Methylation')
dev.off()



myreadfunction = function(path){
  data = read.csv(path,sep = '\t', header = TRUE,stringsAsFactors = FALSE)
  #other data
  colname = strsplit(basename(path),split = 'CGmeth')[[1]][1]
  # LTQ data
  # colname = strsplit(basename(path),split = '_methratio')[[1]][1]
  data[colname] = data['C_count'] / data['CT_count']
  data[colname] = 100* data[colname]
  return(data[colname])
}


file = list.files(path = './WGBS/',pattern = "*OverStableImprints.txt", full.names = TRUE)
file.l = lapply(file, myreadfunction)

#file
readcount.df = do.call(cbind,file.l)
head(readcount.df)

otherdata_file = list.files(path = './WGBS/',pattern = "*OverStableImprints_hg19.txt", full.names = TRUE)
otherdata_file.l = lapply(otherdata_file, myreadfunction)
otherdata = do.call(cbind, otherdata_file.l)
# remove published data
# otherdata$H9_Primed_average = (otherdata$GSM1493983_H9_Primed1 + otherdata$GSM1493984_H9_Primed2 + otherdata$GSM1493985_H9_Primed3) /3
# otherdata$H9_Naive_average = (otherdata$GSM1493986_H9_Naive1 + otherdata$GSM1493987_H9_Naive2 + otherdata$GSM1493988_H9_Naive3) /3
# otherdata$WGBS_UCLA1_Primed_average = (otherdata$GSM2041690_WGBS_UCLA1_Primed_rep1_ + otherdata$GSM2041691_WGBS_UCLA1_Primed_rep2_ + otherdata$GSM2041692_WGBS_UCLA1_Primed_rep3_) / 3
# otherdata$WGBS_UCLA1_SSEA4Neg_average = (otherdata$GSM2041693_WGBS_SSEA4Neg_rep1_ + otherdata$GSM2041694_WGBS_SSEA4Neg_rep2_ + otherdata$GSM2041695_WGBS_SSEA4Neg_rep3_) / 3
# 
# plotdata = cbind(readcount.df, Blastocyst = otherdata$Blastocyst_, 
#                  H9_Primed_average = otherdata$H9_Primed_average, 
#                  H9_Naive_average = otherdata$H9_Naive_average,
#                  WGBS_UCLA1_Primed_average = otherdata$WGBS_UCLA1_Primed_average, 
#                  WGBS_UCLA1_SSEA4Neg_average = otherdata$WGBS_UCLA1_SSEA4Neg_average)

plotdata = cbind(readcount.df, Blastocyst = otherdata$Blastocyst_)
colnames(plotdata) = gsub('_methratio_OverStableImprints.txt','',colnames(plotdata))
plotdata_order = plotdata[,c(8,seq(1,7))]

#readcount.df['mCG'] = 100* readcount.df['mCG']
# colors = c('#B5B5B5','#9C9C9C','#6CA6CD','#87CEFF','#009ACD',
#            '#98F5FF','#8EE5EE','#FFE4C4','#EED5B7')
# colors = c('cadetblue3','cadetblue3','cadetblue3','cadetblue3','lightsalmon3','lightsalmon3','lightsalmon3')
pdf('20210528mCGoverStablePrimaryImprints.pdf',width = 6,height = 5)
#Outvals = boxplot(readcount.df)$out
#readcount.df.nooutvals = readcount.df[! readcount.df %in% Outvals]
boxplot(plotdata_order,border =plot_color,cex.axis=0.3, par(las="2"),
        main = 'Methylation (Stable primary imprints)',
        ylab = 'CG methylation (%)',outline = TRUE,col = 'white')
dev.off()



### 20210512 Methratio Over X ,20210528 heatmap , 20210530 change promoter region to +-1kb###
library(pheatmap)
library(RColorBrewer)
myreadfunction = function(path){
  data = read.csv(path,sep = '\t', header = TRUE,stringsAsFactors = FALSE)
  #other data
  colname = strsplit(basename(path),split = '_methratio')[[1]][1]
  # LTQ data
  # colname = strsplit(basename(path),split = '_methratio')[[1]][1]
  data[colname] = data['C_count'] / data['CT_count']
  data[colname] = 100* data[colname]
  return(data[colname])
}
### color test
test1 = colorRampPalette(c("#002090", "#D4E675"))(50)
test2 = colorRampPalette(c("#D4E675", "#FE0013"))(50)
test_color = c(test1[1:49],test2)
####

file = list.files(path = './WGBS/20210529Promoter2kb/',pattern = "*Promoter_2kb_CpGIsland_chrX_hg38_filter.txt", full.names = TRUE)
file.l = lapply(file, myreadfunction)
readcount.df = do.call(cbind,file.l)
plotdata_order = readcount.df
plot_color= c('cadetblue3','cadetblue3','cadetblue3','cadetblue3','lightsalmon3','lightsalmon3','lightsalmon3')
plot_color2 = c('lightsalmon3','lightsalmon3','lightsalmon3','lightsalmon3','cadetblue3','cadetblue3','cadetblue3')
plotdata_order[is.na(plotdata_order)] = 0

pdf('./WGBS/20210529Promoter2kb/20210530MethylationOverPromoter2kbCpGIsland_ChrX.pdf',width = 6,height = 5)
pheatmap(plotdata_order,
         cluster_cols = FALSE, 
         color = test_color
         )
dev.off()


pdf('./WGBS/20210529Promoter2kb/20210530MethylationOverPromoter2kbCpGIsland2.pdf',width = 6,height = 5)
boxplot(plotdata_order,
        border =plot_color2,
        cex.axis=0.5, par(las="2"),
        main = 'Methylation (Promoter CpG Island on chrX)',
        ylab = 'CG methylation (%)',outline = TRUE,col = 'white')
dev.off()


non_file = list.files(path = './WGBS/20210529Promoter2kb/',pattern = "*Promoter_2kb_nonCpGIsland_chrX_hg38_filter.txt", full.names = TRUE)
non_file.l = lapply(non_file, myreadfunction)
non_readcount.df = do.call(cbind,non_file.l)
plotdata_order = non_readcount.df
plotdata_order[is.na(plotdata_order)] = 0

### color ###
test1 = colorRampPalette(c("#002090", "#D4E675"))(50)
test2 = colorRampPalette(c("#D4E675", "#FE0013"))(50)
test_color = c(test1[1:49],test2[1:11])
test3 = colorRampPalette(c(test2[11],'#EEC900'))(31)
test4 = colorRampPalette(c('#EEC900',test2[50]))(11)
final_color = c(test_color,test3[2:31],test4[2:11])

### 

pdf('./WGBS/20210529Promoter2kb/20210530Methylation2kbOverPromoter_nonCpGIsland_ChrX.pdf',width = 6,height = 5)
pheatmap(plotdata_order,
         cluster_cols = FALSE, 
         color = final_color
)
dev.off()

pdf('./WGBS/20210529Promoter2kb/20210530MethylationOverPromoter2kbnonCpGIsland2.pdf',width = 6,height = 5)
boxplot(plotdata_order,
        border =plot_color2,
        cex.axis=0.5, par(las="2"),
        main = 'Methylation (Promoter non-CpG Island on chrX)',
        ylab = 'CG methylation (%)',outline = TRUE,col = 'white')
dev.off()


shuffle_file = list.files(path = './WGBS/20210529Promoter2kb/',pattern = "*Random_2kb_outside_CGI_non_CGI_promoter_filter.txt", full.names = TRUE)
file.l = lapply(shuffle_file, myreadfunction)
readcount.df = do.call(cbind,file.l)
plotdata_order = readcount.df
plot_color= c('cadetblue3','cadetblue3','cadetblue3','cadetblue3','lightsalmon3','lightsalmon3','lightsalmon3')

plotdata_order[is.na(plotdata_order)] = 0

pdf('./WGBS/20210529Promoter2kb/20210530MethylationOverRandom_2kb_outside_CGI_non_CGI_promoter.pdf',width = 6,height = 5)
pheatmap(plotdata_order,
         cluster_cols = FALSE, 
         color = final_color
)
dev.off()

pdf('./WGBS/20210529Promoter2kb/20210530MethylationOverRandom_2kb_outside_CGI_non_CGI_promoter_boxplot2.pdf',width = 6,height = 5)
boxplot(plotdata_order,
        border =plot_color2,
        cex.axis=0.5, par(las="2"),
        main = 'Methylation (Promoter CpG Island on chrX Shuffle)',
        ylab = 'CG methylation (%)',outline = TRUE,col = 'white')
dev.off()


non_shuffle_file = list.files(path = './WGBS/MethRatioOverX/',pattern = "*Promoter_nonCpGIsland_chrX_hg38_shuffle_seed666_filter.txt", full.names = TRUE)
file.l = lapply(non_shuffle_file, myreadfunction)
readcount.df = do.call(cbind,file.l)
plotdata_order = readcount.df

plot_color= c('cadetblue3','cadetblue3','cadetblue3','cadetblue3','lightsalmon3','lightsalmon3','lightsalmon3')

plotdata_order[is.na(plotdata_order)] = 0

pdf('./20210528MethylationOverPromoter_nonCpGIsland_ChrX_Shuffle2.pdf',width = 6,height = 5)
pheatmap(plotdata_order,
         cluster_cols = FALSE, 
         color = test_color
)
dev.off()

pdf('./Figure/20210512MethylationOverPromoternonCpGIslandShuffle.pdf',width = 6,height = 5)
boxplot(plotdata_order,
        border =plot_color,
        cex.axis=0.5, par(las="2"),
        main = 'Methylation (Promoter non-CpG Island on chrX Shuffle)',
        ylab = 'CG methylation (%)',outline = TRUE,col = 'white')
dev.off()
