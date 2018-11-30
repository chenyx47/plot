library(pheatmap)
library(cowplot)
annotation<-data.frame(Type=sample)
row.names(annotation)<-names(STAD.fpkm_M)[-1]
anno_colors <- list(Type = c("STAD(n=375)" = "red", "Normal(n=32)" = "black"))

STAD.fpkm_M.matrix<-FeaturetoRow(STAD.fpkm_M)
STAD.fpkm_M.scale<-t(scale(t(STAD.fpkm_M.matrix)))
STAD.fpkm_M.scale1<-STAD.fpkm_M.scale[genes,]

rm_outlier<-function(data){
  logiclow<-data<quantile(data,0.05)
  data[logiclow]<-quantile(data,0.05)[1]
  logichigh<-data>quantile(data,0.95)
  data[logichigh]<-quantile(data,0.95)[1]
  return(data)
}


bk <- c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))

p <- pheatmap(STAD.fpkm_M.scale1,
              # cluster_cols = d$hc,
              # cluster_rows = F,
              # cutree_cols = 3,
              # cutree_rows = 3,
              # gaps_row = c(10, 20),
              annotation_col = annotation,
              annotation_colors = anno_colors,
              show_colnames = F, 
              cluster_cols = F,
              cluster_rows = F,
              color = c(colorRampPalette(colors = c("green","black"))(length(bk)/2),colorRampPalette(colors = c("black","red"))(length(bk)/2)),
              legend_breaks=seq(-8,8,2),
              breaks=bk,
              fontsize_row = 15,
              fontsize = 15,
              cellwidth = 0.8,
              # silent = T,
              # breaks = c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01)),legend_breaks=seq(-8,8,2)
              border_color = NA,
              annotation_names_col =F,
              main = "STAD(TCGA database)")
  


annotation<-data.frame(Type=sample)
row.names(annotation)<-names(COAD.fpkm_Mreader)[-1]
anno_colors <- list(Type = c("COAD(n=471)" = "red", "Normal(n=41)" = "black"))
COAD.fpkm_M.matrix<-FeaturetoRow(COAD.fpkm_Mreader)
COAD.fpkm_M.scale<-t(scale(t(COAD.fpkm_M.matrix)))

d_labels = "Row Z-Score"
COAD.fpkm_M.scale1<-COAD.fpkm_M.scale[genes,]

bk <- c(seq(-7,-0.1,by=0.01),seq(0,7,by=0.01))

p <- pheatmap(COAD.fpkm_M.scale1,
              # cluster_cols = d$hc,
              # cluster_rows = F,
              # cutree_cols = 3,
              # cutree_rows = 3,
              # gaps_row = c(10, 20),
              annotation_col = annotation,
              annotation_colors = anno_colors,
              show_colnames = F, 
              cluster_cols = F,
              cluster_rows = F,
              color = c(colorRampPalette(colors = c("green","black"))(length(bk)/2),colorRampPalette(colors = c("black","red"))(length(bk)/2)),
              legend_breaks=seq(-8,8,2),
              breaks=bk,
              fontsize_row = 15,
              fontsize = 15,
              cellwidth = 0.7,
              # silent = T,
              # breaks = c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01)),legend_breaks=seq(-8,8,2)
              border_color = NA,
              annotation_names_col =F,
              main = "COAD(TCGA database)")
