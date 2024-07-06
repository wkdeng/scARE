# /**
#  * @author [Wankun Deng]
#  * @email [dengwankun@gmail.com]
#  * @create date 2023-05-09 10:59:27
#  * @modify date 2023-05-25 16:11:31
#  * @desc [description]
#  */

library(Seurat)
library(stringr)
args = commandArgs(trailingOnly = TRUE)
data_folder<-args[1]
out_folder<-args[2]

log_file<-paste0(out_folder,'/','cell_exp.log')
log<-file(log_file,open = 'w')
for(file_name in list.files(data_folder,pattern = 'rds')){
    dataset<-strsplit(file_name,'\\.')[[1]][1]
    print(dataset)
    seurat_obj<-readRDS(paste0(data_folder,'/',file_name))
    exp_mtx<-as.data.frame(t(as.data.frame(GetAssayData(seurat_obj,slot = 'data'))))
    exp_mtx$UMAP_1<-seurat_obj@reductions[['umap']]@cell.embeddings[,1]
    exp_mtx$UMAP_2<-seurat_obj@reductions[['umap']]@cell.embeddings[,2]
    write.table(exp_mtx,paste0(out_folder,'/',str_replace(file_name,'.rds',''),'.cell_exp.txt'),row.names = T,quote = F,sep = '\t',col.names = T)
    print(paste0('Done writing expression file for ',dataset))
    writeLines(c(file_name),log)
}
close(log)