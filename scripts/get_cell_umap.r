# /**
#  * @author [Wankun Deng]
#  * @email [dengwankun@gmail.com]
#  * @create date 2023-05-09 10:06:03
#  * @modify date 2023-10-08 22:35:05
#  * @desc [description]
#  */

library(Seurat)
library(stringr)

args = commandArgs(trailingOnly = TRUE)
data_folder<-args[1]
out_folder<-args[2]
## creat log file of finishing this step
log_file<-paste0(out_folder,'/','cell_umap.log')
log<-file(log_file,open = 'w')

cell_type_map<-list('Excitatory'='Ex','Inhibitory'='In','OPC'='OPC','Oligodendrocyte'='Oli','Microglia'='Mic','Endothelial'='Endo','Astrocyte'='Ast','VLMC/Per'='VLMC')
for(file_name in list.files(data_folder,pattern = 'rds')){
    cell_meta<-NULL
    dataset<-strsplit(file_name,'\\.')[[1]][1]
    print(dataset)
    seurat_obj<-readRDS(paste0(data_folder,'/',file_name))
    seurat_obj@meta.data$UMAP_1<-seurat_obj@reductions[['umap']]@cell.embeddings[,1]
    seurat_obj@meta.data$UMAP_2<-seurat_obj@reductions[['umap']]@cell.embeddings[,2]
    seurat_obj@meta.data$dataset<-dataset
    seurat_obj@meta.data$stage<-'N/A'
    # seurat_obj@meta.data[seurat_obj@meta.data$predicted.celltype=='Opc','predicted.celltype']<-'OPC'
    if (!('msex' %in% colnames(seurat_obj@meta.data))){
        seurat_obj@meta.data$msex<-'N/A'
    }
    if (!('age_death' %in% colnames(seurat_obj@meta.data))){
        seurat_obj@meta.data$age_death<-'N/A'
    }
    if(is.null(cell_meta)){
        cell_meta<-seurat_obj@meta.data[,c('major.cell.type.prediction.seaad','Diagnosis','stage','msex','age_death','UMAP_1','UMAP_2','dataset')]
    }else{
        cell_meta<-rbind(cell_meta,seurat_obj@meta.data[,c('major.cell.type.prediction.seaad','Diagnosis','stage','msex','age_death','UMAP_1','UMAP_2','dataset')])
    }
    cell_meta["Stage" %in% cell_meta$Diagnosis,'Diagnosis']<-'AD'
    cell_meta["control" %in% cell_meta$Diagnosis,'Diagnosis']<-'Control'
    # cell_meta[cell_meta$predicted.celltype=='Opc','predicted.celltype']<-'OPC'
    cell_meta$predicted.celltype<-unlist(lapply(cell_meta$major.cell.type.prediction.seaad,function(x) cell_type_map[[x]]))
    write.table(cell_meta,paste0(out_folder,'/',str_replace(file_name,'.rds',''),'.cell_umap.txt'),row.names = ,quote = F,sep = '\t',col.names = T)
    writeLines(c(file_name),log)
}
close(log)

