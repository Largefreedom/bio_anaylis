path = 'C:/Users/FREEDOM/Desktop/TCGA_data/mainfest'
floders = list.files(path)#列出总目录下含有文件夹的名称
BRCA_counts = data.frame()#构建dataframe
fd1 = floders[1]#第一个文件；
file_name  =list.files(paste(path,'/',fd1,sep =''))#列出文件夹fd1中的全部文件名
file_list = substr(file_name[1],1,28)#对第一个文件名进行文字截图；
mydata = read.table(gzfile(paste(path,'/',fd1,'/',file_name[1],sep = '')))#读取解压包gz内部的数据；
names(mydata) <- c('ECSG_ID',file_list)#读取第一个gz文件之后，把文件名当作列名；
BRCA_counts = mydata#赋值为BRCA——counts矩阵；
for (fd in floders[2:200]){#循环对200个文件进行处理；
  files_name = list.files(paste(path,'/',fd,sep = ''))
  print(files_name[1])
  file_list =substr(files_name[1],1,28)
  mydata = read.table(gzfile(paste(path,'/',fd,'/',files_name[1],sep = '')))
  names(mydata) <- c('ECSG_ID',file_list)
  BRCA_counts <- merge(BRCA_counts,mydata,by ='ECSG_ID')#进行逐个以ensg编号进行合并；
 
}
write.csv(BRCA_counts,'C:/Users/FREEDOM/Desktop/TCGA_data/BRCA_counts.csv')#把数据写入csv文件中；
group_text =read.csv(file = 'C:/Users/FREEDOM/Desktop/TCGA_data/group_text1.csv',header = T)
library(biomaRt)
library(curl)
#进行基因注释
new_data <- read.csv(file = 'C:/Users/FREEDOM/Desktop/TCGA_data/BRCA_counts1.csv')
rownames(new_data) <- new_data[,1]
new_data <- new_data[c(-1)]
print(rownames(new_data))
# 
char =substr(rownames(new_data),1,15)
print(char)
rownames(new_data) <- substr(rownames(new_data),1,15)
# 
# 
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))  #构建mart
my_ensembl_gene_id<-char #需要转换的exsembl的编码；
my_ensembl_gene_id
# 
# mms_symbols<- getBM(attributes=c('ensembl_gene_id','hgnc_symbol',"description"),filters = 'ensembl_gene_id',values = my_ensembl_gene_id,mart = mart)#基因注释之后的表
# 
# 
ensembl_gen_id <- char
result_diff<-cbind(ensembl_gen_id,new_data)
colnames(new_data)[1]<-c("ensembl_gene_id")
rownames(result_diff) <- NULL
print(colnames(mms_symbols))
print(colnames(result_diff))
colnames(result_diff)[1] <-c("ensembl_gene_id")
colnames(mms_symbols)[1] <- c("ensembl_gene_id")
 
resul_diff <- merge(result_diff,mms_symbols,by = "ensembl_gene_id")  #以ensembl_gene_id为列名进行合并
result_diff <- resul_diff[,1:202]
write.csv(result_diff,'C:/Users/FREEDOM/Desktop/TCGA_data/group_notetext.csv')
#进行差异分析，limma包；
new_data <- read.csv('C:/Users/FREEDOM/Desktop/TCGA_data/group_note.csv',headers <- T)
express_rec <- new_data
rownames(express_rec) <- express_rec[,1]
l
source("https://bioconductor.org/biocLite.R")#载入安装地址；
biocLite("clusterProfiler")#下载程序包；
biocLite("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)#基因注释(ENSG转化为symbol_id用到的两个包；
char <- new_data$ensembl_gene_id

gen_ids <- bitr(char,fromType = 'ENSEMBL',toType = c("SYMBOL", "GENENAME"),OrgDb = 'org.Hs.eg.db')#注释进行处理；
colnames(gen_ids)[1] <-  c("ensembl_gene_id")
gen_ids <- gen_ids[,1:2]
new_data1 <- merge(gen_ids,new_data,by ='ensembl_gene_id')
write.csv(new_data1,'C:/Users/FREEDOM/Desktop/TCGA_data/after_note.csv')

express_rec<- read.csv('C:/Users/FREEDOM/Desktop/TCGA_data/after_note2.csv')#读取数据
group_text <- read.csv('C:/Users/FREEDOM/Desktop/TCGA_data/group_text.csv')

library('DESeq2')#加载包；
install.packages('rpart')#含有这个包可忽略，没有的时候才安装；
express_rec <- express_rec[,-1]
express_rec <- express_rec[,-1]
rownames(express_rec) <-express_rec[,1]
express_rec <- express_rec[(-1)]#表达矩阵的处理；
rownames(group_text) <- group_text[,1]
group_text <- group_text[c(-1)]#分组矩阵的数据处理；
all(rownames(group_text)==colnames(express_rec))#确保表达矩阵的列名与分组矩阵行名相一致；
dds <- DESeqDataSetFromMatrix(countData=express_rec, colData=group_text, design<- ~ group)  #DESeq2的加载
head(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ] #过滤一些low count的数据；
dds <- DESeq(dds)#DESeq进行标准化；
resultsNames(dds)
res <- results(dds)
summary(res)#查看经过标准化矩阵的基本情况；
mcols(res,use.names = TRUE)
res_data <- merge(as.data.frame(res),as.data.frame(counts(dds,normalize = TRUE)),by = 'row.names',sort = FALSE)
res_data <- res_data[,1:7] 
write.csv(res_data,'C:/Users/FREEDOM/Desktop/TCGA_data/result_diff.csv')

#制作plotMA图；
png(file="C:/Users/FREEDOM/Desktop/TCGA_data/plotMA_lfcshrink.png", bg="transparent")
res.shrink <- lfcShrink(dds, contrast = c("group","Tumor","Normal"), res=res)
plotMA(res.shrink, ylim = c(-5,5))#图片边界的限定；
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
dev.off()#图片生成；




res_data <- read.csv('C:/Users/FREEDOM/Desktop/TCGA_data/result_diff.csv',headers <- T)
res_data <- res_data[,-1]#差异分析需要的矩阵提取；
install.packages('ggrepel')#安装程序包；
library(ggplot2)#加载火山图包；
library(ggrepel)
#绘制火山图
rank_data <- res_data[order(res_data[,6]),]#矩阵沿着某一列排序；
volcano_names <- rank_data$Row.names[1:5]#取pvalue最小的五个基因；
res_data$ID2 <- ifelse((res_data$Row.names %in% volcano_names)&abs(res_data$log2FoldChange)>3,gsub('"',''
                ,as.character(res_data$Row.names)),NA)#在矩阵res_data矩阵中创建一个心得列，保存满足|log2folchange|大于3 的基因名，否则保存为NA；
png(file="C:/Users/FREEDOM/Desktop/TCGA_data/myplot_lpval.png", bg="transparent")#先创建一张图片
boundary = ceiling(max(abs(res_data$log2FoldChange)))#确定x轴的边界；
threshold <- ifelse(res_data$pvalue<0.05,ifelse(res_data$log2FoldChange >=3,'UP',ifelse(res_data$log2FoldChange<=(-3),'DW','NoDIFF')),'NoDIFF')#设置分界阀值
ggplot(res_data,aes(x=res_data$log2FoldChange,y =res_data$pvalue,color=threshold))+geom_point(size=1, alpha=0.5) + theme_classic() +
  xlab('log2 fold_change')+ylab(' p-value') +xlim(-1 * boundary, boundary) + theme(legend.position="top", legend.title=element_blank())
+ geom_text_repel(aes(label=res_data$ID2))#火山图中文本标签的注释
dev.off()#保存火山图图片；



deseq2_heatmap <- rank_data[1:30,]#取前30个差异化最大的基因；
write.csv(deseq2_heatmap,'C:/Users/FREEDOM/Desktop/TCGA_data/DEseq2_heatmap.csv')
#绘制热图；


library(pheatmap)#加载热图程序包；


initial_rec <- read.csv('C:/Users/FREEDOM/Desktop/TCGA_data/after_note2.csv')
express_rec <- read.csv('C:/Users/FREEDOM/Desktop/TCGA_data/DEseq2_diffen.csv')
data_rec <-express_rec[1:30,]#取差异性最大的30个基因；
colnames(initial_rec)[2] <- c('gene_name')
data_rec <- merge(initial_rec,data_rec,by <- 'gene_name')
data_rec <-data_rec[,1:202]
data_rec1 <- data_rec[,-2]
group_text <- read.csv('C:/Users/FREEDOM/Desktop/TCGA_data/group_text.csv')
annotation_col = data.frame(#创建annotation_col矩阵，为建立热图的横纵坐标作准备；
  sampleType = factor(group_text$group)
)
data_rec <- data_rec[-16,]
rownames(data_rec) <- data_rec[,1]
data_rec <- data_rec[,-1]
data_rec <- data_rec[,-1]
data_rec[data_rec==0] <-1
data_rec <-log(data_rec,2)#对表达矩阵的值进行标准化；
rownames(annotation_col) <-substr(colnames(data_rec),9,25)#用来截取行名的值；
colnames(data_rec) <-rownames(annotation_col)#确保annotataion_col的行名跟表达矩阵的列名相一致；
#制作热图；
png("C:/Users/FREEDOM/Desktop/TCGA_data/TCGA_GBM_diff.png",height = 800, width = 1600)
pheatmap(data_rec,annotation_col = annotation_col,color <- colorRampPalette(c('red','black','green'))(100),main ='TCGA_BRCA中Tumor和Normal样本差异')
dev.off()

# gene <- deseq2_heatmap['Row.names']#由差异分析结果后提取的基因
# express_rec <- read.csv('C:/Users/FREEDOM/Desktop/TCGA_data/after_note2.csv')
# express_rec <- express_rec[,-1]
# express_rec <- express_rec[,-1]
# colnames(express_rec)[1] <- c('symbol')
# filter_data <- express_rec[gene[,1],]#从表达矩阵筛选出来指定基因的表达情况
# rownames(filter_data) <-  filter_data[,1]
# filter_data <- filter_data[c(-1)]
# filter_data[filter_data==0] <-1#把矩阵中为0的数据转化为1，为后续log做准备；
# filter_data <- log(filter_data,2)#矩阵中的数据进行log2处理；
# data_1 <- as.matrix.data.frame(filter_data)
# data <- matrix(as.numeric(data_1),ncol = 200)#转化为矩阵
# text_group =read.csv('C:/Users/FREEDOM/Desktop/TCGA_data/group_text.csv')#分组矩阵的创建；
# colnames(data) <- as.character(text_group[,1])#分组矩阵列名的设置；
# anno <- data.frame(CellType =factor(text_group[,2]))
# rownames(anno) <- colnames(data)
# anno_colors <- list(CellType = c(Tumor = "#1B9E77", Normal = "#D95F02"))#设置热图lengend(标签)颜色；
# text_sample <- log2(data+1)#数据log化
# rownames(text_sample) <- rownames(data)#热图绘制对象矩阵的行名的设置；
# pheatmap(text_sample,
#          color <- colorRampPalette(c('red','black','green'))(100),cluster_rows = TRUE,
#          cluster_cols = TRUE,
#          main = 'TCGA_BRCA癌症与癌旁基因差异热图',
#          # annotation_col <- anno,
#          )
