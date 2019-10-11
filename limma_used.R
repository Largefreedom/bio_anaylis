library(limma)
# source('https://bioconductor.org/biocLite.R')
# biocLite("edgeR")
library(edgeR)

path = 'C:/Users/FREEDOM/Desktop/TCGA_data/after_note2.csv'
path1 ='C:/Users/FREEDOM/Desktop/TCGA_data/group_text.csv'
express_rec <- read.csv(path,headers <- T)#读取表达矩阵;
group_text <- read.csv(path1,headers <- T)#读取分组矩阵；

rownames(group_text) <- group_text[,1]
group_text <- group_text[c(-1)]

Group <- factor(group_text$group,levels = c('Tumor','Normal'))
design <- model.matrix(~0+Group)
colnames(design) <- c('Tumor','Normal')
rownames(design) <- rownames(group_text)#创建分组矩阵；

express_rec <- express_rec[,-1]
rownames(express_rec) <- express_rec[,1]
express_rec <- express_rec[(-1)]#创建表达矩阵；
express_rec[express_rec == 0] <- 1
express_rec <-log(express_rec,2)

fit <- lmFit(express_rec,design)

contrast.matrix <- makeContrasts(Tumor - Normal,levels=design)
fit2 <- contrasts.fit(fit,contrast.matrix)

fit2 <- eBayes(fit2)

all_diff <- topTable(fit2, adjust.method = 'fdr',coef=1,p.value = 1,lfc <- log(1,2),number = 30000,sort.by = 'logFC')#从高到低排名；

dge <- DGEList(counts = express_rec)
dge <- calcNormFactors(dge)#表达矩阵进行标准化；

v <- voom(dge, design,plot=TRUE)#利用limma_voom方法进行差异分析；
fit <- lmFit(v, design)#线性关系创建；
fit <- eBayes(fit)#贝叶斯算法组建
all <- topTable(fit, coef=ncol(design),n=Inf)#从高到低排名；
sig.limma <- subset(all_diff,abs(all$logFC)>1.5&all$P.Value<0.05)#进行差异基因筛选；
write.csv(sig.limma,'C:/Users/FREEDOM/Desktop/TCGA_data/limm_diff.csv')#写入csv文件中；

write.csv(all,'C:/Users/FREEDOM/Desktop/TCGA_data/limm_rec.csv')


all <- read.csv('C:/Users/FREEDOM/Desktop/TCGA_data/limm_rec.csv')

colnames(all)[1] <- c('name')
#绘制火山图；
library(ggplot2)#加载火山图包；
library(ggrepel)

rank_data <- all[order(all[,5]),]#矩阵沿着pvalue从小到大排序；
rownames(rank_data) <-rank_data[,1]
rank_data <- rank_data[c(-1)]
rank_data$names <- rownames(rank_data)
volcano_names <- rownames(rank_data)[1:5]#取pvalue最小的五个基因；

rank_data$ID2 <- ifelse((rank_data$names %in% volcano_names)&abs(rank_data$logFC)>3
    ,as.character(rank_data$names),NA)#在矩阵res_data矩阵中创建一个心得列，保存满足|log2folchange|大于3 的基因名，否则保存为NA；
png(file="C:/Users/FREEDOM/Desktop/TCGA_data/limma_voloun_log1.png", bg="transparent")#先创建一张图片
boundary = ceiling(max(abs(rank_data$logFC)))#确定x轴的边界；
threshold <- ifelse(rank_data$P.Value<0.05,ifelse(rank_data$logFC >=3,'UP',ifelse(rank_data$logFC<=(-3),'DW','NoDIFF')),'NoDIFF')#设置分界阀值
ggplot(rank_data,aes(x=rank_data$logFC,y =(-1)*log10(rank_data$P.Value),color=threshold),abline(v=c(-log(1.5,2),log(1.5,2))),h =-log10(0.05))+geom_point(size=1, alpha=0.5) + theme_classic() +
  xlab('log2 fold_change')+ylab(' -log10 p-value') +xlim(-1 * boundary, boundary) + theme(legend.position="top", legend.title=element_blank())  + geom_text_repel(aes(label=rank_data$ID2))#火山图中文本标签的注释
dev.off()#保存火山图图片；







