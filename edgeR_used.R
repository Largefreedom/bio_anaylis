path = 'C:/Users/FREEDOM/Desktop/TCGA_data/after_note2.csv'
path1 ='C:/Users/FREEDOM/Desktop/TCGA_data/group_text.csv'
express_rec <- read.csv(path,headers <- T)#读取表达矩阵;
group_text <- read.csv(path1,headers <- T)#读取分组矩阵；
library(edgeR)#加载edgeR包

express_rec <- express_rec[,-1]
rownames(express_rec) <- express_rec[,1]
express_rec <- express_rec[(-1)]#创建表达矩阵；
rownames(group_text) <- group_text[,1]
group_text <- group_text[c(-1)]#加载分组矩阵；

group <-factor(group_text$group)

dge <- DGEList(counts = express_rec,group = group)#构建DEList对象；
y <- calcNormFactors(dge)#利用calcNormFactor函数对DEList对象进行标准化(TMM算法) 

#创建设计矩阵，跟Limma包相似；
rownames(group_text) <- group_text[,1]
group_text <- group_text[c(-1)]
Group <- factor(group_text$group,levels = c('Tumor','Normal'))
design <- model.matrix(~0+Group)
colnames(design) <- c('Tumor','Normal')
rownames(design) <- rownames(group_text)#创建分组矩阵；

y <- estimateDisp(y,design)#估计离散值（Dispersion）
plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE)#进一步通过quasi-likelihood (QL)拟合NB模型
head(fit$coefficients)

TU.vs.NO <- makeContrasts(Tumor-Normal, levels=design)#这一步主要构建比较矩阵；
res <- glmQLFTest(fit, contrast=TU.vs.NO)#用QL F-test进行检验
# ig.edger <- res$table[p.adjust(res$table$PValue, method = "BH") < 0.01, ]#利用‘BH’方法；
result_diff <- res$table#取出最终的差异基因；
write.csv(edge_diff,'C:/Users/FREEDOM/Desktop/TCGA_data/edgeR_diff2.csv')
edge_diff <- subset(result_diff,abs(result_diff$logFC)>1.5&result_diff$PValue<0.05)


#绘制火山图;
library(ggplot2)#加载火山图包；
library(ggrepel)

gene_name <- data.frame(rownames(result_diff))
rank_data <- cbind(gene_name,result_diff)
colnames(rank_data)[1] <- c('gene_name')

rank_data <- rank_data[order(rank_data[,5]),]#矩阵沿着pvalue从小到大排序；
rownames(rank_data) <-rank_data[,1]
rank_data <- rank_data[c(-1)]
rank_data$names <- rownames(rank_data)
volcano_names <- rownames(rank_data)[1:5]#取pvalue最小的五个基因；

rank_data$ID2 <- ifelse((rank_data$names %in% volcano_names)&abs(rank_data$logFC)>3
                        ,as.character(rank_data$names),NA)#在矩阵res_data矩阵中创建一个心得列，保存满足|log2folchange|大于3 的基因名，否则保存为NA；
png(file="C:/Users/FREEDOM/Desktop/TCGA_data/edege_voloun_log2.png", bg="transparent")#先创建一张图片
boundary = ceiling(max(abs(rank_data$logFC)))#确定x轴的边界；
threshold <- ifelse(rank_data$PValue<0.05,ifelse(rank_data$logFC >=3,'UP',ifelse(rank_data$logFC<=(-3),'DW','NoDIFF')),'NoDIFF')#设置分界阀值
ggplot(rank_data,aes(x=rank_data$logFC,y =rank_data$PValue,color=threshold))+geom_point(size=1, alpha=0.5) + theme_classic() +
  xlab('log2 fold_change')+ylab(' log2 p-value') +xlim(-1 * boundary, boundary) + theme(legend.position="top", legend.title=element_blank()) + geom_text_repel(aes(label=rank_data$ID2))#火山图中文本标签的注释
dev.off()#保存火山图图片；







