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