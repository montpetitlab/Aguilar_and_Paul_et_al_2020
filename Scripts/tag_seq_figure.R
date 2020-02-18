csl4ph_pBM766 <- read.csv("~/Paul_et_al_2019/Tag_seq/WT_pBM5vscsl4ph_pBM766.tsv", sep = "\t")[,c(6,1)]
enp1_pBM766 <- read.csv("~/Paul_et_al_2019/Tag_seq/WT_pBM5vsenp1_1_pBM766.tsv", sep="\t")[,c(6,1)]
wt_pBM766 <- read.csv("~/Paul_et_al_2019/Tag_seq/WT_pBM5vsWT_pBM766.tsv", sep="\t")[,c(6,1)]
enp1_1_pBM05 <- read.csv("~/Paul_et_al_2019/Tag_seq/WT_pBM5vsenp1_1_pBM5.tsv", sep="\t")[,c(6,1)]
csl4_ph_pBM05 <- read.csv("~/Paul_et_al_2019/Tag_seq/WT_pBM5vscsl4ph_pBM5.tsv", sep="\t")[,c(6,1)]

temp <- merge(wt_pBM766, csl4_ph_pBM05, by="mrna", all = T)
temp <- merge(temp, csl4ph_pBM766, by="mrna", all = T)
temp <- merge(temp, enp1_1_pBM05, by="mrna", all = T)
temp <- merge(temp, enp1_pBM766, by="mrna", all = T)

#temp <- merge(temp, csl4ph_pBM5vscsl4ph_pBM766, by="mrna", all = T)
#temp <- merge(temp, enp1_1_pBM5vsenp1_1_pBM766, by="mrna", all = T)
colnames(temp) <- c("RNA", "WT_pBM766", "csl4_pBM05", "csl4ph_pBM766", "enp1_1_pBM05", "enp1_1_pBM766")
gene_con <- read.csv("/Users/montpetitlab1/Paul_et_al_2019/Tag_seq/gene_conversion_1.csv", sep ="\t")
colnames(gene_con) <- c("RNA", "Gene")
temp <- merge(gene_con, temp, by="RNA", all = F)
rna_class_trv <- read.csv("~/Paul_et_al_2019/Tag_seq/rna_class2.csv", sep = "\t")
colnames(rna_class_trv) <- c("Gene", "RNA_Class")
temp <- merge(temp, rna_class_trv, by="Gene", all = F)
temp <- temp[,c()]
#temp <- temp[complete.cases(temp),]
df <- melt(temp)
df <- df[complete.cases(df),]
levels(df$RNA_Class)[levels(df$RNA_Class) %in%  c("I","II","III", 'IV')] <- "I-III"
levels(df$RNA_Class)[levels(df$RNA_Class) %in%  c("V","VI","VII", 'VIII', 'IX')] <- "IV-IX"
p <- ggplot(df, aes(variable, value, fill=variable)) + geom_boxplot(width=0.7, outlier.size = 0.1, notch=TRUE) + facet_wrap(~RNA_Class) + theme_classic() + labs(title = "mRNA Class", x="Strains", y="log2FC") 
p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12)) + geom_point(position = "jitter", size=0.1)

rna_class <- read.csv("~/Paul_et_al_2019/Tag_seq/mRNA_list.txt", sep = "\t", header = F)[,c(1,3)]
colnames(rna_class) <- c("Gene", "RNA_Class")
temp <- merge(temp, rna_class, by="RNA", all=F)
CUT <- temp[grep("CUT", temp$RNA),]
df <- melt(CUT)
ggplot(df, aes(variable, value, fill=variable)) + geom_violin(width=1.3) + geom_boxplot(width=0.1) + 
  theme_classic() + labs(title = "CUTs", x="Strains", y="log2FC") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12), axis.text.y = element_text(size=12))

CUT <- temp[grep("SUT", temp$RNA),]
df <- melt(CUT)
ggplot(df, aes(variable, value, fill=variable)) + geom_violin(width=1.3) + geom_boxplot(width=0.1) + 
  theme_classic() +  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12), axis.text.y = element_text(size=12)) +
  labs(title = "SUTs", x="Strains", y="log2FC") + geom_jitter(size=0.1)
CUT <- temp[grep("XUT", temp$RNA),]
df <- melt(CUT)
ggplot(df, aes(variable, value, fill=variable)) + geom_violin(width=1.3) + geom_boxplot(width=0.1) + 
   theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12), axis.text.y = element_text(size=12)) +
  labs(title = "XUTs", x="Strains", y="log2FC") + geom_jitter(size=0.1)
CUT <- temp[grep("NUT", temp$RNA),]
df <- melt(CUT)
ggplot(df, aes(variable, value, fill=variable)) + geom_violin(width=1.3) + geom_boxplot(width=0.1) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12), axis.text.y = element_text(size=12)) + 
  labs(title = "NUTs", x="Strains", y="log2FC") + geom_jitter(size=0.1)

temp <- merge(enp1_1_pBM05, enp1_pBM766, by="mrna", all = T)
colnames(temp) <- c("RNA", "enp1_1_pBM05", "enp1_1_pBM766")
rna_class <- read.csv("~/Desktop/rna3.txt", sep="\t")
colnames(rna_class) <- c("RNA", "RNA_Class")
temp <- merge(temp, rna_class, by="RNA", all = T)
XUT <- temp[grep("XUT", temp$RNA),]
XUT$RNA_Class = "XUT"
temp <- rbind(temp, XUT)
temp <- temp %>% filter(RNA_Class!="mRNA")
df <- melt(temp)
df <- df[complete.cases(df),]
ggplot(df, aes(x=RNA_Class, y=value, fill=variable)) + geom_violin() + ylim(-5,5) +theme(axis.text.x = element_text(angle = 45)) + theme_classic() 
df <- melt(temp)
df <- df[complete.cases(df),]

temp <- merge(enp1_1_pBM05, enp1_pBM766, by="mrna", all = T)
colnames(temp) <- c("RNA", "enp1_1_pBM05", "enp1_1_pBM766")
rna_class_trv <- read.csv("~/Downloads/rna_class2.csv", sep = "\t")
colnames(rna_class_trv) <- c("RNA", "RNA_Class")
temp <- merge(temp, rna_class_trv, by="RNA", all = F)
df <- melt(temp)
df <- df[complete.cases(df),]
levels(df$RNA_Class)[levels(df$RNA_Class) %in%  c("I","II","III", 'IV')] <- "I-IV"
levels(df$RNA_Class)[levels(df$RNA_Class) %in%  c("V","VI","VII", 'VIII', 'IX')] <- "V-IX"
ggplot(df, aes(x=RNA_Class, y=value, fill=variable)) + geom_violin() + ylim(-5,5)+ theme_classic() + 
