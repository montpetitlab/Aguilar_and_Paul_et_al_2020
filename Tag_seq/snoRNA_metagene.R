library(reshape2)
library(ggplot2)
m = read.delim("~/Paul_et_al_2019/Tag_seq/matrix2.mat", skip=2, header=F)
m = as.matrix(m[,-c(1:6)])
wt_pBM5 = m[,c(1:40)]
wt_pBM5=colMeans(wt_pBM5)
enp1_1_pBM5=m[,c(41:80)]
enp1_1_pBM5=colMeans(enp1_1_pBM5)
csl4ph_pBM5=m[,c(81:120)]
csl4ph_pBM5=colMeans(csl4ph_pBM5)
WT_pBM766=m[,c(121:160)]
WT_pBM766=colMeans(WT_pBM766)
enp1_1_pBM766=m[,c(161:200)]
enp1_1_pBM766=colMeans(enp1_1_pBM766)
csl4ph_pBM766=m[,c(201:240)]
csl4ph_pBM766=colMeans(csl4ph_pBM766)
df=data.frame("wt_pBM5"=wt_pBM5,"enp1_1_pBM5"=enp1_1_pBM5,"csl4ph_pBM5"=csl4ph_pBM5,"wt_pBM766"=WT_pBM766,
              "enp1_1_pBM766"=enp1_1_pBM766,"csl4ph_pBM766"=csl4ph_pBM766, "bin"=1:40)

df1 <- melt(df, id=c("bin"))
pBM5 <- df1[grep("pBM5", df1$variable),]
pBM5$plasmid <- "pBM5"
pBM766 <- df1[grep("pBM766", df1$variable),]
pBM766$plasmid <- "pBM766"
df2 <- rbind(pBM5, pBM766)
wt <- df2[grep("wt_pBM", df2$variable),]
wt$strain <- "wt"
enp1 <- df2[grep("enp1_1_pBM", df2$variable),]
enp1$strain <- "enp1-1"
csl4 <- df2[grep("csl4", df2$variable),]
csl4$strain <- "csl4-ph"
df3 <- rbind(wt,enp1,csl4)
ggplot(df3, aes(bin, value, group = interaction(strain, plasmid), 
               color = strain, linetype = plasmid)) +
  geom_line(size=1.2) + xlab("Gene body") + ylab("Number of Reads") + theme_classic()

