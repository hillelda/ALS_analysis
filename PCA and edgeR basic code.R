
          # Libraries etc. #####
library(edgeR) ; library(ggplot2) ; library(ggfortify) ; library(matrixStats)

setwd('./working data - ALS')
meta_trfs<-read.csv('GSE203170_tRF_meta.csv') ; meta_trfs<-meta_trfs[,-1]


#### tRFs ####
CountsDataFrame<-read.csv('GSE203170_tRNA_Exclusive_Combined_data.csv')
rownames(CountsDataFrame)<-CountsDataFrame$X ; CountsDataFrame<-CountsDataFrame[,-1]
ColData<-read.csv('GSE203170_Coldata.csv')
ColData$group<-factor(ColData$group,levels=c('healthy','ALS'))

tmp1<-CountsDataFrame[rowMedians(as.matrix(CountsDataFrame))>1,] ; pcaData1<-as.data.frame(t(tmp1))
pcaData1$id<-rownames(pcaData1) ; pcaData1<-merge(pcaData1,ColData,by='id')
rownames(pcaData1)<-pcaData1$id ; pca_res1<-prcomp(pcaData1[,2:nrow(tmp1)],scale. = T) #prcomp does the PCA

c='genetics' # change c to the parameter you want to lable the PCA with
autoplot(pca_res1, x=1,y=2 , data = pcaData1, colour = c,label=T,size=0.3)
autoplot(pca_res1, x=3,y=4 , data = pcaData1, colour = c,label=T,size=0.3)
autoplot(pca_res1, x=5,y=6 , data = pcaData1, colour = c,label=T,size=0.3)
autoplot(pca_res1, x=7,y=8 , data = pcaData1, colour = c,label=T,size=0.3)
autoplot(pca_res1, x=9,y=10 , data = pcaData1, colour = c,label=T,size=0.3)


cts<-CountsDataFrame ; cold1<-ColData
cold1<-subset(cold1,cold1$id %in% colnames(cts)) ; cts<-cts[,as.character(cold1$id)]
cts<-cts[,order(cold1$id)] ; cold1<-cold1[order(cold1$id),] ; nrow(cold1)==sum(cold1$id==colnames(cts))

# edgeR from here
y <- DGEList(counts=cts,group=cold1$group) # change the group to the name of your group column
keep <- filterByExpr(y) ; y <- y[keep, , keep.lib.sizes=FALSE] ; y$samples$lib.size <- colSums(y$counts) ; y <- calcNormFactors(y)

cts1<-as.data.frame(cpm(y,log = F)) # creates the normalized counts matrix

dsgn <- model.matrix(~group, data = cold1)
y <- estimateDisp(y, dsgn, robust = T)

## if you want, you can run the next code to see if the different coeficients correlate with eachother:
# logFC <- predFC(y,dsgn,prior.count=1,dispersion=0.05) ; cor(logFC) ; plotBCV(y)

head(dsgn)
# you change the coef to the coeficient that intersts you (the number of the column in the dsgn matrix)
fit <- glmQLFit(y, dsgn, robust = T) ; lrt <- glmLRT(fit,coef = 2)

sgGens<-as.data.frame(topTags(lrt,adjust.method = 'fdr',n = nrow(cts1))) ; sgGens$transcript<-rownames(sgGens)


write.csv(sgGens, file='GSE203170_tRF_SgGens_output.csv')


############################
############################
# MiRNAs ####


CountsDataFrame2<-read.csv('GSE203170_unnormcounts_compiled.csv')
rownames(CountsDataFrame2)<-CountsDataFrame2$X.miRNA ; CountsDataFrame2<-CountsDataFrame2[,-1]
CountsDataFrame2<-CountsDataFrame2[,-1]; CountsDataFrame2<-CountsDataFrame2[,-1]
ColData<-read.csv('GSE203170_Coldata.csv')
ColData$group<-factor(ColData$group,levels=c('healthy','ALS'))

tmp2<-CountsDataFrame2[rowMedians(as.matrix(CountsDataFrame2))>1,] ; pcaData2<-as.data.frame(t(tmp2))
pcaData2$id<-rownames(pcaData2) ; pcaData2<-merge(pcaData2,ColData,by='id')
rownames(pcaData2)<-pcaData2$id ; pca_res2<-prcomp(pcaData2[,2:nrow(tmp2)],scale. = T) #prcomp does the PCA

c='genetics' # change c to the parameter you want to lable the PCA with
autoplot(pca_res2, x=1,y=2 , data = pcaData2, colour = c,label=T,size=0.3)
autoplot(pca_res2, x=3,y=4 , data = pcaData2, colour = c,label=T,size=0.3)
autoplot(pca_res2, x=5,y=6 , data = pcaData2, colour = c,label=T,size=0.3)
autoplot(pca_res2, x=7,y=8 , data = pcaData2, colour = c,label=T,size=0.3)
autoplot(pca_res2, x=9,y=10 , data = pcaData2, colour = c,label=T,size=0.3)


cts2<-CountsDataFrame2 ; cold2<-ColData
cold2<-subset(cold2,cold2$id %in% colnames(cts)) ; cts2<-cts2[,as.character(cold2$id)] 
cts2<-cts2[,order(cold2$id)] ; cold2<-cold2[order(cold2$id),] ; nrow(cold2)==sum(cold2$id==colnames(cts2))

# edgeR from here
y <- DGEList(counts=cts,group=cold2$group) # change the group to the name of your group column
keep <- filterByExpr(y) ; y <- y[keep, , keep.lib.sizes=FALSE] ; y$samples$lib.size <- colSums(y$counts) ; y <- calcNormFactors(y)

cts2<-as.data.frame(cpm(y,log = F)) # creates the normalized counts matrix

dsgn <- model.matrix(~group, data = cold2)
y <- estimateDisp(y, dsgn, robust = T)

## if you want, you can run the next code to see if the different coeficients correlate with eachother:
# logFC <- predFC(y,dsgn,prior.count=1,dispersion=0.05) ; cor(logFC) ; plotBCV(y)

head(dsgn)
# you change the coef to the coeficient that intersts you (the number of the column in the dsgn matrix)
fit <- glmQLFit(y, dsgn, robust = T) ; lrt <- glmLRT(fit,coef = 2) 

sgGens2<-as.data.frame(topTags(lrt,adjust.method = 'fdr',n = nrow(cts2))) ; sgGens2$transcript<-rownames(sgGens2)


write.csv(sgGens2, file='GSE203170_miRNA_SgGens_output.csv')


########################################
        #second data base#
#######################################
#### tRFs ####
CountsDataFrame<-read.csv('GSE94888_tRNA_Exclusive_Combined_data.csv')
rownames(CountsDataFrame)<-CountsDataFrame$X ; CountsDataFrame<-CountsDataFrame[,-1]
ColData2<-read.csv('GSE94888_Coldata.csv')
ColData2$group<-factor(ColData2$group,levels=c('healthy','ALS'))

tmp3<-CountsDataFrame[rowMedians(as.matrix(CountsDataFrame))>1,] ; pcaData3<-as.data.frame(t(tmp3))
pcaData3$id<-rownames(pcaData3) ; pcaData3<-merge(pcaData3,ColData2,by='id')
rownames(pcaData3)<-pcaData3$id ; pca_res3<-prcomp(pcaData3[,2:nrow(tmp3)],scale. = T) #prcomp does the PCA

c='genetics' # change c to the parameter you want to lable the PCA with
autoplot(pca_res3, x=1,y=2 , data = pcaData3, colour = c,label=T,size=0.3)
autoplot(pca_res3, x=3,y=4 , data = pcaData3, colour = c,label=T,size=0.3)
autoplot(pca_res3, x=5,y=6 , data = pcaData3, colour = c,label=T,size=0.3)
autoplot(pca_res3, x=7,y=8 , data = pcaData3, colour = c,label=T,size=0.3)
autoplot(pca_res3, x=9,y=10 , data = pcaData3, colour = c,label=T,size=0.3)


cts3<-CountsDataFrame ; cold3<-ColData2
cold3<-subset(cold3,cold3$id %in% colnames(cts)) ; cts3<-cts3[,as.character(cold3$id)]
cts3<-cts3[,order(cold3$id)] ; cold3<-cold3[order(cold3$id),] ; nrow(cold3)==sum(cold3$id==colnames(cts3))

# edgeR from here
y <- DGEList(counts=cts3,group=cold3$group) # change the group to the name of your group column
keep <- filterByExpr(y) ; y <- y[keep, , keep.lib.sizes=FALSE] ; y$samples$lib.size <- colSums(y$counts) ; y <- calcNormFactors(y)

cts3<-as.data.frame(cpm(y,log = F)) # creates the normalized counts matrix

dsgn <- model.matrix(~group, data = cold3)
y <- estimateDisp(y, dsgn, robust = T)

## if you want, you can run the next code to see if the different coeficients correlate with eachother:
# logFC <- predFC(y,dsgn,prior.count=1,dispersion=0.05) ; cor(logFC) ; plotBCV(y)

head(dsgn)
# you change the coef to the coeficient that intersts you (the number of the column in the dsgn matrix)
fit <- glmQLFit(y, dsgn, robust = T) ; lrt <- glmLRT(fit,coef = 2)

sgGens3<-as.data.frame(topTags(lrt,adjust.method = 'fdr',n = nrow(cts3))) ; sgGens3$transcript<-rownames(sgGens3)


write.csv(sgGens, file='GSE94888_tRF_SgGens_output.csv')


############################
############################
# MiRNAs


CountsDataFrame4<-read.csv('GSE94888_unnormcounts_compiled.csv')
rownames(CountsDataFrame4)<-CountsDataFrame4$X.miRNA ; CountsDataFrame4<-CountsDataFrame4[,-1]
CountsDataFrame4<-CountsDataFrame2[,-1]; CountsDataFrame2<-CountsDataFrame2[,-1]
ColData2<-read.csv('GSE203170_Coldata.csv')
ColData2$group<-factor(ColData2$group,levels=c('healthy','ALS'))

tmp4<-CountsDataFrame4[rowMedians(as.matrix(CountsDataFrame4))>1,] ; pcaData4<-as.data.frame(t(tmp4))
pcaData4$id<-rownames(pcaData4) ; pcaData4<-merge(pcaData4,ColData2,by='id')
rownames(pcaData4)<-pcaData4$id ; pca_res4<-prcomp(pcaData4[,2:nrow(tmp4)],scale. = T) #prcomp does the PCA

c='genetics' # change c to the parameter you want to lable the PCA with
autoplot(pca_res4, x=1,y=2 , data = pcaData4, colour = c,label=T,size=0.3)
autoplot(pca_res4, x=3,y=4 , data = pcaData4, colour = c,label=T,size=0.3)
autoplot(pca_res4, x=5,y=6 , data = pcaData4, colour = c,label=T,size=0.3)
autoplot(pca_res4, x=7,y=8 , data = pcaData4, colour = c,label=T,size=0.3)
autoplot(pca_res4, x=9,y=10 , data = pcaData4, colour = c,label=T,size=0.3)


cts4<-CountsDataFrame4 ; cold4<-ColData2
cold4<-subset(cold4,cold4$id %in% colnames(cts)) ; cts4<-cts4[,as.character(cold4$id)] 
cts4<-cts4[,order(cold4$id)] ; cold4<-cold4[order(cold4$id),] ; nrow(cold4)==sum(cold4$id==colnames(cts4))

# edgeR from here
y <- DGEList(counts=cts,group=cold4$group) # change the group to the name of your group column
keep <- filterByExpr(y) ; y <- y[keep, , keep.lib.sizes=FALSE] ; y$samples$lib.size <- colSums(y$counts) ; y <- calcNormFactors(y)

cts4<-as.data.frame(cpm(y,log = F)) # creates the normalized counts matrix

dsgn <- model.matrix(~group, data = cold4)
y <- estimateDisp(y, dsgn, robust = T)

## if you want, you can run the next code to see if the different coeficients correlate with eachother:
# logFC <- predFC(y,dsgn,prior.count=1,dispersion=0.05) ; cor(logFC) ; plotBCV(y)

head(dsgn)
# you change the coef to the coeficient that intersts you (the number of the column in the dsgn matrix)
fit <- glmQLFit(y, dsgn, robust = T) ; lrt <- glmLRT(fit,coef = 2) 

sgGens4<-as.data.frame(topTags(lrt,adjust.method = 'fdr',n = nrow(cts4))) ; sgGens4$transcript<-rownames(sgGens4)


write.csv(sgGens4, file='GSE94888_miRNA_SgGens_output.csv')
