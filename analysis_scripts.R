#Created by Tomas Engelthaler & Thomas Hills.

#This file contains the analyses reported in our manuscript.
#If you have any questions, please contact us on t.engelthaler@warwick.ac.uk

#The file uses three external datasets (not produced by us), these are:
#1. The McRae feature norms ("norms"), retrieved from https://sites.google.com/site/kenmcraelab/norms-data
#2. The Kuperman age of acquisition norms ("aoa"), retrieved from http://crr.ugent.be/archives/806
#3. The MCDI age of acquisition norms ("words"), retrieved from http://mb-cdi.stanford.edu/

#Important results are stored in dataframes named "result1, result2, result3 ... result11"

#Set the working directory to where all the external datasets are.
#Run the code from the beginning to the end. Some analyses rely on running previous bits of code first. It is important not to skip anything.
#Runing most analyses shouldn't take long. Network analyses (result 6 to 8) take about 5 minutes each. Clustering analysis takes 5-10 minutes.

#Install relevant packages
install.packages('gridExtra')
install.packages('dplyr')
install.packages('caret')
install.packages('skmeans')
install.packages('cluster')
install.packages('QuantPsyc')
install.packages('RVAideMemoire')
install.packages('xlsx')
install.packages('mada')

#Load packages
library('gridExtra')
library('dplyr')
library('caret')
library('skmeans')
library('cluster')
library('QuantPsyc')
library('RVAideMemoire')
library('xlsx')
library('mada')

#Load raw data
norms <- read.csv("mcrae.csv", header = T)
words <- read.csv("ToddlerNorms.csv", header = T)
aoab <- read.table('AoA_ratings_Kuperman_et_al_BRM.txt', header = T)

#Load pre-processed feature norms (to save time, can be generated from the McRae feature norms)
load('wordconc')

#Exclude encyclopaedic features.
norms <- subset(norms, norms$BR_Label != "encyclopaedic")
for(word in length(colnames(wordconc)) : 1){
  if((colnames(wordconc)[word] %in% unique(norms$Feature)) == FALSE){
    print(paste("Removing: " , "(" , word , ")", colnames(wordconc)[word]))
    wordconc <- wordconc[,-word]
  }
}

#List of words shared between AoA norms and feature norms
sharedwords <- unique((inner_join(aoab, norms, by = c("Word" = "Concept")))$Word)

#Proportion of visual form and surface features
nrow(subset(norms, norms$BR_Label == "visual-form_and_surface")) / nrow(norms)

#Process MCDI AoA to show first month where over 50% of children learned the word (newaoa)
ww <- words[order(words$Word),]
wwl <- subset(ww, ww$M30 >= 50)
newaoa <- data.frame(matrix(0, nrow = nrow(wwl), ncol = 1))
done <- 0

for (wordn in 1:(nrow(wwl)))
{
  done <- 0
  for (mcol in 3: (ncol(wwl)))
  {
    if (done == 0)
    {
      if (wwl[wordn,mcol] >= 50)
      {
        done <- 1;
        newaoa[wordn,1] <- mcol + 13;
      }
    }
  }
}
newaoa <- cbind(wwl$Word,newaoa)
names(newaoa) <- c("word","aoa")

#Process McRae feature norms to show sum distinctiveness for each word (dist)
dist<-with(norms, tapply(Distinct, list(Concept), sum))
dist <- data.frame(word=rownames(dist), dist)

##################################################################
####ANALYSES A: Overall distinctiveness and age of acquisition####
##################################################################

#Correlate word distinctiveness (dist) against MCDI AoA (newaoa)
test <- inner_join(newaoa, dist)

  ##Normality tests
shapiro.test(test$aoa)
shapiro.test(test$dist)
qqnorm(test$aoa, main = "MCDI Q-Q Plot")
qqline(test$aoa)

density <- density(test$aoa)
plot(density , main = "MCDI Density" , xlab = "AoA(months)")
text(30, y = 0.12, label = "Shapiro-Wilk p = 0.01")

result1 <- with(test, cor.test(aoa, dist, method="spearman")) ##RESULT 1 - more distinctive words learned earlier

#Correlate word distinctiveness (dist) against Brysbaert's AoA (aoab)
names(aoab)<-c('word',names(aoab[,-1]))
dt <- inner_join(dist, aoab)

  ##Normality tests
shapiro.test(dt$dist)

qqnorm(as.numeric(dt$dist), main = "Disctinctiveness Q-Q Plot")
qqline(as.numeric(dt$dist))

density <- density(dt$dist)
plot(density , main = "Distinctiveness Density" , xlab = "Distinctiveness")
text(12, y = 0.12, label = "Shapiro-Wilk p < 0.001")

shapiro.test(dt$Rating.Mean)
result2 <- with(dt, cor.test(dist, Rating.Mean, method="spearman")) ##RESULT 2 - more distinctive words learned earlier (with larger AoA dataset)

#Correlate MCDI AoA (newaoa) against Brysbaert's AoA (aoab) 
dt$dist <- as.numeric(dt$dist)
dtt <- inner_join(newaoa, dt)
result3 <- with(dtt, cor.test(aoa,Rating.Mean, method="spearman")) ##RESULT 3 - cor between brysbaert and mcdi = .58

#############################################
###Analyses B: Distinctiveness regressions###
#############################################

# regressions based on feature types (distinctness as count - 1 or 0)
dist2<-with(norms, tapply(Distinct==1, list(Concept, BR_Label), sum, na.rm=T))
for(i in 1:nrow(dist2)){
  for(j in 1:ncol(dist2)){
    if(is.na(dist2[i,j])) dist2[i,j]<-0
  }
}

dist2 <- data.frame(word=rownames(dist2), dist2)
names(aoab)<-c('word',names(aoab[,-1]))
dt2 <- inner_join(dist2, aoab)
dtx <- dt2[,c(2:11)]
dty <- dt2[,15]
set.seed(100)
indx <- createFolds(dty, returnTrain = TRUE)
ctrl <- trainControl(method = "cv", index = indx)
set.seed(100)
lmTune0 <- train(x = dtx, y = dty,
                 method = "lm",
                 trControl = ctrl,
                 preProc = c("center","scale"))

result4 <- summary(lmTune0) ##RESULT 4 - visual form and surface the most predictive feature against Brysbaert's AoA (aoab)

  #Alternative regression analysis
dt2.z <- Make.Z(dt2[,-1])
dt2.z <- as.data.frame(dt2.z)
result4alt <- lm(formula = dt2.z$Rating.Mean ~ dt2.z$function. + dt2.z$smell + dt2.z$sound + dt2.z$tactile + dt2.z$taste + dt2.z$taxonomic + dt2.z$visual.colour + dt2.z$visual.form_and_surface + dt2.z$visual.motion)

lm.beta(result4alt)
confint(result4alt)

# Regression based on feature types (distinctness as a scale)

dist3<-with(norms, tapply(Distinct, list(Concept, BR_Label), sum, na.rm=T))
for(i in 1:nrow(dist3)){
  for(j in 1:ncol(dist3)){
    if(is.na(dist3[i,j])) dist3[i,j]<-0
  }
}

dist3 <- data.frame(word=rownames(dist3), dist3)
names(aoab)<-c('word',names(aoab[,-1]))
dt3 <- inner_join(dist3, aoab)
dtx2 <- dt3[,c(2:11)]
dty2 <- dt3[,'Rating.Mean']
set.seed(100)
indx2 <- createFolds(dty2, returnTrain = TRUE)
ctrl2 <- trainControl(method = "cv", index = indx2)
set.seed(100)
lmTune1 <- train(x = dtx2, y = dty2,
                 method = "lm",
                 trControl = ctrl2,
                 preProc = c("center","scale"))

result5 <- summary(lmTune1) ##RESULT 5 - visual form and surface the most predictive feature against Brysbaert's AoA (aoab)

###########################################
###Analyses C: Network distance measures###
###########################################

# Prepare dataframes for jaccard index 

dnet <- matrix(0, nrow=nrow(wordconc), ncol=nrow(wordconc))
rownames(dnet)<-rownames(wordconc)
colnames(dnet)<-rownames(wordconc)
thresh<-function(v){as.numeric(v>0)}
jaccard_thresh<-function(v,x){sum(thresh(v)*thresh(x))/sum(thresh(thresh(v)+thresh(x)))}
jaccard<-function(v,x){sum(v*x)/sum(v+x)} ## incorporates reporting %

# Jaccard Index
for(i in 1:nrow(dnet)){
  for(j in 1:nrow(dnet)){
    dnet[i,j] <- 1 - (jaccard_thresh(wordconc[i,],wordconc[j,]))
  }
}

dnett<-dnet[rownames(dnet)%in%dt$word,rownames(dnet)%in%dt$word]

result6 <- cor.test(rowSums(dnett),dty, method="spearman") #RESULT 6 - more similar words learned later (Jaccard index on distinctiveness against Brysbaert's AoA)
result6n <- spearman.ci(rowSums(dnett),dty, nrep = 1000)

## Network based on a count of non-shared features
dnet2 <- matrix(0, nrow=nrow(wordconc), ncol=nrow(wordconc))
rownames(dnet2)<-rownames(wordconc)
colnames(dnet2)<-rownames(wordconc)
thresh<-function(v){as.numeric(v>0)}
non_shared_count<-function(v,x){sum(x+v==1)}

for(i in 1:nrow(dnet2)){
  for(j in 1:nrow(dnet2)){
    dnet2[i,j] <- non_shared_count(thresh(wordconc[i,]),thresh(wordconc[j,]))
  }
}
dnett2<-dnet2[rownames(dnet2)%in%dt$word,rownames(dnet2)%in%dt$word]

result7 <- cor.test(rowSums(dnett2),dty , method = "spearman") #RESULT 7 - more distinct words (count of non-shared features) learned earlier (against Brysbaert's AoA)
result7n <- spearman.ci(rowSums(dnett2),dty, nrep = 1000)

## Manhattan network (weighted)
dnet3 <- matrix(0, nrow=nrow(wordconc), ncol=nrow(wordconc))
rownames(dnet3)<-rownames(wordconc)
colnames(dnet3)<-rownames(wordconc)

for(i in 1:nrow(dnet3)){
  for(j in 1:nrow(dnet3)){
    dnet3[i,j] <- sum(abs(wordconc[i,]-wordconc[j,]))
  }
}
dnett3 <- as.data.frame(cbind(rownames(dnet3),rowSums(dnet3)))
names(dnett3) <- c("word","Manhattan")
dnett3 <- inner_join(as.data.frame(dnett3),as.data.frame(aoab))
cor.test(dnett3$Rating.Mean,as.numeric(dnett3$Manhattan))

result8 <- cor.test(dnett3$Rating.Mean,as.numeric(dnett3$Manhattan) , method = "spearman") # RESULT 8 - more distinct words (Manhattan distance between two feature vectors, where distinctiveness is a count of features weighted against the proportion of people reporting said feature) also learned earlier (against Brysbaert's AoA)
result8n <- spearman.ci(dnett3$Rating.Mean,as.numeric(dnett3$Manhattan), nrep = 1000)

#Regression for different distance measures
regdata <- as.data.frame(cbind(dt2$word,rowSums(dnett), rowSums(dnett2), as.numeric(as.character(dnett3$Manhattan))))
names(regdata) <- c("word","jaccard","raw","manhattan")
regdata <- inner_join(regdata, dnett3)

regdata$manhattan <- as.numeric(as.character(regdata$manhattan))
regdata$jaccard <- as.numeric(as.character(regdata$jaccard))
regdata$raw <- as.numeric(as.character(regdata$raw))

regdata <- regdata[,-7]
regdata <- regdata[,-6]
regdata <- regdata[,-5]
regdata.z <- Make.Z(regdata[,-1])
regdata.z <- as.data.frame(regdata.z)
result_regdata <- lm(formula = regdata.z$Rating.Mean ~ regdata.z$manhattan + regdata.z$raw + regdata.z$jaccard)

lm.beta(result_regdata)
confint(result_regdata)

summary(with(regdata, lm(Rating.Mean~manhattan+jaccard+raw)))

## Calculate Manhattan distance within each feature type and correlate with Brysbaert's AoA.

output <- matrix(0, ncol=10)

for (n in 1: length(unique(norms$BR_Label)) )
{
  norms_s <- subset(norms, norms$BR_Label == unique(norms$BR_Label)[n])  
  
  uniq_feats_s <- unique(norms_s$Feature)
  uniq_words_s <- unique(norms_s$Concept)
  
  wordconc_s <- matrix(0, nrow=length(uniq_words_s), ncol=length(uniq_feats_s))
  
  rownames(wordconc_s) <- uniq_words_s
  colnames(wordconc_s) <- uniq_feats_s
  
  for(i in 1:nrow(norms_s)){
    wordconc_s[toString(norms_s$Concept[i]), toString(norms_s$Feature[i])] <- norms_s$Prod_Freq[i]/30
  }
  
  dnets <- matrix(0, nrow=nrow(wordconc_s), ncol=nrow(wordconc_s))
  rownames(dnets)<-rownames(wordconc_s)
  colnames(dnets)<-rownames(wordconc_s)
  
  for(i in 1:nrow(dnets)){
    for(j in 1:nrow(dnets)){
      dnets[i,j] <- sum(abs(wordconc_s[i,]-wordconc_s[j,]))
    }
  }
  
  dnetts <- as.data.frame(cbind(rownames(dnets),rowSums(dnets)))
  names(dnetts) <- c("word","Manhattan")
  dnetts <- inner_join(as.data.frame(dnetts),as.data.frame(aoab))
  bufff <- spearman.ci(dnetts$Rating.Mean,as.numeric(dnetts$Manhattan), nrep = 1000)
  
  outputs <- as.data.frame (cbind( toString(unique(norms$BR_Label)[n]) , nrow(norms_s) , length(uniq_feats_s) , round(length(uniq_feats_s) / nrow(norms_s), digits = 2) , round(as.numeric(spearman.ci(dnetts$Rating.Mean,as.numeric(dnetts$Manhattan), nrep = 1000)$estimate), digits = 2) , round(as.numeric(cor.test(dnetts$Rating.Mean,as.numeric(dnetts$Manhattan), method = "spearman")$p.value) , digits = 3) , p.adjust(as.numeric(cor.test(dnetts$Rating.Mean,as.numeric(dnetts$Manhattan), method = "spearman")$p.value), method = "bonferroni", n = length(unique(norms$BR_Label))) ,paste("[", formatC(round(bufff$conf.int[1],2),2,format="f"), ", ", formatC(round(bufff$conf.int[2],2),2,format="f"), "]", sep="")  , as.numeric(cor.test(dnetts$Rating.Mean,as.numeric(dnetts$Manhattan))[2]) , nrow(dnetts)))
  output <- rbind (output, outputs)
}

output <- output[-1,]
names(output) <- c("feature_type", "n_of_features" , "n_of_unique features"  , "proportion_of_unique_features" , "cor_with_AoA" , "p" , "p_bonf" , "CI1" , "df" , "n")
rownames(output) <- output[,1]
output <- output[,-1]
result9 <- output

## Calculate Manhattan distance within each for each sub-type of visual-form_and_surface and correlate with Brysbaert's AoA.

output2 <- matrix(0, ncol=10)

norms_s1 <- subset(norms, norms$BR_Label == "visual-form_and_surface")

for (n in 1: length(unique(norms_s1$WB_Label)) )
{
  norms_s2 <- subset(norms_s1, norms_s1$WB_Label == unique(norms_s1$WB_Label)[n])
  if (nrow(norms_s2) > 9)
  {
    uniq_feats_s2 <- unique(norms_s2$Feature)
    uniq_words_s2 <- unique(norms_s2$Concept)
    
    wordconc_s2 <- matrix(0, nrow=length(uniq_words_s2), ncol=length(uniq_feats_s2))
    
    rownames(wordconc_s2) <- uniq_words_s2
    colnames(wordconc_s2) <- uniq_feats_s2
    
    for(i in 1:nrow(norms_s2)){
      wordconc_s2[toString(norms_s2$Concept[i]), toString(norms_s2$Feature[i])] <- norms_s2$Prod_Freq[i]/30
    }
    
    dnets2 <- matrix(0, nrow=nrow(wordconc_s2), ncol=nrow(wordconc_s2))
    rownames(dnets2)<-rownames(wordconc_s2)
    colnames(dnets2)<-rownames(wordconc_s2)
    
    for(i in 1:nrow(dnets2)){
      for(j in 1:nrow(dnets2)){
        dnets2[i,j] <- sum(abs(wordconc_s2[i,]-wordconc_s2[j,]))
      }
    }
    
    dnetts2 <- as.data.frame(cbind(rownames(dnets2),rowSums(dnets2)))
    names(dnetts2) <- c("word","Manhattan")
    dnetts2 <- inner_join(as.data.frame(dnetts2),as.data.frame(aoab))
    bufff2 <- spearman.ci(dnetts2$Rating.Mean,as.numeric(dnetts2$Manhattan), nrep = 1000)
    
    outputs2 <- as.data.frame (cbind( toString(unique(norms_s1$WB_Label)[n]) , nrow(norms_s2) , length(uniq_feats_s2) , round(length(uniq_feats_s2) / nrow(norms_s2), digits = 2) , round(as.numeric(spearman.ci(dnetts2$Rating.Mean,as.numeric(dnetts2$Manhattan), nrep = 1000)$estimate), digits = 2) , round(as.numeric(cor.test(dnetts2$Rating.Mean,as.numeric(dnetts2$Manhattan), method = "spearman")$p.value) , digits = 3) , p.adjust(as.numeric(cor.test(dnetts2$Rating.Mean,as.numeric(dnetts2$Manhattan), method = "spearman")$p.value), method = "bonferroni", n = 5) ,paste("[", formatC(round(bufff2$conf.int[1],2),2,format="f"), ", ", formatC(round(bufff2$conf.int[2],2),2,format="f"), "]", sep="")  , as.numeric(cor.test(dnetts2$Rating.Mean,as.numeric(dnetts2$Manhattan))[2]) , nrow(dnetts2)))
    output2 <- rbind (output2, outputs2)
  }
  else
  {
    outputs2 <- as.data.frame (cbind( toString(unique(norms_s1$WB_Label)[n]) , nrow(norms_s2) , "NA" , "NA" , "NA" , "NA" , "NA"))
    output2 <- rbind (output2, outputs2)
  }
}

output2 <- output2[-1,]
names(output2) <- c("feature_type", "n_of_features" , "n_of_unique features"  , "proportion_of_unique_features" , "cor_with_AoA" , "p" , "p_bonf" , "CI" , "df" , "n" )
rownames(output2) <- output2[,1]
output2 <- output2[,-1]
result10 <- output2

############################
###Analyses D: Clustering###
############################

#Establish an empty output dataframe
clustoutput <- matrix(0, ncol=8)

#Set the min and max number of clusters you want to run.
#Currecntly set to 5 - 50, as reported in the paper (might take a while to run).

for(x in 5:50)
{
  
  print(paste("Running clustering analysis for " , x , " clusters."))
  
  #Establish output variables
  maxcor <- 0
  mincor <- 1
  avgcor <- 0
  psig <- 0
  maxpvalue <- 0
  minpvalue <- 0
  avgpvalue <- 0
  clustern <- x
  km <- skmeans(as.matrix(wordconc), clustern, method="genetic")
  kmt <- as.data.frame(table(names(km$cluster),km$cluster)) # makes an output table
  
  #for each cluster, do
  for(y in 1:x)
  { 
    print(paste("Cluster [" , y , "/" , x , "]"))
    subsetwordconc <- subset(wordconc, rownames(wordconc) %in% names(subset(km$cluster, km$cluster == y)))
    dnetclust <- matrix(0, nrow=nrow(subsetwordconc), ncol=nrow(subsetwordconc))
    rownames(dnetclust)<-rownames(subsetwordconc)
    colnames(dnetclust)<-rownames(subsetwordconc)
    
    for(i in 1:nrow(dnetclust)){
      for(j in 1:nrow(dnetclust)){
        dnetclust[i,j] <- sum(abs(subsetwordconc[i,]-subsetwordconc[j,]))
      }
    }
    dnettclust <- as.data.frame(cbind(rownames(dnetclust),rowSums(dnetclust)))
    names(dnettclust) <- c("word","Manhattan")
    dnettclust <- inner_join(as.data.frame(dnettclust),as.data.frame(aoab))
    
    if (nrow(dnettclust) > 5)
    {
      
      coroutcome <- as.numeric(cor.test(dnettclust$Rating.Mean,as.numeric(dnettclust$Manhattan))[4])
      pval <- as.numeric(cor.test(dnettclust$Rating.Mean,as.numeric(dnettclust$Manhattan))[3])
      
      if (pval < 0.05)
      {
        psig <- psig + 1
      }
      
      if(abs(coroutcome) > abs(maxcor))
      {
        maxcor <- coroutcome
        maxpvalue <- pval
      }
      
      if(abs(coroutcome) < abs(mincor))
      {
        mincor <- coroutcome
        minpvalue <- pval
      }
      
      avgcor <- avgcor + coroutcome
      avgpvalue <- avgpvalue + pval
      
    }
    
  }
  
  psigprop <- psig / x
  avgcor <- avgcor / x
  avgpvalue <- avgpvalue / x
  
  clustoutputs <- as.data.frame (cbind( as.numeric(x) , as.numeric(mincor) , as.numeric(maxcor) , as.numeric(avgcor), as.numeric(minpvalue) , as.numeric(maxpvalue) , as.numeric(avgpvalue) , as.numeric(psigprop) ))
  clustoutput <- rbind (clustoutput, clustoutputs)
  
}

clustoutput <- clustoutput[-1,]
names(clustoutput) <- c("Cluster n", "weakest cor" , "strongest cor" , "average cor" , "p for weakes cor" , "p for strongest cor" , "average p" , "% of p < .05" )

result11 <- clustoutput