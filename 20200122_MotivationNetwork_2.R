install.packages("qgraph") # version 1.6.4
install.packages("biclust")
install.packages("IsingFit")
install.packages("plyr")
install.packages("Hmisc")
install.packages("ggm")
install.packages("ggplot2")
# for vif: 
#install.packages("tidyverse")
#install.packages("caret")
#install.packages("car")
install.packages("ltm") # für Mair et al. (2015) Code
install.packages("homals") # für Mair et al. (2015) Code
install.packages("lattice")
install.packages("simex")
install.packages("memisc")
install.packages("corpcor")

library(corpcor)
library(ppcor)
library(memisc)
library(MASS)
library(memisc)
library(effects)
library(simex)
library(car)
library(tidyverse)
library(caret)
library(Hmisc)
library(ggm)
library(ggplot2)
library(polycor)
#library(plyr)
library(IsingFit)
library(biclust)
library(qgraph)
library(ltm)
library(homals)
# library(mice)
# library(VIM)
library(lattice)
library(igraph)

###
###   Daten einlesen
###

datcompleteraw <- read.csv("Mair.pnas.1506047112.sd01.csv", sep = ";", header = TRUE)
View(datcompleteraw)
dim(datcompleteraw)  ## 1448 Personen insgesamt

## ------ Missing Values 
nmiss <- apply(datcompleteraw, 1, function(x) sum(is.na(x)))        ## check missings
table(nmiss)  
## - 310 persons quit the questionnaire immediately
## - 51 scrolled through without answering

## Eliminieren von Personen mit all NA
d <- datcompleteraw[which(nmiss < 176),]
dim(d) 
# noch 1087 Personen

###
###   Daten zur Motivation und Partizipation vorbereiten und Kontrolle fehlender Werte
###

# motpartdatcompleteraw1 <- datcompleteraw[1:4]
motdatcompleteraw2 <- d[72:107] #Motivationitems
partdatcompleteraw3 <- d["v_192"] #Participationitem1
partdatcompleteraw4 <- d[, c("v_210", "v_211")] #Participationsitems2&3 (v_210, v_211)
motpartdatwithoutfullmissings1 <- cbind(motdatcompleteraw2, partdatcompleteraw3, partdatcompleteraw4)
dim(motpartdatwithoutfullmissings1)
motpartdatwithoutfullmissings2 <- cbind(motdatcompleteraw2, partdatcompleteraw3)
dim(motpartdatwithoutfullmissings2)
motpartdatwithoutfullmissings3 <- cbind(motdatcompleteraw2, partdatcompleteraw4)

motpartdatwithoutfullmissings <- motpartdatwithoutfullmissings1[which(colSums(apply(motpartdatwithoutfullmissings1, 1, is.na)) != 39),]
dim(motpartdatwithoutfullmissings) # jetzt bleiben noch 906 Personen, die nicht Nichts ausgefüllt haben in Mot und/oder Partizipation

omittedall <- na.omit(motpartdatwithoutfullmissings)
dim(omittedall) #720 Personen

###
###   Checken der Zahl (720 Personen), weil Mair et al. (2015) angegeben, dass 767 Personen alles eingetragen haben.   
###

omittedmot <- na.omit(motdatcompleteraw2)
dim(omittedmot) #767 Personen - Motivationitems

omittedpart1 <- na.omit(partdatcompleteraw3)
dim(omittedpart1) #895 Personen - Number of Packages

omittedpart2 <- na.omit(partdatcompleteraw4)
dim(omittedpart2) #802 Personen - Mailing lists 

omittedmotv_192 <- na.omit(motpartdatwithoutfullmissings2)
dim(omittedmotv_192) #763 Personen - Motivation & Number of Packages

omittedmotv_210 <- na.omit(cbind(motdatcompleteraw2, d[141]))
dim(omittedmotv_210) #747 Personen - Motivation & Mailing lists

omittedmotv_211 <- na.omit(cbind(motdatcompleteraw2, d[142]))
dim(omittedmotv_211) #727 Personen - Motivation & R conferences

omittedmotv_210_v_211 <- na.omit(motpartdatwithoutfullmissings3)
dim(omittedmotv_210_v_211) #720 Personen - Motivation & Mailing Lists & R Conferences

########################### Regressionsmodell ##################################################################

###
###  Code aus Mair et al. (2015) mit einigen Notizen von mir, 
###  um herauszufinden, woher die Differenz zu den 720 Individuen kommt.
###  (Siehe Fußnote in Subsection 4.1.2)
###  Kann übersprungen werden bis Zeile 437 (--> Network Analysis)
###

###
###    Dependent Variables
###

RMotivation <- d[, c("v_210", "v_211", "v_192")] # RMotivation ist erstmal nur die drei dependent Variables
# RMotivation 
names(RMotivation) <- c("lists", "meet", "npkgs") # Umbenennen der Spalten
# RMotivation
## ?factor # vector wird zu factor 
RMotivation$lists <- factor(RMotivation$lists, levels = 1:2, labels = c("no", "yes"))
# RMotivation$lists # aus 1 und 2 wird "no" und "yes" 
RMotivation$meet <- factor(RMotivation$meet, levels = 1:2, labels = c("no", "yes"))
# RMotivation
Dep_RMotivation <- cbind(RMotivation, d[72:107])
omittedDep_RMotivation <- na.omit(Dep_RMotivation)
dim(omittedDep_RMotivation) # Noch immmer 720 Individuals


## ---------------------------------------- IRT Analysis ---------------------------------------
## Unidimensionality is checked with a categorical PCA (homals package)
## Adjustment for multiple itemfit correction since all of the hypothesis to be tested 
## are connected/related within the psychological measures to be used. 
## We use alpha = 0.05 for each psychological measure:
## alpha_adj = alpha / total # of items of the corresponding subscale.
## For each itemfit analysis we use the Monte Carlo simulation approach to approximate the p-values. 

?colSums
## ------ Reinholt Motivation Scale ------
## remove full NAs in Motivation scale 
reinholt <- d[,72:107] - 1
reinholt <- reinholt[which(colSums(apply(reinholt, 1, is.na)) != 36),]
dim(reinholt)                                     ## 852 persons, 36 items

## --- extreme extrinsic motivation
## - external regulation: 01-08
## - introjection: 09-12
mextrinsic <- reinholt[, c(
  paste0("eim_0", 1:8),
  c("eim_09", paste0("eim_1", 0:2))
)]
dim(mextrinsic)
## 12 items

plot(homals(na.omit(mextrinsic), level = "ordinal"))             ## unidimensionality: eim_05 and eim_06 extremely suspicious
mextrinsic <- mextrinsic[, -grep("eim_05", names(mextrinsic))]   ## eliminate eim_05 and eim_06 right away
mextrinsic <- mextrinsic[, -grep("eim_06", names(mextrinsic))] 
fit.mextrinsic <- ltm(mextrinsic ~ z1)                 ## 2-PL fit
set.seed(123)
# ifit.mextrinsic <- item.fit(fit.mextrinsic, simulate.p.value = TRUE)  ## compute Q1 statistic
alpha_reinholt <- 0.05/dim(mextrinsic)[2]              ## alpha correction
# which(ifit.mextrinsic$p.value < alpha_reinholt)        ## itemfit OK

## --- well-internalized extrinsic motivation/moderated intrinsic motivation
## - identification: 13-17
## - integration: 18-22
## - obligation based motivation: 23-27
## - self-reinforcement: 28-31
mhybrid <- reinholt[, c(
  paste0("eim_1", 3:7),
  c(paste0("eim_1", 8:9), paste0("eim_2", 0:2)),
  paste0("eim_2", 3:7),
  c(paste0("eim_2", 8:9), paste0("eim_3", 0:1))
)]
dim(mhybrid)                                       ## 19 items

plot(homals(na.omit(mhybrid), level = "ordinal"))  ## unidimensionality check
fit.mhybrid <- ltm(mhybrid ~ z1)                   ## 2-PL fit
set.seed(123)
# ifit.mhybrid <- item.fit(fit.mhybrid, simulate.p.value = TRUE)  ## compute Q1 statistic
alpha_reinholt <- 0.05/dim(mhybrid)[2]             ## alpha correction
# which(ifit.mhybrid$p.value < alpha_reinholt)       ## itemfit OK

## --- extreme intrinsic motivation
## - enjoyment based intrinsic motivation: 32-36
mintrinsic <- reinholt[, paste0("eim_3", 2:6)]
dim(mintrinsic)                              ## 5 items

plot(homals(na.omit(mintrinsic), level = "ordinal"))         ## unidimensionality check
fit.mintrinsic <- ltm(mintrinsic ~ z1)                       ## 2-PL fit
set.seed(222)
# ifit.mintrinsic <- item.fit(fit.mintrinsic, G = 3, simulate.p.value = TRUE) ## Q1 statistic; G lowered due to convergence problems
alpha_reinholt <- 0.05/dim(mintrinsic)[2]                    ## alpha correction
# which(ifit.mintrinsic$p.value < alpha_reinholt)              ## itemfit OK

?apply
?paste
## ------------------------------ Compute IRT Person Parameters -------------------------- 
## some convenience functions for computing the person parameters
# für unsere Daten M (Motivation) werden die Reihen immer durch ein " " getrennt. 
tagify <- function(M) apply(M, 1, paste, collapse = " ") 
##?factor.scores # --> Was ist der einzelnen Person Wert der latenten Variable? 
## fs: Gibt mir eine Liste mit den "score.dat" --> ein data.frame mit bspw. factor scores und Standardabweichungen
## 
## method = "EB" - Außerdem verwenden sie für die Factor Scores Empirical Bayes
#?ltm
## ltm(data ~ z1) - Eine latente Variable beschreibt mir die factor scores in den Daten
#?match
fscore <- function(data) {
  fs <- factor.scores(ltm(data ~ z1), method = "EB")
  ix <- match(tagify(data), tagify(fs$score.dat[, 1:nrow(fs$coef)]))
  list(
    scores = fs$score.dat$z1[ix],
    serrs  = fs$score.dat$se.z1[ix],
    names  = rownames(data)
  )
}
## Hier werden jetzt die fscores ausgerechnet, aber es sind noch NA Values mit drin. 
## Wieso ist das möglich? 
## Weil die ltm function unter MAR durch na.action = NULL die observed data miteinbezieht. 
## compute person parameters and associated standard errors
pp <- list(
  mextrinsic = fscore(mextrinsic),
  mhybrid = fscore(mhybrid),
  mintrinsic = fscore(mintrinsic)
)
View(pp)
summary(pp)

## ---------------------------------- Final Preparation Step for Data used in GLMs -------------------------
#View(pp)
## person parameters for scales and associated standard errors
for(i in names(pp)) {
  RMotivation[pp[[i]]$names, i] <- pp[[i]]$scores
  RMotivation[pp[[i]]$names, paste0(i, ".se")] <- pp[[i]]$serrs
}
# vorherige Zeilen nehmen wir für jede latent Variable im pp Element (früher e.g.: wtask; jetzt: mextrinsic, mhybrid und mintrinsic) 
# die scores Spalte aus dem pp (person parameter) Element die laufende Nummer dieser Variable 
# und die Spalte (also das "i" --> # die latente Variable) 
# und fügt dort die scores und serrs (standarderrors) ein 
# RMotivation
## socio-demographic variables
RMotivation$gender    <- factor(d$v_216, levels = 1:2, labels = c("male", "female"))   ## gender
RMotivation$phd       <- factor(d$v_220, levels = 0:1, labels = c("no", "yes"))        ## PhD 
RMotivation$statseduc <- factor(d$v_222, levels = 0:1, labels = c("no", "yes"))        ## statistical education
RMotivation$fulltime  <- factor(d$v_231, levels = 0:1, labels = c("no", "yes"))        ## working full-time
RMotivation$academia  <- factor(d$v_237, levels = 0:1, labels = c("no", "yes"))        ## working in academia
RMotivation$statswork <- factor(d$v_245, levels = 0:1, labels = c("no", "yes"))        ## working in the statistical area
#RMotivation # jetzt haben wir das halt auch noch alles in meiner Übersicht. 

## save(RMotivation, file = "RMotivation.rda")

## ------------------------------ Part II: This code reproduces all the analyses in the article -------------------------------------------
## load("RMotivation.rda")  

RMotivation <- na.omit(RMotivation)
dim(RMotivation) # Hier andere Anzahl an Reihen (Individuen) aber auch an Spalten (6 Kontrollvariablen, 
# 3 dependent, 3 motivation + 3 motivation.se), weil bereits factor scores und nicht mehr die einzelnen 
# Indikatoren in der list sind. 
psychometrics <- c("mextrinsic", "mhybrid", "mintrinsic")  
psychometricsSE <- paste0(psychometrics, ".se")
demographics <- c("phd", "statseduc", "fulltime", "academia", "statswork")

## --------------------------- descriptive analysis ----------------------------
## histogram: number of packages (Fig. S2)
hist(RMotivation$npkgs, main = "Histogram R Packages", xlab = "Number of R Packages", breaks = -1:max(RMotivation$npkgs, na.rm = TRUE) + 0.5)
plot(hist(RMotivation$npkgs, main = "Histogram R Packages", xlab = "Number of R Packages", breaks = -1:max(RMotivation$npkgs, na.rm = TRUE) + 0.5))

## --------------------------- regression models -------------------------------

## --- negative-binomial regression: number of packages
if(file.exists("models-npkgs.rda")) {
  load("models-npkgs.rda")
} else {
  formulaGLMnpkgs <- as.formula(paste("npkgs ~", paste(c(psychometrics, demographics), collapse = " +")))
  fitnpkgsNB <- glm.nb(formulaGLMnpkgs, RMotivation, x = TRUE, y = TRUE)
  fitnpkgs <- glm(formulaGLMnpkgs, RMotivation, family = negative.binomial(fitnpkgsNB$theta), x = TRUE, y = TRUE)
  
  ## stepwise selection
  fitnpkgsStep <- step(fitnpkgs, formulaGLMnpkgs, trace = 0)
  fitnpkgsStepplot <- fitnpkgsStep
  
  ## SIMEX versions
  ME <- RMotivation[, psychometricsSE]
  fitnpkgsSimex <- simex(fitnpkgs, SIMEXvariable = psychometrics, measurement.error = ME, asymptotic = FALSE) 
  psychoind <- psychometrics %in% names(coef(fitnpkgsStep))
  ME1 <- RMotivation[, psychometricsSE[psychoind]]
  fitnpkgsStepSimex <- simex(fitnpkgsStep, SIMEXvariable = psychometrics[psychoind], measurement.error = ME1, asymptotic = FALSE) 
  fitnpkgs$xlevels <- NULL
  fitnpkgsStep$xlevels <- NULL
  
  ## save results
  save(fitnpkgs, fitnpkgsSimex, fitnpkgsStep, fitnpkgsStepSimex, fitnpkgsStepplot, file = "models-npkgs.rda")
}
View(fitnpkgs) # entspricht meine Full-ML Werten aus Mair et al. (2015)
View(fitnpkgsStepplot) # entspricht meinen Step (ML) Werten 
glm(formulaGLMnpkgs, RMotivation, family = negative.binomial(fitnpkgsNB$theta), x = TRUE, y = TRUE)

## regression table (Table S1) --> mtable funktioniert leider nicht mit meiner R-Version. 
tabnpkgs <- mtable(
  "Full (ML)" = fitnpkgs,
  "Full (SIMEX)" = fitnpkgsSimex,
  "Step (ML)" = fitnpkgsStep,
  "Step (SIMEX)" = fitnpkgsStepSimex,
  summary.stats = FALSE) 
toLatex(tabnpkgs)

## effect plots 
plot(allEffects(fitnpkgsStepplot), ylab = "Number of packages", type = "response", ylim = c(1.6, 4))

## --- logistic regression (participation in lists)
if(file.exists("models-lists.rda")) {
  load("models-lists.rda")
} else {
  ## formula
  formulaGLMlists <- as.formula(paste("lists ~", paste(c(psychometrics, demographics), collapse = " +")))
  
  ## full logistic regression
  fitlists <- glm(formulaGLMlists, RMotivation, family = binomial(), x = TRUE, y = TRUE)
  
  ## stepwise selection
  fitlistsStep <- step(fitlists, formulaGLMlists, trace = 0)
  fitlistsStepplot <- fitlistsStep
  
  ## SIMEX versions
  ME <- RMotivation[, psychometricsSE]
  fitlistsSimex <- simex(fitlists, SIMEXvariable = psychometrics, measurement.error = ME, asymptotic = FALSE) 
  psychoind <- psychometrics %in% names(coef(fitlistsStep))
  ME1 <- RMotivation[, psychometricsSE[psychoind]]
  fitlistsStepSimex <- simex(fitlistsStep, SIMEXvariable = psychometrics[psychoind], measurement.error = ME1, asymptotic = FALSE) 
  fitlists$xlevels <- NULL
  fitlistsStep$xlevels <- NULL
  
  ## save results
  save(fitlists, fitlistsSimex, fitlistsStep, fitlistsStepSimex, fitlistsStepplot, file = "models-lists.rda")
}

## regression table (Table S2)
tablists <- mtable(
  "Full (ML)" = fitlists,
  "Full (SIMEX)" = fitlistsSimex,
  "Step (ML)" = fitlistsStep,
  "Step (SIMEX)" = fitlistsStepSimex,
  summary.stats = FALSE)  
toLatex(tablists)
View(tablists)

## effect plots 
plot(allEffects(fitlistsStepplot), ylab = "Probability of mailing list participation", type = "response", ylim = c(0.35, 0.75))


## --- logistic regression (conference participation)
if(file.exists("models-meet.rda")) {
  load("models-meet.rda")
} else {
  ## formula
  formulaGLMmeet <- as.formula(paste("meet ~", paste(c(psychometrics, demographics), collapse = " +")))
  
  ## full logistic regression
  fitmeet <- glm(formulaGLMmeet, RMotivation, family = binomial(), x = TRUE, y = TRUE)
  
  ## stepwise selection
  fitmeetStep <- step(fitmeet, formulaGLMmeet, trace = 0)
  fitmeetStepplot <- fitmeetStep
  
  ## SIMEX versions
  ME <- RMotivation[, psychometricsSE]
  fitmeetSimex <- simex(fitmeet, SIMEXvariable = psychometrics, measurement.error = ME, asymptotic = FALSE) 
  psychoind <- psychometrics %in% names(coef(fitmeetStep))
  ME1 <- RMotivation[, psychometricsSE[psychoind]]
  fitmeetStepSimex <- simex(fitmeetStep, SIMEXvariable = psychometrics[psychoind], measurement.error = ME1, asymptotic = FALSE) 
  fitmeet$xlevels <- NULL
  fitmeetStep$xlevels <- NULL
  
  ## save results
  save(fitmeet, fitmeetSimex, fitmeetStep, fitmeetStepSimex, fitmeetStepplot, file = "models-meet.rda")
}


## regression table (Table S3)
tabmeet <- mtable(
  "Full (ML)" = fitmeet,
  "Full (SIMEX)" = fitmeetSimex,
  "Step (ML)" = fitmeetStep,
  "Step (SIMEX)" = fitmeetStepSimex,
  summary.stats = FALSE) 
toLatex(tabmeet)


## effect plots 
plot(allEffects(fitmeetStepplot), ylab = "Probability of conference participation", type = "response", ylim = c(0.13, 0.45))


###
###   Negative Binomial Regression aus Mair Code --> Motivation auf Participation. 
###

plot(fitnpkgsNB)
vif(fitnpkgsNB)

###
###     Extrahieren der einzelnen ltm-models für die latenten Variablen mextrinsic, mhybrid, mintrinsic 
###

tagify <- function(M) apply(M, 1, paste, collapse = " ") 
count <- 0
fscore_MZ <- function(data) {
  #print(count)
  fs <- factor.scores(ltm(data ~ z1), method = "EB")
  print(paste("ltm", count, sep = ""))
  assign(paste("ltm", count, sep = ""), ltm(data ~ z1), envir = .GlobalEnv)
  ix <- match(tagify(data), tagify(fs$score.dat[, 1:nrow(fs$coef)]))
  list(
    scores = fs$score.dat$z1[ix],
    serrs  = fs$score.dat$se.z1[ix],
    names  = rownames(data)
  )
  count <<- count + 1
  #print(count)
}

pp_MZ <- list(
  mextrinsic = fscore_MZ(mextrinsic),
  mhybrid = fscore_MZ(mhybrid),
  mintrinsic = fscore_MZ(mintrinsic)
)
View(pp_MZ)

##### END ######### Mair et al. (2015) Code ################


######################## Network analysis #################################################################

###
### Erstellen des Netzwerks
###

omittedall <- na.omit(motpartdatwithoutfullmissings) # Datenbasis ohne Zeilen mit NA values
dim(omittedall) # 720 Personen - All Motivation items, and all three participation items
omittedall

coromittedall <-cor(omittedall) # Pearson Correlations meiner Daten

coromittedalltable <- as.data.frame(as.table(coromittedall)) # 
View(coromittedalltable)

# Korrelationstabelle komplett: 
coromittedall_clean <- coromittedalltable[which(coromittedalltable[,3] < 1, TRUE), ] 
View(coromittedall_clean)
arrangedomittedall_clean = arrange(coromittedall_clean, -Freq)

# Nur noch jede zweite Zeile --> Weil Korrelationen immer doppelt: z.B.: eim_30 mit eim_33 und eim_33 mit eim_30
toDelete <- seq(2, nrow(arrangedomittedall_clean), 2)

arrangedomittedall_clean_2 <- arrangedomittedall_clean[-toDelete,]
View(arrangedomittedall_clean_2)
export <- arrangedomittedall_clean_2 # Exportieren der Datei für den Appendix
write.table(export, "20200124_PearsonCorrelations.txt", row.names=FALSE)
write.csv2(export, "20200124_PearsonCorrelations.csv")


##
## Mergen der nodes mit Korrelationen über .5
##

omittedall <- na.omit(motpartdatwithoutfullmissings)
dim(omittedall) # Sicherheitscheck
eim21_18 <- cbind(omittedall$eim_21, omittedall$eim_18)
eim15_16 <- cbind(omittedall$eim_15, omittedall$eim_16)
eim36_34 <- cbind(omittedall$eim_36, omittedall$eim_34)
eim30_33 <- cbind(omittedall$eim_30, omittedall$eim_33)
m21_18 <- rowMeans(eim21_18)
m15_16 <- rowMeans(eim15_16)
m36_34 <- rowMeans(eim36_34)
m30_33 <- rowMeans(eim30_33)

# Löschen der Knotenpunkte, die durch einen merged Node im Netzwerk repräsentiert werden. 
omittedall$eim_21 <- NULL #1
omittedall$eim_18 <- NULL #2
omittedall$eim_15 <- NULL #3
omittedall$eim_16 <- NULL #4
omittedall$eim_36 <- NULL #5
omittedall$eim_34 <- NULL #6
omittedall$eim_30 <- NULL #7
omittedall$eim_33 <- NULL #8

View(omittedall)
omittedmerged <- cbind(omittedall, m21_18, m15_16, m36_34, m30_33) # Merged Nodes wieder in Daten einbinden
dim(omittedmerged)

###
### Nochmal Korrelationen checken
###

coromittedmerged <- cor(omittedmerged)
coromittedmerged_table <- as.data.frame(as.table(coromittedmerged))
View(coromittedmerged_table) # Passt. Keine Korrelationen über .5 
sd(coromittedmerged_table$Freq)
summary(coromittedmerged_table$Freq)

###
### --> coromittedmerged als finale Korrelationsmatrix für die Konstruktion meines Netzwerks
###

qgraph(coromittedmerged, 
       layout = "spring", 
       graph = "pcor") # graph = pcor macht es zu einer partial correlation matrix

qgraph(coromittedmerged, 
       layout = "spring") # graph auf Basis von Pearson Korrelations

qgraph(coromittedmerged, 
       layout = "spring", 
       sampleSize = 720,
       graph = "glasso") # graph nutzt Glasso zur Bestimmung der edges

# Als PDF speichern für bessere Auflösung
qgraph(coromittedmerged, 
       layout = "spring", 
       graph = "pcor", 
       filetype = "pdf")

qgraph(coromittedmerged, 
       layout = "spring", 
       filetype = "pdf") 

qgraph(coromittedmerged, 
       layout = "spring", 
       sampleSize = 720,
       graph = "glasso", 
       filetype = "pdf")

PcorGraph <- qgraph(coromittedmerged, 
                    layout = "spring", 
                    graph = "pcor")


GlassoGraph <- qgraph(coromittedmerged, 
                      layout = "spring", 
                      sampleSize = 720,
                      graph = "glasso")

###
### Nun filtere ich aus der Qgraph Datei meine Partial Correlations Matrix heraus
###

Pcormatrix <- as.data.frame(PcorGraph$Edgelist) # Problem: Variablen nicht wie urpsrünglich eim_12 etc., sondern durchlaufend benannt, wobei 1 nicht eim_01 ist. 
PcorVarNames <- PcorGraph$graphAttributes$Nodes$names # Hier finde ich aber meine Variablen Namen zu den jeweiligen Zeilen.
Needed <- cbind(c(1:35), PcorVarNames) # Einfügen einer Spalte, die das Mergen anhand der Spaltennummer ermöglicht. 
merge1 <- merge(Needed, Pcormatrix, by.x = "V1", by.y = "from")
pcor_usable <- merge(Needed, merge1, by.x = "V1", by.y = "to") # Jetzt habe ich die Variablennamen und die Partial Correlations jeweils passend zusammen und nicht mehr nichtssagende durchlaufende Nummern

# Löschen der Spalten mit den durchlaufenden Nummern, weil nun ja die richtige Bezeichnung in meiner pcor_usable Matrix ist
pcor_usable$V1 <- NULL
pcor_usable$V1.y <- NULL
pcor_usable$directed <- NULL
pcor_usable$bidirectional <- NULL
pcor_usable
orderedpcor_usable <- pcor_usable[order(-pcor_usable$weight), ] # Sortieren nach Weight
View(orderedpcor_usable)
summary(orderedpcor_usable)
quantile(orderedpcor_usable$weight, c(.98))
?quantile
write.csv2(orderedpcor_usable, "20200130_PartialCorrelations.csv") # Exportieren meiner Partial Correlations 

Pcor_table <- as.data.frame(as.table(Pcormatrix$weight))
View(Pcor_table)
sd(Pcor_table$Freq) # Standardabweichung
summary(Pcor_table$Freq)

###
### Centrality Maße im Netzwerk
###

dim(omittedmerged)
centralityPlot(GlassoGraph, include = "All")
centralityPlot(GlassoGraph, include = "ExpectedInfluence")
EI_Plot <- centralityPlot(GlassoGraph, include = "ExpectedInfluence")
View(EI_Plot$data) # Expected Influence - Werte
View(centralityTable(GlassoGraph))
View(centralityTable(PcorGraph))

###
### Ermitteln des Durchschnitts der Centrality Werte der einzelnen Formen von Motivation
###

EI_name_value <- data.frame(data=EI_Plot$data$node, EI_Plot$data$value) # Dataframe erstellen mit Knotenpunkten und den jeweiligen EI-Werten dazu
EI_name_value
EM_EI <- EI_name_value[c(1:12), ] # Nur EM Werte
EM_EI
wEM_mIM_EI <- rbind(EI_name_value[c(13:26), ], EI_name_value[c(29:30), ]) # Nur wEM-mIM Werte
wEM_mIM_EI
IM_EI <- rbind(EI_name_value[c(27:28), ], EI_name_value[c(31:32), ]) # Nur IM Werte
IM_EI
# Jeweils den Durschnitt berechnen:
sum(EM_EI$EI_Plot.data.value)
avsumEM <- sum(EM_EI$EI_Plot.data.value)/nrow(EM_EI)
avsumEM
sum(wEM_mIM_EI$EI_Plot.data.value)
avsumwEM_mIM <-  sum(wEM_mIM_EI$EI_Plot.data.value)/nrow(wEM_mIM_EI)
avsumwEM_mIM
sum(IM_EI$EI_Plot.data.value)
avsumIM <- sum(IM_EI$EI_Plot.data.value)/nrow(IM_EI)
avsumIM

avsumEM
avsumwEM_mIM
avsumIM

summary(EI_name_value$EI_Plot.data.value)
sd(EI_name_value$EI_Plot.data.value)
EI_name_value$EI_Plot.data.value

#sd(matrix(1:2, nrow = 5000000, ncol = 1)) # Kontrolle, ob die sd function auch echt funktioniert

###
### Test, ob eine Edge zwischen 24 und 210 im GLASSO Graph erscheint, wenn ich nodes 23 und 27 rausnehmen. 
###

omittedmerged_minus2427 <- omittedmerged

omittedmerged_minus2427$eim_24 <- NULL
omittedmerged_minus2427$eim_27 <- NULL
coromittedmerged_minus2427 <- cor(omittedmerged_minus2427)

qgraph(coromittedmerged_minus2427, 
       layout = "spring", 
       sampleSize = 720,
       graph = "glasso", 
       filetype = "pdf") # Ja, sie erscheint. 

###
### Netzwerkkonstruktion mit nur dem wichtigsten Participation item: Packages
###
               
omittedmerged # Meine Ausgangsbasis
omittedmerged_onlypack <- omittedmerged
omittedmerged_onlypack$v_210 <- NULL
omittedmerged_onlypack$v_211 <- NULL
dim(omittedmerged_onlypack) # Check

coromittedmerged_onlypack <- cor(omittedmerged_onlypack)

qgraph(coromittedmerged_onlypack, 
       layout = "spring", 
       sampleSize = 720,
       graph = "glasso", 
       filetype = "pdf") 

###
### Netzwerkkonstruktion mit Nodes merged auf Basis meiner Results und den theoretischen Überlegungen
###

# Alles analog zu obigem Mergen der Nodes im Netzwerk 
omittedall <- na.omit(motpartdatwithoutfullmissings) # zuvor in omittedall ja bereits Knotenpunkte rausgelöscht und durch gemergte ersetzt - daher wieder die anfängliche Omittedall Matrix wiederhergestellt
dim(omittedall) # Sicherheitscheck
eim_09_11 <- cbind(omittedmerged$eim_09, omittedmerged$eim_11)
eim_08_04_12 <- cbind(omittedmerged$eim_08, omittedmerged$eim_04, omittedmerged$eim_12)
eim_21_20 <- cbind(omittedmerged$m21_18, omittedmerged$m21_18, omittedmerged$eim_20) # 2x merged nodes, weil diese ja bereits aus zwei Elementen bestehen (um den Durschnitt richtig zu berechnen) 
eim_29_31 <- cbind(omittedmerged$eim_31, omittedmerged$eim_29)
eim_32_36 <- cbind(omittedmerged$m36_34, omittedmerged$m36_34, omittedmerged$eim_32) # siehe oben
eim_23_24_27 <- cbind(omittedmerged$eim_23, omittedmerged$eim_24, omittedmerged$eim_27)

m09_11 <- rowMeans(eim_09_11)
m08_04_12 <- rowMeans(eim_08_04_12)
m21_20 <- rowMeans(eim_21_20)
m29_31 <- rowMeans(eim_29_31)
m32_36 <- rowMeans(eim_32_36)
m23_24_27 <- rowMeans(eim_23_24_27)

omittedmerged2 <- omittedmerged

omittedmerged2$eim_09 <- NULL #9
omittedmerged2$eim_11 <- NULL #10
omittedmerged2$eim_08 <- NULL #11
omittedmerged2$eim_04 <- NULL #12
omittedmerged2$eim_12 <- NULL #13
omittedmerged2$m21_18 <- NULL #14
omittedmerged2$eim_20 <- NULL #15
omittedmerged2$eim_31 <- NULL #16
omittedmerged2$eim_29 <- NULL #17
omittedmerged2$m36_34 <- NULL #18
omittedmerged2$eim_32 <- NULL #19
omittedmerged2$eim_23 <- NULL #20
omittedmerged2$eim_24 <- NULL #21
omittedmerged2$eim_27 <- NULL #22

dim(omittedall)
omittedmerged_th <- cbind(omittedmerged2, m09_11, m08_04_12, m21_20, m29_31, m32_36, m23_24_27)
dim(omittedmerged_th)

coromittedmerged_th <- cor(omittedmerged_th)

qgraph(coromittedmerged_th, 
       layout = "spring", 
       sampleSize = 720,
       graph = "glasso")

qgraph(coromittedmerged_th, 
       layout = "spring", 
       sampleSize = 720,
       graph = "glasso", 
       filetype = "pdf")


###
### Trying to test the Connectivity in less and more participating individuals
###

omittedmerged
summary(omittedmerged$v_192)
names(omittedmerged)
omittedmerged_gr2 <- matrix(data = NA, ncol = 35) # Erstellen einer leeren Matrix mit der richtigen Anzahl Columns für den späteren For-Loop
omittedmerged_gr2
omittedmerged_max2 <- matrix(data = NA, ncol = 35) # Erstellen einer leeren Matrix mit der richtigen Anzahl Columns für den späteren For-Loop
omittedmerged_gr2 <- as.data.frame(omittedmerged_gr2)
omittedmerged_max2 <- as.data.frame(omittedmerged_max2)
colnames(omittedmerged_gr2) <- names(omittedmerged) # Richtige Namen hinzufügen
colnames(omittedmerged_max2) <- names(omittedmerged) # Richtige Namen hinzufügen
index <- 0 # Meine Count Variable

# Der folgende For-Loop geht mir durch meine Spalte v_192 durch und fügt die Zeilen, 
# für die v_192 größer als 2 ist / kleiner als 2 ist in die jeweiligen Matrizen ein
for (i in omittedmerged$v_192) {
  if (i > 2) {
    index <- index + 1
    omittedmerged_gr2 <- rbind(omittedmerged_gr2, omittedmerged[index, ])
    print("bigger!")
  }
  else {
    index <- index + 1
    omittedmerged_max2 <- rbind(omittedmerged_max2, omittedmerged[index, ])
    print("nope")
  }
}
dim(omittedmerged)
index # test, ob durch alle Reihen durchgegangen wurde 
dim(omittedmerged_gr2) # Hier jeweils zusätzlich eine leere Zeile mit in meiner leeren Matrix
dim(omittedmerged_max2) # Hier jeweils zusätzlich eine leere Zeile mit in meiner leeren Matrix
omittedmerged_gr2 <- na.omit(omittedmerged_gr2) # Hier lösche ich diese 
omittedmerged_max2 <- na.omit(omittedmerged_max2) # Hier lösche ich diese 
dim(omittedmerged_gr2) # 260 Partizipierende mit mehr als 1 Packages --> Definitiv zu wenig
dim(omittedmerged_max2) # 460 Partizipierende mit 1 (oder weniger Packages) --> Wohl auch zu grenzwertig


# Konstruktion der Netzwerke

coromittedmerged_gr2 <- cor(omittedmerged_gr2)
coromittedmerged_max2 <- cor(omittedmerged_max2)
coromittedmerged_gr2_table <- as.data.frame(as.table(coromittedmerged_gr2))
coromittedmerged_max2_table <- as.data.frame(as.table(coromittedmerged_max2))

GlassoGraph_gr2 <- qgraph(coromittedmerged_gr2, 
                      layout = "spring", 
                      sampleSize = 260,
                      graph = "glasso", 
                      filetype = "pdf")

GlassoGraph_max2 <- qgraph(coromittedmerged_max2, 
                           layout = "spring", 
                           sampleSize = 460,
                           graph = "glasso", 
                           filetype = "pdf")

GlassoGraph_gr2 <- qgraph(coromittedmerged_gr2, 
                          layout = "spring", 
                          sampleSize = 260,
                          graph = "pcor") # Unterschied wird in Pcor Darstellung nicht gut deutlich

GlassoGraph_max2 <- qgraph(coromittedmerged_max2, 
                           layout = "spring", 
                           sampleSize = 460,
                           graph = "pcor") # Unterschied wird in Pcor Darstellung nicht gut deutlich

GlassoGraph_gr2 <- qgraph(coromittedmerged_gr2, 
                          layout = "spring", 
                          sampleSize = 260,
                          graph = "glasso")

GlassoGraph_max2 <- qgraph(coromittedmerged_max2, 
                           layout = "spring", 
                           sampleSize = 460,
                           graph = "glasso")

# Hier sind meine trans_target Werte auf Basis des GLASSO Graphs für die connectivity Analyse
# Da der Graph auf Basis von partial correlations ein saturated Graph ist, ist es nicht sinnvoll, die Werte aus diesem zu nehmen. 
smallworldness(GlassoGraph_gr2)  
smallworldness(GlassoGraph_max2)



