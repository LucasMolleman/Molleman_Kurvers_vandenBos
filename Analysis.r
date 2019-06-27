#########################

# NOTES

# this R script contains the analyses presented in the paper 'Unleashing the BEAST: a brief measure for human social information use'
# Overall, this script follows the main text in terms of the order in which the analyses are presented. Where appropriate, this script will refer to the place where the results are presented (e.g. Figures and Tables).

# NB: this script assumes that it is executed in R in the order presented here; that is, later commands or functions might depend on previous variable definitions

# assumptions: have the packages 'lme4', 'rptR' and 'visreg' installed


##########################

# set path
#setwd('')

# read in the data for waves 1,2 and 3 from Experiment 1, and store in variable 'a'
# note that the variables for wave 2 and 3 have the respective suffices "_R" and "_RR"
a<-read.table('data_waves123.txt', sep='\t', header=TRUE)


# FIGURE 1E (distribution of mean adjustment S in Wave 1; Experiment 1)

#calculate adjustment 's' in each of the rounds
#the (somewhat clunky) variable names correspond to the animals shown in the stimuli (rounds 1-5: ants-bees-flamingos-cranes-crickets). Social info (in the paper referred to as 'X') has the name 'obsSocialInfo' with the round number appended to it.
a$s1<- (a$ants_revised - a$ants) / (a$obsSocialInfo1 - a$ants)
a$s2<- (a$bees_revised - a$bees) / (a$obsSocialInfo2 - a$bees)
a$s3<- (a$flamingos_revised - a$flamingos) / (a$obsSocialInfo3 - a$flamingos)
a$s4<- (a$cranes_revised - a$cranes) / (a$obsSocialInfo4 - a$cranes)
a$s5<- (a$crickets_revised - a$crickets) / (a$obsSocialInfo5 - a$crickets)

# remove (qualitatively different) cases that violated 0 <= s <= 1
# here, these columns are referred by their column number
for (i in 116:120) {
	a[,i]<-ifelse(a[,i] < 0, NA, a[,i])
	a[,i]<-ifelse(a[,i] > 1, NA, a[,i])
}

# calculate mean S (for wave 1, for all participants)
a$S<-NA;
for (i in 1:nrow(a)){
	a$S[i]<-mean(c(a$s1[i],a$s2[i],a$s3[i],a$s4[i],a$s5[i]), na.rm=T)
}

#calculate frequency distribution and store in vector 'f'
f<-rep(0,11)
for (k in 0:10){
	d<-length(which(round(a$S*10)==k))
	f[k+1]<-d
}
# normalise vector
f<-f/sum(f)

# plot to reproduce Figure 1e
par(cex.lab=1.5, cex.axis=1.5, las=1, yaxs='i')
plot(0, type='n', xlim=c(-0.05,1.05), ylim=c(0,0.25), xlab='', ylab='', axes=FALSE)
for (k in 0:10){
	rect(k/10-0.04, 0, k/10+0.04, f[k+1], col='#72a553')
}
axis(1)
axis(1, at=0:10/10)
axis(2)

# show distribution characteristics
summary(a$S)
# frequency of S=0
length(which(a$S==0))/length(a$S)
# frequency of 0<S<=0.5
(length(which(a$S>0)) - length(which(a$S>0.5)) ) / length(a$S)
# frequency of 0.5<S<=1
length(which(a$S>=0.5))/length(a$S)


##### INTERNAL CONSISTENCY OF RESPONSES IN THE BEAST #####

# internal correlations between rounds (reported in Table S1)
# again, columns are referred to by their number (these are the by-round adjustments ('s') calculated above)
internalCorMat<-matrix(nrow=5,ncol=5)
for (i in 1:5){
	for (j in 1:5){
		if (i<j) internalCorMat[i,j]<-cor.test(a[,115+i], a[,115+j])$estimate
		if (i>j) internalCorMat[i,j]<-cor.test(a[,115+i], a[,115+j])$p.value
	}
}
round(internalCorMat,3)

# another way to look at the same thing - Principal Component Analysis: variance summarized by Principal Components
pca1<-prcomp(~cbind(a$s1, a$s2, a$s3, a$s4, a$s5), na.action = na.omit)
summary(pca1)

# yet another way to look at the same thing, based on regression analysis, allowing to control for additional factors. Question here is: what proportion of variance is explained by the random effect of 'individual' in a regression fitted to adjustments?
# create a new matrix with each trial on a new row (for GLMM analysis)
mat<-matrix(nrow=0, ncol=8)
cntInd<-1; # counter
for (i in 1:nrow(a)){
	# for each decision store a unique identifier for 'participant',
	# the round number (i.e. the stimulus), the estimation error (how far was first estimate E1 off the true value T), the distance between E1 and social information X, and the adjustment s
	r<-c(cntInd, 1, a$ants[i], a$obsSocialInfo1[i], a$ants_revised[i], abs(a$ants[i]-93), abs(a$obsSocialInfo1[i] - a$ants[i]), a$s1[i])
	mat<-rbind(mat,r)
	r<-c(cntInd, 2, a$bees[i], a$obsSocialInfo2[i], a$bees_revised[i],abs(a$bees[i]-78), abs(a$obsSocialInfo2[i] - a$bees[i]), a$s2[i])
	mat<-rbind(mat,r)
	r<-c(cntInd, 3, a$flamingos[i], a$obsSocialInfo3[i], a$flamingos_revised[i],abs(a$flamingos[i]-59), abs(a$obsSocialInfo3[i] - a$flamingos[i]), a$s3[i])
	mat<-rbind(mat,r)
	r<-c(cntInd, 4, a$cranes[i], a$obsSocialInfo4[i], a$cranes_revised[i],abs(a$cranes[i]-74), abs(a$obsSocialInfo4[i] - a$cranes[i]), a$s4[i])
	mat<-rbind(mat,r)
	r<-c(cntInd, 5, a$crickets[i], a$obsSocialInfo5[i], a$crickets_revised[i],abs(a$crickets[i]-69), abs(a$obsSocialInfo5[i] - a$crickets[i]), a$s5[i])
	mat<-rbind(mat,r)
	cntInd<-cntInd+1;
}

mat<-data.frame(mat)
names(mat)<-c('ind', 'stimulus', 'E1',  'X', 'E2', 'estError', 'distSocInfo', 'adjustment');

# assume that estimates that were more than 100 off are a typo: remove that trial (not actually consequential for any of the main results, but this seems a reasonable assumption; consequences can be checked by just commenting out this line)
nrow(mat)
mat<-mat[mat$estError<100,]
nrow(mat)

# create factors from categorical variables
mat$stimulus<-factor(mat$stimulus)
mat$ind<-factor(mat$ind)

# first check repeatability of decisions to stay (s=0) or move (s>0; reported in Table S2)
mat$move<-ifelse(mat$adjustment==0,0,1)

# load libraries for regressions and repeatability analysis
library('rptR'); 
library('lme4');

# run analyses - exact outcomes may slightly vary with the optimizer used - all the ones I checked lead to very similar outcomes (see the manuals of rptR and lme4 for details).
m1<-glmer(move ~ estError + distSocInfo + stimulus + (1|ind), family='binomial', data=mat)
summary(m1)
rptBinary(move ~ estError + distSocInfo + stimulus + (1 | ind), grname="ind", data=mat, 
    nboot=10, npermut=10)

# then check repeatability of degrees of adjustment, conditional upon moving at all (s>0; reported in Table S3)
m1<-lmer(adjustment ~ estError + distSocInfo + stimulus + (1|ind), data=subset(mat, mat$adjustment>0))
summary(m1)
rpt(adjustment ~ estError + distSocInfo + stimulus + (1 | ind), grname = "ind", data=subset(mat, mat$adjustment>0), 
    nboot = 10, npermut = 10)

	
##### Accuracy of first and second estimates (Figure S3)
# normalise the estimates by dividing them over the true value in the trials
mat$E1n<-NA; mat$E2n<-NA;
trueValues<-c(93,78,59,74,69)
for (i in 1:nrow(mat)){
	for (j in 1:5) {
		if (mat$stimulus[i]==j) {
			mat$E1n[i]<-mat$E1[i]/trueValues[j]; mat$E2n[i]<-mat$E2[i]/trueValues[j]; 
		}
	}
}
head(mat)

# a simplified version of Figure S3 for Wave 1
par(mfrow=c(1,2))
hist(mat$E1n)
hist(mat$E2n)

# statistics on these distributions
summary(mat$E1n)
summary(mat$E2n)


# Does individuals' accuracy predict their mean adjustment (Figure S4)?
accMat<-matrix(nrow=0,ncol=2)
for (id in unique(mat$ind)){
	b<-subset(mat, mat$ind==id)
	meanError<-mean(b$estError, na.rm=TRUE)
	meanS<-mean(b$adjustment, na.rm=TRUE)
	accMat<-rbind(accMat, c(meanError, meanS))
}
accMat<-data.frame(accMat)
names(accMat)<-c('meanError','S')

# plot the relation
cor.test(accMat$meanError, accMat$S)
m1<-lm(S~meanError, data=accMat)
library('visreg')
visreg(m1,line=list(col='#72a553'), points=list(col=adjustcolor('#72a553',alpha=0), pch=16, cex=0.0), xlim=c(0, 50), ylim=c(0,1),xlab='', ylab='', axes=FALSE)
box()
axis(1, at=0:10*5, labels=FALSE)
axis(2, at=0:10/10, labels=FALSE)

# add points for individuals to the plot
for (i in 1:nrow(accMat)){
	points(accMat$meanError[i], accMat$S[i],pch=16, col='#72a553', cex=1)	
}

########### TIME CONSISTENCY OF SOCIAL INFORMATION USE IN THE BEAST (compare individuals' mean adjustment S across waves 1,2,3) #########

# calculate mean S for each of the Waves
a$s1_R<- (a$ants_revised_R - a$ants_R) / (a$obsSocialInfo1_R - a$ants_R)
a$s2_R<- (a$bees_revised_R - a$bees_R) / (a$obsSocialInfo2_R - a$bees_R)
a$s3_R<- (a$flamingos_revised_R - a$flamingos_R) / (a$obsSocialInfo3_R - a$flamingos_R)
a$s4_R<- (a$cranes_revised_R - a$cranes_R) / (a$obsSocialInfo4_R - a$cranes_R)
a$s5_R<- (a$crickets_revised_R - a$crickets_R) / (a$obsSocialInfo5_R - a$crickets_R)

a$s1_RR<- (a$ants_revised_RR - a$ants_RR) / (a$obsSocialInfo1_RR - a$ants_RR)
a$s2_RR<- (a$bees_revised_RR - a$bees_RR) / (a$obsSocialInfo2_RR - a$bees_RR)
a$s3_RR<- (a$flamingos_revised_RR - a$flamingos_RR) / (a$obsSocialInfo3_RR - a$flamingos_RR)
a$s4_RR<- (a$cranes_revised_RR - a$cranes_RR) / (a$obsSocialInfo4_RR - a$cranes_RR)
a$s5_RR<- (a$crickets_revised_RR - a$crickets_RR) / (a$obsSocialInfo5_RR - a$crickets_RR)

# as above remove the cases that violate 0<=s<=1
#remove rounds that violated 0 <= s <= 1
for (i in 122:131) {
	a[,i]<-ifelse(a[,i] < 0, NA, a[,i])
	a[,i]<-ifelse(a[,i] > 1, NA, a[,i])
}
# create a matrix 'Smat' that calculates S for each individuals across the three waves
Smat<-matrix(nrow=102, ncol=3)
for (i in 1:nrow(a)){
	a$S1[i]<-mean(c(a$s1[i],   a$s2[i],   a$s3[i],   a$s4[i],   a$s5[i]), na.rm=T)
	a$S2[i]<-mean(c(a$s1_R[i], a$s2_R[i], a$s3_R[i], a$s4_R[i], a$s5_R[i]), na.rm=T)
	a$S3[i]<-mean(c(a$s1_RR[i],a$s2_RR[i],a$s3_RR[i],a$s4_RR[i],a$s5_RR[i]), na.rm=T)
}
# calculate correlations
cor(a$S1, a$S2, use='complete', method='pearson')
cor.test(a$S1, a$S2)
cor(a$S1, a$S3, use='complete', method='pearson')
cor.test(a$S1, a$S3)
cor(a$S2, a$S3, use='complete', method='pearson')
cor.test(a$S2, a$S3)

# plot this (Figure 2)
par(mfrow=c(1,3), cex.lab=1.5, cex.axis=1.5, las=1, mar=c(4,7,1,1), xaxs='i', yaxs='i')
# wave 1 vs 2
m1<-lm(S2~S1, data=a)
visreg(m1, xlim=c(-0.05,1.05), ylim=c(-0.05,1.05), xlab='', ylab='', pch=16, line=list(col='#72a553'), points=list(col='#72a553', pch=16, cex=1), axes=FALSE)
axis(1)
axis(2)
box()

# wave 1 vs 3
m2<-lm(S3~S1, data=a)
visreg(m2, xlim=c(-0.05,1.05), ylim=c(-0.05,1.05), xlab='', ylab='', pch=16, line=list(col='#72a553'), points=list(col='#72a553', pch=16, cex=1), axes=FALSE)
axis(1)
axis(2)
box()

# wave 2 vs 3 [not shown in main paper]
m3<-lm(S3~S2, data=a)
visreg(m3, xlim=c(-0.05,1.05), ylim=c(-0.05,1.05), xlab='', ylab='', pch=16, line=list(col='#72a553'), points=list(col='#72a553', pch=16, cex=1), axes=FALSE)
axis(1)
axis(2)
box()


####### CROSS-TASK CONSISTENCY OF SOCIAL INFORMATION USE (Experiment 1) ########

# the data for this task is formatted somewhat differently (not each individual on one line, as the data set above, but rather trial-by-trial, with each trial on a new line)
# below we will slightly reorganise this in a simple way with some code (all commented as we go along).

# read the data from Experiment 2
a<-read.table('data_tasks123.txt', sep='\t', header=TRUE)

# create a new matrix 'indMat', storing summary statistics for each individual across the three tasks (BEAST, Moving Dots and Bandit task)
indMat<-matrix(ncol=11, nrow=0)

# loop over each each individual to calculate these summary statistics
for (ind in unique(a$uniquePlayerNr)){
	b<-subset(a, a$uniquePlayerNr==ind)
	
	##### BEAST #####
	# trails where the variable 'animalName' was defined are from the BEAST 
	d<-subset(b, b$animalName!='')
	d$s<-NA
	for (i in 1:5){
		d$s[i]<-(d$secondEstimate[i] - d$firstEstimate[i]) / (d$observedEstimate[i] - d$firstEstimate[i])
	}
	# remove overshoots and moving away
	d$s<-ifelse(d$s>1,NA,d$s)
	d$s<-ifelse(d$s<0,NA,d$s)
	meanS<-mean(d$s, na.rm=TRUE)

	# register which 'type' of the BEAST a participant played. These types are exactly the same, except for that in type 1, they entered their estimates in a numeric input field, and in type 2, they entered their estimates using a slider running from 1-150
	BEASTtype<-d$BEASTtype[1]

	##### MOVING DOTS #####
	# trials where the variable 'actualAngle' was defined are from the MOVING DOTS task
	d<-subset(b, !is.na(b$actualAngle))
	d$M<-NA
	for (i in 1:5){
		distanceMoved<- min(abs(d$secondAngle[i] - d$firstAngle[i]), abs(d$secondAngle[i]-360 - d$firstAngle[i]), abs(d$secondAngle[i]+360 - d$firstAngle[i]) )
		distancePeer<-min(abs(d$observedAngle[i] - d$firstAngle[i]), abs(d$observedAngle[i]-360 - d$firstAngle[i]), abs(d$observedAngle[i]+360 - d$firstAngle[i]) )
		d$M[i]<-distanceMoved/distancePeer
	}
	# remove overshoots and moving away
	d$M<-ifelse(d$M>1,NA,d$M)
	d$M<-ifelse(d$M<0,NA,d$M)
	y<-d$M
	# in cases where first angle matched the peer info (so division by zero occurred). this excludes that trial
	y<-y[!is.infinite(y)] 
	# calculate mean M
	meanM<-mean(y, na.rm=TRUE) 
	
	#### BANDIT TASK - calculate relative interest for social (rather than individual) information #####
	# trials where the variable 'infoChoice' was defined are from the BANDIT TASK
	d<-subset(b, !is.na(b$infoChoice))
	socChoice<-0;
	for (i in 1:14){
		if (d$infoChoice[i]==1) socChoice<-socChoice+1;
	}

	##### we also add age and gender to this matrix, as well as questionnaire scores (for analysis below, regarding Figure 4 of the main text)
		
	# questionnaire items
	d<-subset(b, b$age!='NA')
	IOS<-d$IOS[1]
	age<-d$age[1]
	gender<-d$female[1]
	
	# conformity
	x1<-which(names(d)=='conf1');
	conformity<-0;
	for (i in x1:(x1+11)){
		xx<-d[1,i];
		if (i==2 || i==7 || i==9 || i==11) xx<-10-xx # reversed items
		sc<-xx-5; # normalise so that 'neutral' is 0
		conformity<-conformity+sc
	}
	
	# individualism
	
	individualism<-d$indiv1[1] + d$indiv2[1] + d$indiv3[1] + d$indiv4[1]
	collectivism<- d$coll1[1] + d$coll2[1]+ d$coll3[1]+d$coll4[1]

	
	indRow<-c(ind, meanS, meanM, socChoice, age, gender, IOS, conformity, individualism, collectivism, BEASTtype)
	indMat<-rbind(indMat, indRow)
		
}

indMat<-data.frame(indMat)
# remove rownames for legibility and add titles to the columns
rownames(indMat)<-c()
names(indMat)<-c('participant', 'BEAST', 'MovingDots','Bandits', 'age', 'female', 'IOS', 'conformity', 'individualism', 'collectivism', 'BEASTtype')
indMat


# start plotting to create Figure 3, panels b and d (where data is shown)
par(mfrow=c(1,2),cex.lab=1.5, cex.axis=1.5, las=1, lend=1)

####### PLOT MOVING DOTS vs BEAST ########

# create cohorts and calculate means and SEs
f<-rep(0,11); se<-rep(0,11); n<-rep(0,11)
for (i in 0:10){
	d<-subset(indMat, round(indMat$MovingDots*10)==i)
	f[i+1]<-mean(d$BEAST)
	se[i+1]<-sd(d$BEAST)/sqrt(nrow(d))
	n[i+1]<-nrow(d)
}

m1<-lm(BEAST~MovingDots, data=indMat)
# calculate correlations
cor.test(indMat$MovingDots, indMat$BEAST)

# plot this, using visreg
visreg(m1,line=list(col='#72a553'), points=list(col=adjustcolor('#72a553',alpha=0), pch=16, cex=0.0), xlim=c(0, 1), ylim=c(0,1),xlab='', ylab='', axes=FALSE)
box()
axis(1, at=0:10/10, labels=FALSE)
axis(2, at=0:10/10, labels=FALSE)

# add by-cohort summary and individual points (with some jitter in the x direction) to the plot
for (i in 0:10){
	d<-subset(indMat, round(indMat$MovingDots*10)==i)
	if (nrow(d)>0){
		#add individual data points
		for (j in 1:nrow(d)){
			y<-d$BEAST[j];
			# add some jitter
			x<-i/10-0.3+runif(1)*0.6
			points(x,y,pch=16, col='#72a553', cex=0.3)
		}
		# add cohort mean + SE
		me<-mean(d$BEAST, na.rm=TRUE)
		se<-sd(d$BEAST, na.rm=TRUE)/sqrt(nrow(d))
		arrows(i/10,me-se, i/10, me+se, lwd=2, code=0)
		points(i/10, f[i+1],pch=15, col='#72a553', cex=1.3)	
	}
}


##### PLOT BANDITS TASK vs BEAST #######

# same approach as above, use cohorts (in this case just the counts of choosing social information)
f<-rep(0,14); se<-rep(0,14); n<-rep(0,14)
for (i in 0:14){
	d<-subset(indMat, indMat$Bandits==i)
	f[i+1]<-mean(d$BEAST, na.rm=TRUE)
	se[i+1]<-sd(d$BEAST, na.rm=TRUE)/sqrt(nrow(d))
	n[i+1]<-nrow(d)
}

cor.test(indMat$Bandits, indMat$BEAST)
m1<-lm(BEAST~Bandits, data=indMat)
visreg(m1,line=list(col='#72a553'), points=list(col=adjustcolor('#72a553',alpha=0), pch=16, cex=0.0), xlim=c(0.1,14.1), ylim=c(0,1),xlab='', ylab='', axes=FALSE)
box()
axis(1, at=0:14, labels=FALSE)
axis(2, at=0:10/10, labels=FALSE)

# add by-cohort summary and individual points (with some jitter in the x direction) to the plot
for (i in 0:14){
	d<-subset(indMat, indMat$Bandits==i)
	if (nrow(d)>0){
		#add individual data points
		for (j in 1:nrow(d)){
			y<-d$BEAST[j];
			x<-i-0.3+runif(1)*0.6
			points(x,y,pch=16, col='#72a553', cex=0.3)
		}
		# add cohort mean + SE
		me<-mean(d$BEAST, na.rm=TRUE)
		se<-sd(d$BEAST, na.rm=TRUE)/sqrt(nrow(d))
		arrows(i,me-se, i, me+se, lwd=2, code=0)
		points(i, f[i+1],pch=15, col='#72a553', cex=1.3)
	}
}


####### CONVERGENT VALIDITY (based on Experiment 1 (wave 1) and Experiment 2) ########
# In both Experiments, we had participants fill out questionnaires on the Inclusion of the Other in the Self (social proximity), Conformity, Individualism and Collectivism

# note that for Experiment 2, we already have all relevant data stored in the matrix 'indMat'
# for Experiment 1, we now create a similar matrix called 'indMat1'
# we then concatenate the two matrices to investigate the links between social information use in the BEAST and these questionnaire-based constructs

# first, read the data again (remember, the data of each participant is on 1 row)
a<-read.table('data_waves123.txt', sep='\t', header=TRUE)

# this repeats code from above (my apologies)
a$s1<- (a$ants_revised - a$ants) / (a$obsSocialInfo1 - a$ants)
a$s2<- (a$bees_revised - a$bees) / (a$obsSocialInfo2 - a$bees)
a$s3<- (a$flamingos_revised - a$flamingos) / (a$obsSocialInfo3 - a$flamingos)
a$s4<- (a$cranes_revised - a$cranes) / (a$obsSocialInfo4 - a$cranes)
a$s5<- (a$crickets_revised - a$crickets) / (a$obsSocialInfo5 - a$crickets)

# remove (qualitatively different) cases that violated 0 <= s <= 1
# here, these columns are referred by their column number
for (i in 116:120) {
	a[,i]<-ifelse(a[,i] < 0, NA, a[,i])
	a[,i]<-ifelse(a[,i] > 1, NA, a[,i])
}

# calculate mean S (for wave 1, for all participants)
a$S<-NA;
for (i in 1:nrow(a)){
	a$S[i]<-mean(c(a$s1[i],a$s2[i],a$s3[i],a$s4[i],a$s5[i]), na.rm=T)
}

# create matrix
indMat1<-matrix(ncol=11, nrow=0)
indCnt<-1;
for (i in 1:nrow(a)){

	# calculate conformity scale
	conformity<-
		a$conf1[i]-5 + 
		(10-a$conf2[i])-5 + 
		a$conf3[i]-5 + 
		a$conf4[i]-5 + 
		a$conf5[i]-5 + 
		a$conf6[i]-5 + 
		(10-a$conf7[i])-5 + 
		a$conf8[i]-5 + 
		(10-a$conf9[i])-5 + 
		a$conf10[i]-5 + 
		(10-a$conf11[i])-5

	individualism<-	a$indiv1[i] + a$indiv2[i] + a$indiv3[i] +a$indiv4[i]
	collectivism<-	a$coll1[i] + a$coll2[i] + a$coll3[i] + a$coll4[i]

	r<-c(indCnt, a$S[i], NA, NA, a$age[i], 2-a$gender[i], a$IOS[i], conformity, individualism, collectivism,1)
	indMat1<-rbind(indMat1, r)
	
	indCnt<-indCnt+1;
}

indMat1<-data.frame(indMat1)
# remove rownames for legibility and add titles to the columns
rownames(indMat1)<-c()
names(indMat1)<-c('participant', 'BEAST', 'MovingDots','Bandits', 'age', 'female', 'IOS', 'conformity', 'individualism', 'collectivism', 'BEASTtype')
indMat1

# Concatenate the matrices with individual-level summary statistics
mat<-rbind(indMat1, indMat)

# calculate and plot each of the correlations with the BEAST one by one (Figure 4)
dev.off()
dev.new(width=12.13, height=3.42)
par(mfrow=c(1,4), cex.lab=1.5, cex.axis=1.5, las=1, lend=1) 

# conformity
cor.test(mat$conformity, mat$BEAST)

m1<-lm(BEAST~conformity, data=mat)
visreg(m1,line=list(col='#72a553'), points=list(col=adjustcolor('#72a553',alpha=0), pch=16, cex=0.0), xlim=c(-34,42), ylim=c(0,1),xlab='', ylab='', axes=FALSE)
box()
axis(1, at=-5:5*10, labels=FALSE)
axis(2, at=0:10/10, labels=FALSE)

# add cohort means + individual data points
for (i in 1:10){
	xx<-(i-5);
	d<-subset(mat, round(mat$conformity / 10)==xx)
	# individual data points
	if (nrow(d)>0){
		for (j in 1:nrow(d)){
			y<-d$BEAST[j];
			x<-d$conformity[j]
			points(x,y,pch=16, col='#72a553', cex=0.3)
		}
		#cohort means
		me<-mean(d$BEAST)
		se<-sd(d$BEAST)/sqrt(nrow(d))
		arrows((i-5)*10,me-se, (i-5)*10, me+se, lwd=2, code=0)
		points((i-5)*10,me,pch=15, col='#72a553', cex=1.3)
	}
}


# inclusion of the other in the self
cor.test(mat$IOS, mat$BEAST)

m1<-lm(BEAST~IOS, data=mat)
visreg(m1,line=list(col='#72a553'), points=list(col=adjustcolor('#72a553',alpha=0), pch=16, cex=0.0), xlim=c(1,7), ylim=c(0,1),xlab='', ylab='', axes=FALSE)
box()
axis(1, at=1:7, labels=FALSE)
axis(2, at=0:10/10, labels=FALSE)

for (i in 1:7){
	d<-subset(mat, mat$IOS==i)
	# individual data points
	if (nrow(d)>0){
		for (j in 1:nrow(d)){
			y<-d$BEAST[j];
			x<-i-0.2+runif(1)*0.4
			points(x,y,pch=16, col='#72a553', cex=0.3)
		}
		# cohort means
		me<-mean(d$BEAST)
		se<-sd(d$BEAST)/sqrt(nrow(d))
		arrows(i,me-se, i, me+se, lwd=2, code=0)
		points(i,me,pch=15, col='#72a553', cex=1.3)
		
	}
}

# collectivism
cor.test(mat$collectivism, mat$BEAST)

m1<-lm(BEAST~collectivism, data=mat)
visreg(m1,line=list(col='#72a553'), points=list(col=adjustcolor('#72a553',alpha=0), pch=16, cex=0.0), xlim=c(3,36), ylim=c(0,1),xlab='', ylab='', axes=FALSE)
box()
axis(1, at=0:10*5, labels=FALSE)
axis(2, at=0:10/10, labels=FALSE)

for (i in 0:20){
	d<-subset(mat, round(mat$collectivism / 2)==i)
	if (nrow(d)>0){
		#individual data points
		for (j in 1:nrow(d)){
			y<-d$BEAST[j];
			x<-d$collectivism[j]
			points(x,y,pch=16, col='#72a553', cex=0.3)
		}
		# cohort means
		me<-mean(d$BEAST)
		se<-sd(d$BEAST)/sqrt(nrow(d))
		arrows(i*2,me-se, i*2, me+se, lwd=2, code=0)
		points(i*2,me,pch=15, col='#72a553', cex=1.3)
	}
}

# individualism
cor.test(mat$individualism, mat$BEAST)

m1<-lm(BEAST~individualism, data=mat)
visreg(m1,line=list(col='#72a553'), points=list(col=adjustcolor('#72a553',alpha=0), pch=16, cex=0.0), xlim=c(12,36), ylim=c(0,1),xlab='', ylab='', axes=FALSE)
box()
axis(1, at=0:10*5, labels=FALSE)
axis(2, at=0:10/10, labels=FALSE)

for (i in 0:20){
	d<-subset(mat, round(mat$individualism / 2)==i)
	if (nrow(d)>0){
		# individual data points
		for (j in 1:nrow(d)){		
			y<-d$BEAST[j];
			x<-d$individualism[j]
			points(x,y,pch=16, col='#72a553', cex=0.3)
			
		}
		
		# cohort means
		me<-mean(d$BEAST)
		se<-sd(d$BEAST)/sqrt(nrow(d))
		arrows(i*2,me-se, i*2, me+se, lwd=2, code=0)
		points(i*2,me,pch=15, col='#72a553', cex=1.3)
		
	}
}


#### for Table S4 (summarizing relations between BEAST and questionnaires)
# overall regression (Table S4)
m1<-lm(BEAST ~ conformity + IOS + individualism + collectivism, data=mat)
summary(m1)

# partial regressions
m1a<-lm(BEAST ~ conformity, data=mat)
summary(m1a)

m1b<-lm(BEAST ~ IOS, data=mat)
summary(m1b)

m1c<-lm(BEAST ~ individualism, data=mat)
summary(m1c)

m1d<-lm(BEAST ~ collectivism, data=mat)
summary(m1d)