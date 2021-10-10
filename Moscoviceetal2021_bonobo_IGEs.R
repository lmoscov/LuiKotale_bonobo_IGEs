# R scripts for the manuscript: “Dominance or Tolerance? Causes and consequences of a period of increased inter-community encounters among bonobos (Pan paniscus) at LuiKotale” 

# Authors: L. R. Moscovice, G. Hohmann, B. C. Trumble, B. Fruth & A. V. Jaeggi  

library(brms)
library(rethinking)

# Do FAI scores per plot differ between East and West? (Eco 1 prediction, see results section, Table 2)

FAI=read.table(file="bnb_fai_per_plot.txt",header=T, sep="\t",stringsAsFactors = T)

FAI$Habitat=relevel(FAI$Habitat,ref="west")

get_prior(bnb.fai~Habitat+(1|Plot), data=FAI, family=hurdle_gamma())

fai.model<- brm(bnb.fai~Habitat+(1|Plot), data=FAI, family=hurdle_gamma(), prior=
                  c(prior(normal(0,50), class=Intercept),
                    prior(normal(0,5), class=b),
                    prior(gamma(0.1,0.1), class=shape),
                    prior(beta(1,1), class=hu),
                    prior(normal(0,5), class=b, coef="Habitateast"),                                 prior(cauchy(0,2), class=sd)),
                sample_prior = TRUE, chains = 2, cores = 2, 
                iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99))


plot(fai.model)
summary(fai.model)

post.fai=posterior_samples(fai.model)

logistic(mean(post.fai$b_Intercept)); logistic(HPDI(post.fai$b_Intercept, prob=0.95)) 
# [1] 0.2392785
# |0.95      0.95| 
#  0.09694763 0.47150588 


logistic(mean(post.fai$b_Intercept+post.fai$b_Habitateast)); logistic(HPDI(post.fai$b_Intercept+post.fai$b_Habitateast, prob=0.95)) 
# [1] 0.4702691
# |0.95     0.95| 
#  0.2626439 0.6840419

sum(post.fai$b_Habitateast>0)/length(post.fai$b_Habitateast)
#[1] 0.9521667

#######################################

# Do number of bonobo trees per plot differ between East and West? (Eco 1 prediction, see Electronic Supplementary Materials (ESM Table 3)

trees=read.table(file="bnb_trees_per_plot.txt",header=T, sep="\t",stringsAsFactors = T)

trees$Habitat=relevel(trees$Habitat,ref="west")

get_prior(tot.bnb.trees~Habitat+(1|Plot), data=trees, family=poisson())

prior <- c(prior(normal(0,5), class = b), 
           prior(normal(0,50), class = Intercept),
           prior(cauchy(0,2), class=sd))

trees.model=brm(tot.bnb.trees~Habitat+(1|Plot), data=trees, family=poisson(), prior=prior, sample_prior = TRUE, chains = 2, cores = 2, iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99)) 

plot(trees.model)
summary(trees.model)

post.tree<-posterior_samples(trees.model)

logistic(mean(post.tree$b_Intercept)); logistic(HPDI(post.tree$b_Intercept, prob=0.95)) 
# [1] 0.7557946
# |0.95     0.95| 
#  0.7187034 0.7857594 


logistic(mean(post.tree$b_Intercept+post.tree$b_Habitateast)); logistic(HPDI(post.tree$b_Intercept+post.tree$b_Habitateast, prob=0.95)) 
# [1] 0.7484141
# |0.95     0.95| 
#  0.7126813 0.7828447

sum(post.tree$b_Habitateast>0)/length(post.tree$b_Habitateast)
# [1] 0.3781667

###############################################

# Does feeding on fleshy fruits differ between East and West (Eco 2 prediction, see results section, Table 3):

prop.fleshy=read.table(file="prop_scans_fleshy_fruit.txt", header=T, sep="\t",stringsAsFactors = T)

prop.fleshy$location=relevel(prop.fleshy$location,ref="West")

prop.fleshy$sex=relevel(prop.fleshy$sex,ref="M")

prop.fleshy$loc.code=as.numeric(prop.fleshy$location==levels(prop.fleshy$location)[2])

get_prior(num_fleshy_fruit_scans~location + sex + prop.dom.by.sex + offset(log(tot_food_scans)) + (1|subj)+ (1|date) + (0+loc.code|subj), data=prop.fleshy, family=zero_inflated_poisson("log"))

prop.fleshy.model=brm(num_fleshy_fruit_scans~ location + sex+ prop.dom.by.sex + offset(log(tot_food_scans)) + (1|subj)+ (1|date) + (0+loc.code|subj), data=prop.fleshy, family=zero_inflated_poisson("log"), prior=c(
  prior(normal(0,50), class=Intercept),
  prior(normal(0,5), class=b),
  prior(beta(1,1), class=zi),
  prior(cauchy(0,2), class=sd)),
  sample_prior = TRUE, chains = 2, cores = 2, 
  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99)) 

summary(prop.fleshy.model)

post.fleshy=posterior_samples(prop.fleshy.model)

# estimate for eating fleshy on West:

logistic(mean(post.fleshy$b_Intercept)); logistic(HPDI(post.fleshy$b_Intercept, prob=0.95)) 
# [1] 0.3036273
# |0.95     0.95| 
#  0.2563984 0.3582665 

# estimate for eating fleshy on East:

logistic(mean(post.fleshy$b_Intercept+post.fleshy$b_locationEast)); logistic(HPDI(post.fleshy$b_Intercept+post.fleshy$b_locationEast, prob=0.95)) 
# [1] 0.4345101
# |0.95     0.95| 
#  0.3783935 0.4901988 

# Odds Ratio for eating fleshy on west vs. east:

exp(mean(post.fleshy$b_locationEast)); exp(HPDI(post.fleshy$b_locationEast, prob=0.95)) 
# [1] 1.762284
# |0.95    0.95| 
#  1.451691 2.155324 

sum(post.fleshy$b_locationEast> 0)/length(post.fleshy$b_locationEast)
# [1] 1

##############################

# do West activity scans differ among baseline, incursion and IGE contexts? (ECO 3 and ICD/ICT 4, see results section Fig 2 and ESM Table 4)

d<- read.table(file="activity_scans.txt", header=T, sep="\t",stringsAsFactors = T) 

d$rank.z<- scale(d$prop.dom.by.sex) 
d$west.tot.c<- d$west.tot-mean(d$west.tot) 
d$hour.mn.prop.c<- d$hour.mn.prop-mean(d$hour.mn.prop)
d$location<- relevel(d$location, ref="West")
d$activity<- relevel(d$activity, ref="Rest")
d$sex=relevel(d$sex,ref="M")

get_prior(activity ~ location + ic.yn + rank.z + sex + west.tot.c + hour.mn.prop.c+(1+rank.z + west.tot.c + hour.mn.prop.c|subj), data = d, family = categorical())

prior <- c(prior(normal(0,5), class = b), 
           prior(normal(0,10), class = Intercept),
           prior(cauchy(0,2), class=sd, dpar=muFeed),
           prior(cauchy(0,2), class=sd, dpar=muMove),
           prior(cauchy(0,2), class=sd, dpar=muSocial),
           prior(lkj(2), class=cor)
)

m.time<- brm(activity ~ location + ic.yn + rank.z + sex + west.tot.c + hour.mn.prop.c+(1+rank.z + west.tot.c + hour.mn.prop.c|subj), data = d, family = categorical(), prior=prior, chains = 3, cores = 3, iter = 8000, warmup = 4000,control=list(adapt_delta= 0.99))

plot(m.time)
summary(m.time)
post.time<- posterior_samples(m.time)

## predicted probabilities of each behavior by west/east/IGE

location<- c("West", "East", "East")
newdata<- as.data.frame(location)
newdata$ic_yn<- c(0,0,1)
newdata$rank.z<- 0
newdata$west.tot.c<- 0
newdata$hour.mn.prop.c<- 0
newdata$sex<- "F"

pred.act<- fitted(m.time, newdata, re_formula=NA, summary=FALSE)

# means and HPDI's by context
rest.west.mean<- mean(pred.act[,1,1])
rest.east.mean<- mean(pred.act[,2,1])
rest.ige.mean<- mean(pred.act[,3,1])
rest.west.HPDI<- HPDI(pred.act[,1,1], prob=0.95)
rest.east.HPDI<- HPDI(pred.act[,2,1], prob=0.95)
rest.ige.HPDI<- HPDI(pred.act[,3,1], prob=0.95)

feed.west.mean<- mean(pred.act[,1,2])
feed.east.mean<- mean(pred.act[,2,2])
feed.ige.mean<- mean(pred.act[,3,2])
feed.west.HPDI<- HPDI(pred.act[,1,2], prob=0.95)
feed.east.HPDI<- HPDI(pred.act[,2,2], prob=0.95)
feed.ige.HPDI<- HPDI(pred.act[,3,2], prob=0.95)

move.west.mean<- mean(pred.act[,1,3])
move.east.mean<- mean(pred.act[,2,3])
move.ige.mean<- mean(pred.act[,3,3])
move.west.HPDI<- HPDI(pred.act[,1,3], prob=0.95)
move.east.HPDI<- HPDI(pred.act[,2,3], prob=0.95)
move.ige.HPDI<- HPDI(pred.act[,3,3], prob=0.95)

social.west.mean<- mean(pred.act[,1,4])
social.east.mean<- mean(pred.act[,2,4])
social.ige.mean<- mean(pred.act[,3,4])
social.west.HPDI<- HPDI(pred.act[,1,4], prob=0.95)
social.east.HPDI<- HPDI(pred.act[,2,4], prob=0.95)
social.ige.HPDI<- HPDI(pred.act[,3,4], prob=0.95)

## testing specific hypotheses 

# Feed during east < west?
sum(pred.act[,1,2]-pred.act[,2,2]>0)/length(pred.act[,2,2]) 
# 1

# Feed during ige > east? 
sum(pred.act[,3,2]-pred.act[,2,2]>0)/length(pred.act[,2,2]) 
# 0.79

# Move during east < west?
sum(pred.act[,1,3]-pred.act[,2,3]>0)/length(pred.act[,2,2]) 
# 0.0006

# Move during ige > east?
sum(pred.act[,3,3]-pred.act[,2,3]>0)/length(pred.act[,2,2]) 
# 0.34

# R during east > west?
sum(pred.act[,2,1]-pred.act[,1,1]>0)/length(pred.act[,2,2]) 
# 1

# R during ige < east?
sum(pred.act[,2,1]-pred.act[,3,1]>0)/length(pred.act[,2,2]) 
# 0.90

######################################

# Likelihood of maintaining close proximity to members of joining parties in the first hour following fusions when fusions involve in-group (base) vs. out-group (IGE) members (ICT 1, see results section, Table 5)

fuse=read.table(file="post_fuse_friendly.txt",header=T, sep="\t", fill=T,stringsAsFactors = T)

fuse$num.bef.z=as.vector(scale(fuse$num.bef))
fuse$num.join.z=as.vector(scale(fuse$num.join))

fuse$sex.res=factor(fuse$sex.res,levels=c("m","f"))

get_prior(friendly.yn~fusion.type*sex.res+dur.event+num.bef.z+num.join.z+(1|resident)+(1|fusion.id)+(0+fusion.type|resident)+(0+num.join.z|resident), data=fuse, family=bernoulli())

prior <- c(prior(normal(0,5), class = b), 
           prior(normal(0,50), class = Intercept),
           prior(cauchy(0,2), class=sd))

friendly.fuse.model_2way<- brm(friendly.yn~fusion.type*sex.res+dur.event+num.bef.z+num.join.z+(1|resident)+(1|fusion.id)+(0+fusion.type|resident)+(0+num.join.z|resident), data=fuse, family=bernoulli(), prior = prior, sample_prior = TRUE, chains = 2, cores = 2, iter = 8000, warmup = 4000, control = list(adapt_delta = 0.99))

twowayresults= summary(friendly.fuse.model_2way, waic = TRUE) # no influence of sex on maintenance of proximity with out-group vs in-group members

# re-run model without interaction term:

friendly.fuse.model<- brm(friendly.yn~fusion.type+sex.res+dur.event+num.bef.z+num.join.z+(1|resident)+(1|fusion.id)+(0+fusion.type|resident)+(0+num.join.z|resident), data=fuse, family=bernoulli(), prior = prior, sample_prior = TRUE, chains = 2, cores = 2, iter = 8000, warmup = 4000, control = list(adapt_delta = 0.99))

## how strong is the support for reduced friendly during IGE?

post.friendly<- posterior_samples(friendly.fuse.model)

sum(post.friendly$b_fusion.typeIGE<1)/length(post.friendly$b_fusion.typeIGE) 
# [1] 0.999

# likelihood of close proximity with in-group post-fusion

logistic(mean(post.friendly$b_Intercept)); logistic(HPDI(post.friendly$b_Intercept, prob=0.95)) 
# [1] 0.4143632
# |0.95      0.95| 
#  0.01997974 0.95707351 

# likelihood of close proximity with out-group post-fusion:

logistic(mean(post.friendly$b_Intercept+post.friendly$b_fusion.typeIGE)); logistic(HPDI(post.friendly$b_Intercept+post.friendly$b_fusion.typeIGE, prob=0.95)) 
# [1] 0.01285357
# |0.95      0.95| 
#  0.00048927 0.20838763 

#odds ratio

exp(mean(post.friendly$b_fusion.typeIGE)); exp(HPDI(post.friendly$b_fusion.typeIGE, prob=0.95)) 
# [1] 0.01840303
# |0.95        0.95| 
#  0.0008028031 0.3771642851 

# how many times less likely?

1/exp(mean(post.friendly$b_fusion.typeIGE))
# [1] 54.33889

# support for OR<1

sum(exp(post.friendly$b_fusion.typeIGE)<1)/length(post.friendly$b_fusion.typeIGE)
# [1] 0.994625

###########################################
