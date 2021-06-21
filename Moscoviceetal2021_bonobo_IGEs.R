# R scripts for the manuscript: “Dominance or Tolerance? Causes and consequences of a period of increased inter-community encounters among bonobos (Pan paniscus) at LuiKotale” 

# Authors: L. R. Moscovice, G. Hohmann, B. C. Trumble, B. Fruth & A. V. Jaeggi  

library(brms)
library(rethinking)

# Do FAI scores per plot differ between East and West? (Eco 1 prediction, see results section, Table 3)

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

# Does feeding on fleshy fruits differ between East and West (Eco 2 prediction, see results section, Table 4):

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

# do West activity scans differ among baseline, incursion and IGE contexts? (ECO 3 and ICD/ICT 3, see results section Fig 2 and ESM Table 4)

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

 
# Likelihood of maintaining close proximity to members of joining parties in the first hour following fusions when fusions involve in-group (base) vs. out-group (IGE) members (ICD/ICT 1, see results section, Table 5)

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

# Analysis of log cortisol concentrations (ng/ml SG) by context (ICD/ICT 5, see results section, Fig 4 and Table 6)

cort=read.table(file="cort_responses.txt", header=T, sep="\t", stringsAsFactors = T) 

get_prior(cort.log~context-1+time.z + (1|subj)+(1|date), data=cort, family=gaussian())

prior <- c(prior(normal(0,1), class = b), 
           prior(cauchy(0,2), class=sd),
           prior(cauchy(0,2), class=sigma))

cort.mod<-brm(cort.log~context-1+time.z + (1|subj)+(1|date), data=cort, family=gaussian(), prior=prior, sample_prior = TRUE, chains = 3, cores = 3, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))

plot(cort.mod)
summary(cort.mod)

#############################################################