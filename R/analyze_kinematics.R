library(lme4)
library(car)
library(emmeans)
library(PResiduals)

#subject<-'betta'
subject<-'samovar'
#n_states<-6
n_states<-5

# Read data
df<-read.csv(paste0(subject,'/kinematic_data.csv'))
df$day<-as.factor(df$day)
df$trial<-as.factor(df$trial)
df$condition<-as.factor(df$condition)

m1<-lmer(rt ~ condition + (1 | day), data=df)
Anova(m1)

m2<-lmer(mt ~ condition + (1 | day), data=df)
Anova(m2)
emmeans(m2,pairwise~condition)


m3<-lmer(pt ~ condition + (1 | day), data=df)
Anova(m3)
emmeans(m3,pairwise~condition)
