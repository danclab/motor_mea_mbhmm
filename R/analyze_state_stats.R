library(lme4)
library(car)
library(emmeans)

subject<-'betta'
#subject<-'samovar'
n_states<-6
#n_states<-5

# Read data
df<-read.csv(paste0(subject,'/model_tv_plnorm_',n_states,'states_1_activations.csv'))
df$day<-as.factor(df$day)
df$trial<-as.factor(df$trial)
df$condition<-as.factor(df$condition)
df$state<-as.factor(df$state)
m1<-glmer(activations ~ condition*state + (1 | day/trial), family = "poisson", data=df)
print(Anova(m1))
print(emmeans(m1,pairwise~condition|state))

df<-read.csv(paste0(subject,'/model_tv_plnorm_',n_states,'states_1_lifetime.csv'))
df$trial<-as.factor(df$trial)
df$condition<-as.factor(df$condition)
df$state<-as.factor(df$state)
m2<-lmer(lifetime ~ condition*state + (1 | day/trial), data=df)
print(Anova(m2))
print(emmeans(m2,pairwise~condition|state))
