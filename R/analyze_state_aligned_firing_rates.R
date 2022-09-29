library(lme4)
library(car)
library(emmeans)
library(PResiduals)

#subject<-'betta'
subject<-'samovar'
#n_states<-6
n_states<-5

# Read data
df<-read.csv(paste0(subject,'/state_aligned_firing_rates.csv'))
df$state<-as.factor(df$state)
df$day<-as.factor(df$day)
df$trial<-as.factor(df$trial)
df$electrode<-as.factor(df$electrode)
df$condition<-as.factor(df$condition)

n_interaction<-0
n_state<-0
n_condition<-0
electrodes<-unique(df$electrode)
state_changes<-c()
for(state in 1:n_states) {
  state_changes<-c(state_changes,0)
}
state_diffs<-matrix(0,n_states,n_states)
for(electrode in electrodes) {
  sub_df<-df[df$electrode==electrode & (df$condition=='baseline' | df$condition=='first'),]
  m<-lmer(rate ~ state*condition + (1 | day:trial), data=sub_df)
  a<-Anova(m)
  print(a)
  if(a$`Pr(>Chisq)`[3]<0.05) {
    pw<-emmeans(m,pairwise~condition|state)
    print(pw)
    contrasts<-as.data.frame(pw$contrasts)
    n_interaction<-n_interaction+1
    for(state in 1:n_states) {
      if(contrasts$p.value[contrasts$state==state]<0.05) {
        state_changes[state]<-state_changes[state]+1
      }
    }
  }
  else if(a$`Pr(>Chisq)`[1]<0.05) {
    n_state<-n_state+1
    pw<-emmeans(m,pairwise~state)
    print(pw)
    contrasts<-as.data.frame(pw$contrasts)
    for(s1 in 1:n_states) {
      for (s2 in 1:n_states) {
        if(s1!=s2) {
          if(contrasts$p.value[contrasts$contrast==paste0(s1,' - ',s2)]<0.05) {
            state_diffs[s1,s2]<-state_diffs[s1,s2]+1
          }
        }
      }
    }
  }
  else if(a$`Pr(>Chisq)`[2]<0.05) {
    n_condition<-n_condition+1
  }
}

n_interaction<-0
n_state<-0
n_condition<-0
electrodes<-unique(df$electrode)
state_changes<-c()
for(state in 1:n_states) {
  state_changes<-c(state_changes,0)
}
state_diffs<-matrix(0,n_states,n_states)
for(electrode in electrodes) {
  sub_df<-df[df$electrode==electrode & (df$condition=='second' | df$condition=='post'),]
  m<-lmer(rate ~ state*condition + (1 | day/trial), data=sub_df)
  a<-Anova(m)
  print(a)
  if(a$`Pr(>Chisq)`[3]<0.05) {
    pw<-emmeans(m,pairwise~condition|state)
    print(pw)
    contrasts<-as.data.frame(pw$contrasts)
    n_interaction<-n_interaction+1
    for(state in 1:n_states) {
      if(contrasts$p.value[contrasts$state==state]<0.05) {
        state_changes[state]<-state_changes[state]+1
      }
    }
  }
  else if(a$`Pr(>Chisq)`[1]<0.05) {
    n_state<-n_state+1
    pw<-emmeans(m,pairwise~state)
    print(pw)
    contrasts<-as.data.frame(pw$contrasts)
    for(s1 in 1:n_states) {
      for (s2 in 1:n_states) {
        if(s1!=s2) {
          if(contrasts$p.value[contrasts$contrast==paste0(s1,' - ',s2)]<0.05) {
            state_diffs[s1,s2]<-state_diffs[s1,s2]+1
          }
        }
      }
    }
  }
  else if(a$`Pr(>Chisq)`[2]<0.05) {
    n_condition<-n_condition+1
  }
}

n_interaction<-0
n_state<-0
n_condition<-0
electrodes<-unique(df$electrode)
state_changes<-c()
for(state in 1:n_states) {
  state_changes<-c(state_changes,0)
}
state_diffs<-matrix(0,n_states,n_states)
for(electrode in electrodes) {
  sub_df<-df[df$electrode==electrode & (df$condition=='first' | df$condition=='second'),]
  m<-lmer(rate ~ state*condition + (1 | day/trial), data=sub_df)
  a<-Anova(m)
  print(a)
  if(a$`Pr(>Chisq)`[3]<0.05) {
    pw<-emmeans(m,pairwise~condition|state)
    print(pw)
    contrasts<-as.data.frame(pw$contrasts)
    n_interaction<-n_interaction+1
    for(state in 1:n_states) {
      if(contrasts$p.value[contrasts$state==state]<0.05) {
        state_changes[state]<-state_changes[state]+1
      }
    }
  }
  else if(a$`Pr(>Chisq)`[1]<0.05) {
    n_state<-n_state+1
    pw<-emmeans(m,pairwise~state)
    print(pw)
    contrasts<-as.data.frame(pw$contrasts)
    for(s1 in 1:(n_states-1)) {
      for (s2 in (s1+1):n_states) {
        if(contrasts$p.value[contrasts$contrast==paste0(s1,' - ',s2)]<0.05) {
          state_diffs[s1,s2]<-state_diffs[s1,s2]+1
        }
      }
    }
  }
  else if(a$`Pr(>Chisq)`[2]<0.05) {
    n_condition<-n_condition+1
  }
}