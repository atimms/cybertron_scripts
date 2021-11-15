##using r4.0.4
library(ggplot2)
library(tidyr)

setwd('/archive/mirzaa_g/misc/PPP2R5D_graphs_1021')

##initial data

##read in data for cbcl
data_wide <- read.table('cbcl_data_100521.txt', header=T, sep='\t')

##convert to long format
#http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/
data_long <- gather(data_wide, score, t_score, ach_anxdep_t_score:dsm5_conduct_t_score, factor_key=TRUE)
data_long

##graph as dotplot
#http://www.sthda.com/english/wiki/ggplot2-dot-plot-quick-start-guide-r-software-and-data-visualization
dp <-ggplot(data_long, aes(x=score, y=t_score, fill=score, binwidth = 1)) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 1,show.legend = FALSE)+
  labs(y = "T score")
dp + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("cbcl_dotplot.100521.pdf", width = 9, height = 6)


##bar charts for counts
pheno_count_data <- read.table('pheno_data_100521.txt', header=T, sep='\t')
head(pheno_count_data)
##to have same order as input df... https://stackoverflow.com/questions/12774210/how-do-you-specifically-order-ggplot2-x-axis-instead-of-alphabetical-order
#Turn your 'treatment' column into a character vector
pheno_count_data$Condition <- as.character(pheno_count_data$Condition)
#Then turn it back into a factor with the levels in the correct order
pheno_count_data$Condition <- factor(pheno_count_data$Condition, levels=unique(pheno_count_data$Condition))
##graph
ggplot(data=pheno_count_data, aes(x=Condition, y=count, fill=type)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") + theme_classic()  + theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("pheno_barchart.100521.pdf", width = 18, height = 12)

##second set of data .... 1014

##CBCL6-18_T
##read in data for cbcl
data_wide <- read.table('CBCL6-18_T.txt', header=T, sep='\t')
##convert to long format
keycol <- "score"
valuecol <- "t_score"
gathercols <- c('Anxious_Depressed', 'Withdrawn', 'Sleep_Problems', 'Social_Problems', 'Thought_Problems', 'Attention_Problems', 'Rulebreaking_Behavior', 'Aggressive_Behavior', 'Internalizing_Problems', 'Externalizing_Problems', 'Total_Problems', 'DSM5_Depressive_Problems', 'DSM5_Anxiety_Problems', 'DSM5_Somatic_Problems', 'Attention_Deficit_Hyperactivity_Problems', 'DSM5_Oppositional_Defiant_Problems', 'DSM5_Conduct_Problems')
data_long <- gather_(data_wide, keycol, valuecol, gathercols)
data_long
##have same order as input
#Turn your 'treatment' column into a character vector
data_long$score <- as.character(data_long$score)
#Then turn it back into a factor with the levels in the correct order
data_long$score <- factor(data_long$score, levels=unique(data_long$score))
##graph as dotplot, no color i.e. removed fill=score
dp <-ggplot(data_long, aes(x=score, y=t_score, binwidth = 1)) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 1,show.legend = FALSE)+
  labs(y = "T score")
dp + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1)) + geom_hline(yintercept=70) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.5, color = 'red') +
  scale_x_discrete(breaks=c('Anxious_Depressed', 'Withdrawn', 'Sleep_Problems', 'Social_Problems', 'Thought_Problems', 'Attention_Problems', 'Rulebreaking_Behavior', 'Aggressive_Behavior', 'Internalizing_Problems', 'Externalizing_Problems', 'Total_Problems', 'DSM5_Depressive_Problems', 'DSM5_Anxiety_Problems', 'DSM5_Somatic_Problems', 'Attention_Deficit_Hyperactivity_Problems', 'DSM5_Oppositional_Defiant_Problems', 'DSM5_Conduct_Problems'),
                     labels=c('Anxious/Depressed', 'Withdrawn', 'Sleep Problems', 'Social Problems', 'Thought Problems', 'Attention Problems', 'Rulebreaking Behavior', 'Aggressive Behavior', 'Internalizing Problems', 'Externalizing Problems', 'Total Problems', 'DSM5 Depressive Problems', 'DSM5 Anxiety Problems', 'DSM5 Somatic Problems', 'Attention Deficit/Hyperactivity Problems', 'DSM5 Oppositional Defiant Problems', 'DSM5 Conduct Problems'))
ggsave("CBCL6-18_T.101421.pdf", width = 9, height = 6)

##CBCL6-18_percentile
data_wide <- read.table('CBCL6-18_percentile.txt', header=T, sep='\t')
##convert to long format
keycol <- "score"
valuecol <- "t_score"
gathercols <- c('Anxious_Depressed', 'Withdrawn', 'Sleep_Problems', 'Social_Problems', 'Thought_Problems', 'Attention_Problems', 'Rulebreaking_Behavior', 'Aggressive_Behavior', 'Internalizing_Problems', 'Externalizing_Problems', 'Total_Problems', 'DSM5_Depressive_Problems', 'DSM5_Anxiety_Problems', 'DSM5_Somatic_Problems', 'Attention_Deficit_Hyperactivity_Problems', 'DSM5_Oppositional_Defiant_Problems', 'DSM5_Conduct_Problems')
data_long <- gather_(data_wide, keycol, valuecol, gathercols)
data_long
##have same order as input
#Turn your 'treatment' column into a character vector
data_long$score <- as.character(data_long$score)
#Then turn it back into a factor with the levels in the correct order
data_long$score <- factor(data_long$score, levels=unique(data_long$score))
##graph as dotplot, no color i.e. removed fill=score
dp <-ggplot(data_long, aes(x=score, y=t_score, binwidth = 1)) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 1,show.legend = FALSE)+
  labs(y = "percentile")
dp + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1)) + geom_hline(yintercept=98)+
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.5, color = 'red')+
  scale_x_discrete(breaks=c('Anxious_Depressed', 'Withdrawn', 'Sleep_Problems', 'Social_Problems', 'Thought_Problems', 'Attention_Problems', 'Rulebreaking_Behavior', 'Aggressive_Behavior', 'Internalizing_Problems', 'Externalizing_Problems', 'Total_Problems', 'DSM5_Depressive_Problems', 'DSM5_Anxiety_Problems', 'DSM5_Somatic_Problems', 'Attention_Deficit_Hyperactivity_Problems', 'DSM5_Oppositional_Defiant_Problems', 'DSM5_Conduct_Problems'),
                     labels=c('Anxious/Depressed', 'Withdrawn', 'Sleep Problems', 'Social Problems', 'Thought Problems', 'Attention Problems', 'Rulebreaking Behavior', 'Aggressive Behavior', 'Internalizing Problems', 'Externalizing Problems', 'Total Problems', 'DSM5 Depressive Problems', 'DSM5 Anxiety Problems', 'DSM5 Somatic Problems', 'Attention Deficit/Hyperactivity Problems', 'DSM5 Oppositional Defiant Problems', 'DSM5 Conduct Problems'))
ggsave("CBCL6-18_percentile.101421.pdf", width = 9, height = 6)

##CBCL2-5_T
##read in data for cbcl
data_wide <- read.table('CBCL2-5_T.txt', header=T, sep='\t')
##convert to long format
keycol <- "score"
valuecol <- "t_score"
gathercols <- c('Emotionally_Reactive', 'Anxious_Depressed', 'Somatic_Compaints', 'Withdrawn', 'Sleep_Problems', 'Attention_Problems', 'Aggressive_Behavior', 'Internalizing_Problems_', 'Externalizing_Problems', 'Total_Problems', 'Stress_Problems', 'DSM5_Depressive_Problems', 'DSM5_Anxiety_Problems', 'DSM5_Autism_Spectrum_Problems', 'DSM5_Attention_Deficit_Hyperactivity_Problems', 'DSM5_Oppositional_Defiant_Problems')
data_long <- gather_(data_wide, keycol, valuecol, gathercols)
data_long
##have same order as input
#Turn your 'treatment' column into a character vector
data_long$score <- as.character(data_long$score)
#Then turn it back into a factor with the levels in the correct order
data_long$score <- factor(data_long$score, levels=unique(data_long$score))
##graph as dotplot, no color i.e. removed fill=score
dp <-ggplot(data_long, aes(x=score, y=t_score, binwidth = 1)) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 1,show.legend = FALSE)+
  labs(y = "T score")
dp + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1)) + geom_hline(yintercept=70)+
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.5, color = 'red')+
scale_x_discrete(breaks=c('Emotionally_Reactive', 'Anxious_Depressed', 'Somatic_Compaints', 'Withdrawn', 'Sleep_Problems', 'Attention_Problems', 'Aggressive_Behavior', 'Internalizing_Problems_', 'Externalizing_Problems', 'Total_Problems', 'Stress_Problems', 'DSM5_Depressive_Problems', 'DSM5_Anxiety_Problems', 'DSM5_Autism_Spectrum_Problems', 'DSM5_Attention_Deficit_Hyperactivity_Problems', 'DSM5_Oppositional_Defiant_Problems'),
                 labels=c('Emotionally Reactive', 'Anxious/Depressed', 'Somatic Compaints', 'Withdrawn', 'Sleep Problems', 'Attention Problems', 'Aggressive Behavior', 'Internalizing Problems ', 'Externalizing Problems', 'Total Problems', 'Stress Problems', 'DSM5 Depressive Problems', 'DSM5 Anxiety Problems', 'DSM5 Autism Spectrum Problems', 'DSM5 Attention Deficit/Hyperactivity Problems', 'DSM5 Oppositional Defiant Problems'))

ggsave("CBCL2-5_T.101421.pdf", width = 9, height = 6)

##CBCL2-5_percentile
data_wide <- read.table('CBCL2-5_percentile.txt', header=T, sep='\t')
##convert to long format
keycol <- "score"
valuecol <- "t_score"
gathercols <- c('Emotionally_Reactive', 'Anxious_Depressed', 'Somatic_Compaints', 'Withdrawn', 'Sleep_Problems', 'Attention_Problems', 'Aggressive_Behavior', 'Internalizing_Problems_', 'Externalizing_Problems', 'Total_Problems', 'Stress_Problems', 'DSM5_Depressive_Problems', 'DSM5_Anxiety_Problems', 'DSM5_Autism_Spectrum_Problems', 'DSM5_Attention_Deficit_Hyperactivity_Problems', 'DSM5_Oppositional_Defiant_Problems')
data_long <- gather_(data_wide, keycol, valuecol, gathercols)
data_long
##have same order as input
#Turn your 'treatment' column into a character vector
data_long$score <- as.character(data_long$score)
#Then turn it back into a factor with the levels in the correct order
data_long$score <- factor(data_long$score, levels=unique(data_long$score))
##graph as dotplot, no color i.e. removed fill=score
dp <-ggplot(data_long, aes(x=score, y=t_score, binwidth = 1)) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 1,show.legend = FALSE)+
  labs(y = "percentile")
dp + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1)) + geom_hline(yintercept=98)+
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.5, color = 'red') +
  scale_x_discrete(breaks=c('Emotionally_Reactive', 'Anxious_Depressed', 'Somatic_Compaints', 'Withdrawn', 'Sleep_Problems', 'Attention_Problems', 'Aggressive_Behavior', 'Internalizing_Problems_', 'Externalizing_Problems', 'Total_Problems', 'Stress_Problems', 'DSM5_Depressive_Problems', 'DSM5_Anxiety_Problems', 'DSM5_Autism_Spectrum_Problems', 'DSM5_Attention_Deficit_Hyperactivity_Problems', 'DSM5_Oppositional_Defiant_Problems'),
                     labels=c('Emotionally Reactive', 'Anxious/Depressed', 'Somatic Compaints', 'Withdrawn', 'Sleep Problems', 'Attention Problems', 'Aggressive Behavior', 'Internalizing Problems ', 'Externalizing Problems', 'Total Problems', 'Stress Problems', 'DSM5 Depressive Problems', 'DSM5 Anxiety Problems', 'DSM5 Autism Spectrum Problems', 'DSM5 Attention Deficit/Hyperactivity Problems', 'DSM5 Oppositional Defiant Problems'))
ggsave("CBCL2-5_percentile.101421.pdf", width = 9, height = 6)

##vineland
data_wide <- read.table('vineland.txt', header=T, sep='\t')
##convert to long format
keycol <- "score_type"
valuecol <- "score"
gathercols <- c('Communication_domain', 'Daily_living_skills_domain', 'Socialization_domain', 'Motor_skills_domain', 'Adaptive_Behavior_Composite')
data_long <- gather_(data_wide, keycol, valuecol, gathercols)
data_long
##have same order as input
#Turn your 'treatment' column into a character vector
data_long$score_type <- as.character(data_long$score_type)
#Then turn it back into a factor with the levels in the correct order
data_long$score_type <- factor(data_long$score_type, levels=unique(data_long$score_type))
##graph as dotplot, no color i.e. removed fill=score_type
dp <-ggplot(data_long, aes(x=score_type, y=score, binwidth = 1)) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 1,show.legend = FALSE)+
  labs(y = "score")
dp + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1)) + geom_hline(yintercept = c(70,85,114))+ expand_limits(y=c(0,120)) +
  scale_x_discrete(breaks=c('Communication_domain', 'Daily_living_skills_domain', 'Socialization_domain', 'Motor_skills_domain', 'Adaptive_Behavior_Composite'),
                   labels=c('Communication domain', 'Daily living skills domain', 'Socialization domain', 'Motor skills domain', 'Adaptive Behavior Composite'))
ggsave("vineland.101421.pdf", width = 9, height = 6)

##bar charts for counts
pheno_count_data <- read.table('pheno_data_101521.txt', header=T, sep='\t')
head(pheno_count_data)
##to have same order as input df... https://stackoverflow.com/questions/12774210/how-do-you-specifically-order-ggplot2-x-axis-instead-of-alphabetical-order
#Turn your 'treatment' column into a character vector
pheno_count_data$Condition <- as.character(pheno_count_data$Condition)
#Then turn it back into a factor with the levels in the correct order
pheno_count_data$Condition <- factor(pheno_count_data$Condition, levels=unique(pheno_count_data$Condition))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
##graph
ggplot(data=pheno_count_data, aes(x=Condition, y=count, fill=type)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") + theme_classic()  + theme(axis.text.x = element_text(angle = 45, hjust=1)) +  scale_fill_manual(values=cbPalette)
ggsave("pheno_barchart.101521.pdf", width = 18, height = 12)

##bar charts for counts
pheno_count_data <- read.table('pheno_data_101821.txt', header=T, sep='\t')
head(pheno_count_data)
##to have same order as input df... https://stackoverflow.com/questions/12774210/how-do-you-specifically-order-ggplot2-x-axis-instead-of-alphabetical-order
#Turn your 'treatment' column into a character vector
pheno_count_data$Condition <- as.character(pheno_count_data$Condition)
#Then turn it back into a factor with the levels in the correct order
pheno_count_data$Condition <- factor(pheno_count_data$Condition, levels=unique(pheno_count_data$Condition))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
##graph
ggplot(data=pheno_count_data, aes(x=Condition, y=count, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") + theme_classic()  + theme(axis.text.x = element_text(angle = 45, hjust=1)) +  scale_fill_manual(values=cbPalette)
ggsave("pheno_barchart.101821.pdf", width = 18, height = 12)

