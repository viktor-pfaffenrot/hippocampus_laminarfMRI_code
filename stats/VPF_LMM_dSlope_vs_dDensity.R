#----
pacman::p_load(pacman,R.matlab,car,reshape,dplyr,GGally,ggplot2,ggthemes,ggvis,
               httr,lubridate,plotly,rio,rmarkdown,shiny,stringr,tidyr,lme4,nlme,
               foreach,doParallel,fitdistrplus,DescTools,buildmer,pbkrtest,lmerTest,
               rstatix,performance,partR2)


data <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/breathhold/density_change_vs_slope_change.mat")

subfields_levels <- c('Subiculum','CA1','CA2','CA3','CA4/DG')
N_subfields <- length(subfields_levels)
N_subjects <- dim(data[[1]])[1]/N_subfields

df <- data.frame(matrix(data[[1]], nrow = N_subjects*N_subfields, ncol=2))
colnames(df) <- c("density_change","slope_change")

subjects <-c(1:N_subjects)
subjects <- rep(subjects,N_subfields)
subfields <- rep(subfields_levels,each=N_subjects)

df$subfield <- subfields
df$subject <- subjects

#----
df$subfield <-factor(df$subfield,levels=subfields_levels)
df$subject <- factor(df$subject,levels=subjects[1:N_subjects])
df$slope_change <- as.numeric(df$slope_change)
df$density_change <- as.numeric(df$density_change)
#----

m1 <- lmer(slope_change ~ density_change + (1|subfield) + (1|subject),data=df,REML = TRUE,
           control = lmerControl(optimizer = "Nelder_Mead"))
m2 <- lmer(slope_change ~ (1|subfield) + (1|subject),data=df,REML = TRUE,
           control = lmerControl(optimizer = "Nelder_Mead"))
browser()
m_comp_KR<-KRmodcomp(m1, m2)
model_summary <- summary(m1)
r2(m1)
#----

results <- list(KR = list(p_value = m_comp_KR$test$"p.value"[2], statistic = m_comp_KR$test$"stat"[2]),
                Estimates = model_summary$coefficients[,1],SE = model_summary$coefficients[,2])

outname <- "/media/pfaffenrot/Elements/postdoc/projects/data/avg/breathhold/results_stat_density_change_vs_slope_change.mat"

writeMat(outname,KR=results$KR,Estimates=results$Estimates,SE=results$SE)