#------
pacman::p_load(pacman,R.matlab,car,reshape,dplyr,GGally,ggplot2,ggthemes,ggvis,
               httr,lubridate,plotly,rio,rmarkdown,shiny,stringr,tidyr,lme4,nlme,
               foreach,doParallel,fitdistrplus,DescTools,buildmer,pbkrtest,lmerTest,
               rstatix)


VPF_LMM_compare_contrasts <- function(df,subfield_levels){
  
  REML <- TRUE
  
  results_LME <- vector("list",length=N_subfields)
  names(results_LME) <- subfields_levels
  p_values <- list()
  F_value <- list()
  coeff <- list()
  err <- list()
  for (subfield in subfields_levels){
    inpdata <- df[df$subfield==subfield,]
    
    form <- "z~layer + contrast + layer:contrast + (1|subject)"
    
    m <- lmer(formula=form, data = inpdata, REML = TRUE,
              control = lmerControl(optimizer = "Nelder_Mead"))
    
    a <- anova(m,ddf="Kenward-Roger")
    s <- summary(m)
    
    p_values[[subfield]] <- list(p_contrast = a$'Pr(>F)'[2],
                                    p_interaction = a$'Pr(>F)'[3])
    
    F_value[[subfield]] <- list(F_contrast = a$'F value'[2],
                                   F_interaction = a$'F value'[3])
    
    coeff[[subfield]] <- list(coeff_contrast = s$coefficients[3,1],
                              coeff_interaction= s$coefficients[4,1])
    
    err[[subfield]] <- list(err_contrast = s$coefficients[3,2],
                            err_interaction= s$coefficients[4,2])
    
    
    results_LME[[subfield]] <- list(p_value = p_values[[subfield]], 
                                    F_value = F_value[[subfield]],
                                    coeff = coeff[[subfield]],
                                    err = err[[subfield]]) 
    
  }
  
  #correct for multiple comparisons
  for (ii in c(1,2)){
    tmp <- lapply(p_values,'[[',ii)
    tmp <- p.adjust(tmp, method="BH")
    
    #feed back corrected p-values
    for (subfield in subfields_levels) {
      results_LME[[subfield]][[1]][ii] <- tmp[subfield]
    }
  }
  
  return(results_LME)
}



s_pre <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/pre_vs_post_all_layers_masked.mat")

s_mem <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/memory_vs_math_all_layers_masked.mat")

#-------
subfields_levels <- c('Subiculum','CA1','CA2','CA3','CA4/DG')
depths <- c(1:30)
N_subfields <- length(subfields_levels)
N_subjects <- dim(s_pre[[1]])[1]/N_subfields


df_pre <- data.frame(matrix(s_pre[[1]], nrow = N_subjects*N_subfields, ncol=length(depths)))
df_mem <- data.frame(matrix(s_mem[[1]], nrow = N_subjects*N_subfields, ncol=length(depths)))
colnames(df_pre) <- depths
colnames(df_mem) <- depths

subjects <-c(1:N_subjects)
subjects <- rep(subjects,each=N_subfields)
subfields <- rep(subfields_levels,N_subjects)

df_pre$subfield <- subfields
df_pre$subject <- subjects
df_pre$contrast <- rep("pre_vs_post",N_subjects)
df_mem$subfield <- subfields
df_mem$subject <- subjects
df_mem$contrast <- rep("memory_vs_math",N_subjects)

# transform to long ------
df_pre_long <- melt(df_pre,id.vars=c("subfield","subject","contrast")) 
names(df_pre_long)[names(df_pre_long)=="variable"] <- "layer"
names(df_pre_long)[names(df_pre_long)=="value"] <- "z"

df_mem_long <- melt(df_mem,id.vars=c("subfield","subject","contrast")) 
names(df_mem_long)[names(df_mem_long)=="variable"] <- "layer"
names(df_mem_long)[names(df_mem_long)=="value"] <- "z"

#-------
df <-rbind(df_pre_long,df_mem_long)
#------
df$subfield <- factor(df$subfield,levels=subfields_levels)
df$subject <- factor(df$subject,levels=subjects[seq(1, length(subjects), by = N_subfields)])
df$layer <- as.numeric(df$layer)
df$contrast <- factor(df$contrast,levels=c("pre_vs_post","memory_vs_math"))

results_LME <- VPF_LMM_compare_contrasts(df,subfield_levels)


