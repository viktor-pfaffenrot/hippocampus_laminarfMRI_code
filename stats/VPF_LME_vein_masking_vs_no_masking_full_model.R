#------
pacman::p_load(pacman,R.matlab,car,reshape,dplyr,GGally,ggplot2,ggthemes,ggvis,
               httr,lubridate,plotly,rio,rmarkdown,shiny,stringr,tidyr,lme4,nlme,
               foreach,doParallel,fitdistrplus,DescTools,buildmer,pbkrtest,lmerTest,
               rstatix)

pre_vs_post <- TRUE

VPF_LMM_compare_masking <- function(df,subfield_levels){
  
  REML <- TRUE
  
  results_LME <- vector("list",length=N_subfields)
  names(results_LME) <- subfields_levels
  p_values <- list()
  F_value <- list()
  coeff <- list()
  err <- list()
  for (subfield in subfields_levels){
    inpdata <- df[df$subfield==subfield,]
    
    form <- "z~layer + masking + layer:masking + (1|subject)"
    
    m <- lmer(formula=form, data = inpdata, REML = TRUE,
              control = lmerControl(optimizer = "Nelder_Mead"))
    
    a <- anova(m,ddf="Kenward-Roger")
    s <- summary(m)
    
    p_values[[subfield]] <- list(p_masking = a$'Pr(>F)'[2],
                                 p_interaction = a$'Pr(>F)'[3])
    
    F_value[[subfield]] <- list(F_masking = a$'F value'[2],
                                F_interaction = a$'F value'[3])
    
    coeff[[subfield]] <- list(coeff_masking = s$coefficients[3,1],
                              coeff_interaction= s$coefficients[4,1])
    
    err[[subfield]] <- list(err_masking = s$coefficients[3,2],
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



# read matrix and create data frame  ------------------------------------------
if (pre_vs_post==TRUE) {
  s_masked <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/pre_vs_post_all_layers_masked.mat")
  s_no_masked <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/pre_vs_post_all_layers_no_masked.mat")
  
} else {
  s_masked <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/memory_vs_math_all_layers_masked.mat")
  s_no_masked <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/memory_vs_math_all_layers_no_masked.mat")
  
}

subfields_levels <- c('Subiculum','CA1','CA2','CA3','CA4/DG')
depths <- c(1:30)
masking <- c("masked","not_masked")
N_subfields <- length(subfields_levels)
N_subjects <- dim(s_masked[[1]])[1]/N_subfields

df_masked    <- data.frame(matrix(s_masked[[1]], nrow = N_subjects*N_subfields, ncol=length(depths)))
df_no_masked <- data.frame(matrix(s_no_masked[[1]], nrow = N_subjects*N_subfields, ncol=length(depths)))
colnames(df_masked) <- depths
colnames(df_no_masked) <- depths

subjects <-c(1:N_subjects)
subjects <- rep(subjects,each=N_subfields)
subfields <- rep(subfields_levels,N_subjects)

df_masked$subfield <- subfields
df_masked$subject <- subjects
df_masked$masking <- rep("masked",N_subjects)

df_no_masked$subfield <- subfields
df_no_masked$subject <- subjects
df_no_masked$masking <- rep("not_masked",N_subjects)



# transform to long ------
df_masked_long <- melt(df_masked,id.vars=c("subfield","subject","masking")) 
names(df_masked_long)[names(df_masked_long)=="variable"] <- "layer"
names(df_masked_long)[names(df_masked_long)=="value"] <- "z"

df_no_masked_long <- melt(df_no_masked,id.vars=c("subfield","subject","masking")) 
names(df_no_masked_long)[names(df_no_masked_long)=="variable"] <- "layer"
names(df_no_masked_long)[names(df_no_masked_long)=="value"] <- "z"

df_long <- rbind(df_masked_long,df_no_masked_long)
# create factors --------
df_long$subfield <- factor(df_long$subfield,levels=subfields_levels)
df_long$subject <- factor(df_long$subject,levels=subjects[seq(1, length(subjects), by = N_subfields)])
df_long$layer <- as.numeric(df_long$layer)
df_long$masking <-factor(df_long$masking,level=masking)


results_LME_memory_vs_math <-VPF_LMM_compare_masking(df_long,subfields_levels)

