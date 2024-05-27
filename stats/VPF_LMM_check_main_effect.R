pacman::p_load(pacman,R.matlab,car,reshape,dplyr,GGally,ggplot2,ggthemes,ggvis,
               httr,lubridate,plotly,rio,rmarkdown,shiny,stringr,tidyr,lme4,nlme,
               foreach,doParallel,fitdistrplus,DescTools,buildmer,pbkrtest,lmerTest,
               rstatix,redres)

#when calling this script via Rscript, a bool is passed indicating whether pre_vs_post
#or memory_vs_math should be loaded
args <- commandArgs(trailingOnly = TRUE)
pre_vs_post <- tolower(args[1]) == "true"

#pre_vs_post = TRUE

# function to do the LME model comparison, i.e. to test whether there is a main
# effect of layer. If yes, this implies that the laminar profiles 'go in the same
# direction'
VPF_LMM_check_main_effect <- function(df_long,subfield_levels,pre_vs_post){
  
  REML <- TRUE
  
  results_LME <- vector("list",length=N_subfields)
  names(results_LME) <- subfields_levels
  p_values_KR <- list()
  F_value_KR <- list()
  coeff <- list()
  err <- list()
  for (subfield in subfields_levels){
    inpdata <- df_long[df_long$subfield==subfield,]
    
    #model with quadratic predictor
    m1_squared <- lmer(z~ layer + I(layer^2) +(1|subject), data = inpdata, REML = REML,
                control = lmerControl(optimizer = "Nelder_Mead"))      

    #model with layers
    m1_linear <- lmer(z~ layer +(1|subject), data = inpdata, REML = REML,
                             control = lmerControl(optimizer = "Nelder_Mead"))
    
    #simpler model, just with intercept
    m1_nolayer <- lmer(z~ 1 + (1|subject), data = inpdata, REML = REML,
                       control = lmerControl(optimizer = "Nelder_Mead"))
    
    #compare linear model with no predictor model and quadratic with linear model
    m_comp_KR <- list(lin_vs_no = KRmodcomp(m1_linear, m1_nolayer),
                      quad_vs_lin = KRmodcomp(m1_squared, m1_linear))
    
    s <- list(summary_lin = summary(m1_linear),
              summary_quad = summary(m1_squared))
    
    #some quality control metrics
    #raw_cond <- compute_redres(m1_linear)
    #hist(raw_cond)
    #plot_redres(m1_squared)+
    #  geom_smooth(method = "loess") +
    #  theme_classic()
    
    #bla <- inpdata %>% friedman_test(z ~ layer |subject)
    #if (subfield == "CA1"){browser()}
    
    
    p_values_KR[[subfield]] <- list(p_lin_vs_no = m_comp_KR$lin_vs_no$test$"p.value"[2],
                                    p_quad_vs_lin = m_comp_KR$quad_vs_lin$test$"p.value"[2])
    
    F_value_KR[[subfield]] <- list(F_lin_vs_no = m_comp_KR$lin_vs_no$test$"stat"[2],
                                   F_quad_vs_lin = m_comp_KR$quad_vs_lin$test$"stat"[2])
    coeff[[subfield]] <- list(coeff_lin = s$summary_lin$coefficients[,1],
                              coeff_quad= s$summary_quad$coefficients[,1])

    err[[subfield]] <- list(err_lin = s$summary_lin$coefficients[,2],
                              err_quad= s$summary_quad$coefficients[,2])
    
    
    results_LME[[subfield]] <- list(p_value = p_values_KR[[subfield]], 
                                    F_value = F_value_KR[[subfield]],
                                    coeff = coeff[[subfield]],
                                    err = err[[subfield]]) 
    
  }

  #correct for multiple comparisons
  for (ii in c(1,2)){
  tmp <- lapply(p_values_KR,'[[',ii)
  tmp <- p.adjust(tmp, method="BH")
  
  #feed back corrected p-values
  for (subfield in subfields_levels) {
      results_LME[[subfield]][[1]][ii] <- tmp[subfield]
      }
  }
  
  return(results_LME)
}

if (pre_vs_post==TRUE){
  data <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/pre_vs_post_all_layers_masked.mat")
} else {
  data <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/memory_vs_math_all_layers_masked.mat")
}
subfields_levels <- c('Subiculum','CA1','CA2','CA3','CA4/DG')
depths <- c(1:30)
N_subfields <- length(subfields_levels)
N_subjects <- dim(data[[1]])[1]/N_subfields

df <- data.frame(matrix(data[[1]], nrow = N_subjects*N_subfields, ncol=length(depths)))
colnames(df) <- depths

subjects <-c(1:N_subjects)
subjects <- rep(subjects,each=N_subfields)
subfields <- rep(subfields_levels,N_subjects)

df$subfield <- subfields
df$subject <- subjects

# transform to long
df_long <- melt(df,id.vars=c("subfield","subject")) 
names(df_long)[names(df_long)=="variable"] <- "layer"
names(df_long)[names(df_long)=="value"] <- "z"

# create factors
df_long$subfield <- factor(df_long$subfield,levels=subfields_levels)
df_long$subject <- factor(df_long$subject,levels=subjects[seq(1, length(subjects), by = N_subfields)])
df_long$layer <- as.numeric(df_long$layer)
#----
results <- VPF_LMM_check_main_effect(df_long,subfield_levels,pre_vs_post)
browser()
if (pre_vs_post==TRUE){
  outname <- "/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/res_LMM_check_main_effect_pre_vs_post.mat"
} else {
  outname <- "/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/res_LMM_check_main_effect_memory_vs_math.mat"
}
writeMat(outname,Subiculum_p=results$Subiculum$p_value,Subiculum_F=results$Subiculum$F_value,Subiculum_coeff=results$Subiculum$coeff,Subiculum_err = results$Subiculum$err, 
                 CA1_p=results$CA1$p_value,CA1_F=results$CA1$F_value,CA1_coeff=results$CA1$coeff,CA1_err = results$CA1$err, 
                 CA2_p=results$CA2$p_value,CA2_F=results$CA2$F_value,CA2_coeff=results$CA2$coeff,CA2_err = results$CA2$err, 
                 CA3_p=results$CA3$p_value,CA3_F=results$CA3$F_value,CA3_coeff=results$CA3$coeff,CA3_err = results$CA3$err, 
                 CA4_p=results$"CA4/DG"$p_value,CA4_F=results$"CA4/DG"$F_value,CA4_coeff=results$"CA4/DG"$coeff,CA4_err = results$"CA4/DG"$err)
