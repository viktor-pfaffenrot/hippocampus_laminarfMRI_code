pacman::p_load(pacman,R.matlab,car,reshape,dplyr,GGally,ggplot2,ggthemes,ggvis,
               httr,lubridate,plotly,rio,rmarkdown,shiny,stringr,tidyr,lme4,nlme,
               flexplot,foreach,doParallel,fitdistrplus,DescTools,buildmer,pbkrtest,lmerTest,
               rstatix)

#when calling this script via Rscript, a bool is passed indicating whether pre_vs_post
#or memory_vs_math should be loaded
args <- commandArgs(trailingOnly = TRUE)
pre_vs_post <- tolower(args[1]) == "true"

# function to do the LME model comparison, i.e. to test whether there is a main
# effect of layer. If yes, this implies that the laminar profiles 'go in the same
# direction'
VPF_LMM_check_main_effect <- function(df_long,subfield_levels,pre_vs_post){
  
  REML <- TRUE
  
  results_LMM <- vector("list",length=N_subfields)
  names(results_LMM) <- subfields_levels
  p_values_KR <- list()
  for (subfield in subfields_levels){
    if (subfield=="CA1"){
      #browser()
    }
    inpdata <- df_long[df_long$subfield==subfield,]
    
    #for CA1 and only for memory vs math, build a quadratic model as th shape
    #of the profiles suggests that this would fit them better
    if (subfield=="CA1" & pre_vs_post==FALSE){
      m1 <- lmer(z~ layer + I(layer^2) +(1|subject), data = inpdata, REML = REML,
                 control = lmerControl(optimizer = "Nelder_Mead"))      
    }

    #model with layers
    m1_nosquarelayer <- lmer(z~ layer +(1|subject), data = inpdata, REML = REML,
                             control = lmerControl(optimizer = "Nelder_Mead"))
    #simpler model, just with intercept
    m1_nolayer <- lmer(z~ (1|subject), data = inpdata, REML = REML,
                       control = lmerControl(optimizer = "Nelder_Mead"))
    
    #compare linear model with simple model
    m_comp_KR<-KRmodcomp(m1_nosquarelayer, m1_nolayer)
    
    if (m_comp_KR$stats$p.value<0.05) {
      s <- summary(m1_nosquarelayer)
    } else {
      s <- summary(m1_nolayer)
    }
    
    #in case of CA1, also compare linear model with quadratic model
    if (subfield=="CA1" & pre_vs_post==FALSE & m_comp_KR$stats$p.value<0.05){
      #browser()
      m_comp_KR<-KRmodcomp(m1, m1_nosquarelayer)
      
      if (m_comp_KR$stats$p.value<0.05) {
        s <- summary(m1)
      } else {
        s <- summary(m1_nosquarelayer)
      }
    }
    
    p_values_KR[[subfield]] <- m_comp_KR$test$"p.value"[2]
    
    results_LMM[[subfield]] <- list(p_value = m_comp_KR$test$"p.value"[2], 
                                    F_value = m_comp_KR$test$"stat"[2],
                                    coeff = s$coefficients[,1])
  }
  #correct for multiple comparisons
  p_values_KR = p.adjust(p_values_KR, method="BH")
  
  #feed back corrected p-values
  for (subfield in subfields_levels) {
    results_LMM[[subfield]][[1]][1] <- p_values_KR[[subfield]]
  }
  
  return(results_LMM)
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

if (pre_vs_post==TRUE){
  outname <- "/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/res_LMM_check_main_effect_pre_vs_post.mat"
} else {
  outname <- "/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/res_LMM_check_main_effect_memory_vs_math.mat"
}
writeMat(outname,Subiculum=results$Subiculum,CA1=results$CA1,CA2=results$CA2,
         CA3=results$CA3,CA4DG=results$`CA4/DG`)
