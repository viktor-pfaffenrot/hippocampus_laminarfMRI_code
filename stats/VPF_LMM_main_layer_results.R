# init ------
pacman::p_load(pacman,R.matlab,car,reshape,dplyr,GGally,ggplot2,ggthemes,ggvis,
               httr,lubridate,plotly,rio,rmarkdown,shiny,stringr,tidyr,lme4,nlme,
               flexplot,foreach,doParallel,fitdistrplus,DescTools,buildmer,pbkrtest,lmerTest,
               rstatix)
N_cores <- 5
cl <- makeCluster(N_cores,type = "FORK")
registerDoParallel(cores=N_cores)

pre_vs_post <- FALSE
# read matrix and create data frame  ------------------------------------------
if (pre_vs_post==TRUE) {
  s <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/pre_vs_post_aggregated_masked.mat")
} else {
  s <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/memory_vs_math_aggregated_masked.mat")
}

subfields_levels <- c('Subiculum','CA1','CA2','CA3','CA4/DG')
depths <- c('SRLM','inner','midthickness','outer')
N_subfields <- length(subfields_levels)
N_subjects <- dim(s[[1]])[1]/N_subfields


df <- data.frame(matrix(s[[1]], nrow = N_subjects*N_subfields, ncol=length(depths)))
colnames(df) <- depths

subjects <-c(1:N_subjects)
subjects <- rep(subjects,each=N_subfields)
subfields <- rep(subfields_levels,N_subjects)

df$subfield <- subfields
df$subject <- subjects

# transform to long ------
df_long <- melt(df,id.vars=c("subfield","subject")) 
names(df_long)[names(df_long)=="variable"] <- "layer"
names(df_long)[names(df_long)=="value"] <- "z"

# create factors --------
df_long$subfield <- factor(df_long$subfield,levels=subfields_levels)
df_long$subject <- factor(df_long$subject,levels=subjects[seq(1, length(subjects), by = N_subfields)])
df_long$layer <- factor(df_long$layer,levels=depths,ordered = TRUE)

#----
# run omnibus model per subfield: Is there a significant main effect of layer?
# I do this by building models with and without layers and then compare them
# with the Kenward-Roger approach for small samples
#https://www.jstatsoft.org/article/view/v059i09


formula_full <-as.formula("z ~ layer + (1|subject)")
formula_nolayer <-as.formula("z ~ (1|subject)")
REML <- TRUE
results_omnibus <- vector("list",length=N_subfields)
names(results_omnibus) <- subfields_levels
p_values_KR <- list()
p_values_PB <- list()

for (subfield in subfields_levels){
  inpdata <- df_long[df_long$subfield==subfield,]
  
  #build full model and model without layer as factor
  m1 <- lmer(formula_full, data = inpdata, REML = REML,
             control = lmerControl(optimizer = "Nelder_Mead"))
  m1_nolayer <- lmer(formula_nolayer, data = inpdata, REML = REML,
                     control = lmerControl(optimizer = "Nelder_Mead"))
  
  #Kenward-Roger approach to compare the models
  m_comp_KR<-KRmodcomp(m1, m1_nolayer)
  
  #parametric bootstrapping as alternative. Works only with complete set of layers, so drop nan:
  inpdata <- na.omit(df_long[df_long$subfield==subfield,])
  #build full model and model without layer as factor
  m1 <- lmer(formula_full, data = inpdata, REML = REML,
            control = lmerControl(optimizer = "Nelder_Mead"))
  m1_nolayer <- lmer(formula_nolayer, data = inpdata, REML = REML,
                    control = lmerControl(optimizer = "Nelder_Mead"))
  
  m_comp_PB<-PBmodcomp(m1, m1_nolayer,nsim=50)
  
  p_values_KR[[subfield]] <- m_comp_KR$test$"p.value"[2]
  p_values_PB[[subfield]] <- m_comp_PB$test$"p.value"[2]
  
  results_omnibus[[subfield]] <- list(
    KR = list(p_value = m_comp_KR$test$"p.value"[2], statistic = m_comp_KR$test$"stat"[2]),
    PB = list(p_value = m_comp_PB$test$"p.value"[2], statistic = m_comp_PB$test$"stat"[2]))
}

#correct for multiple comparisons
p_values_KR = p.adjust(p_values_KR, method="BH")
p_values_PB = p.adjust(p_values_PB, method="BH")


#save corrected p-values
for (subfield in subfields_levels) {
  results_omnibus[[subfield]][[1]][1] <- p_values_KR[[subfield]]
  results_omnibus[[subfield]][[2]][1] <- p_values_PB[[subfield]]
}

#----
#for those subfields with significant main effect, compare layers using wilkoxon signed
#rank test

results_pairs <- vector("list",length=N_subfields)
names(results_pairs) <- subfields_levels
for (subfield in subfields_levels){
  if (results_omnibus[[subfield]][[1]][[1]] > 0.05) {
    results_pairs[[subfield]] <- "ns"
  }
  else {
    inpdata <- na.omit(df_long[df_long$subfield==subfield,])
    N_layer_levels <- length(unique(inpdata$layer))
    
    level_pairs <- combn(N_layer_levels,2)
    
    for (ii in 1:ncol(level_pairs)){
      level1 <- level_pairs[1,ii]
      level1 <- unique(inpdata$layer)[level1]
      level2 <- level_pairs[2,ii]
      level2 <- unique(inpdata$layer)[level2]
      
      x <- inpdata[inpdata$layer==level1, ]$z
      y <- inpdata[inpdata$layer==level2, ]$z
      
     res <- pairwise_wilcox_test(inpdata,z~layer,comparisons=list(c(level1,level2)),
                                 paired=TRUE,exact=TRUE,p.adjust.method="none")
     eff_size <- wilcox_effsize(inpdata,z~layer,comparisons=list(c(level1,level2)),
                                paired=TRUE,exact=TRUE,p.adjust.method="none")
     z_out <- qnorm(res$p/2)
     
     cname <- paste(level1,level2, sep = " vs. ")
     results_pairs[[subfield]][[cname]] <- list(p_value=res$p, 
                                                effect_size=as.numeric(eff_size$effsize),
                                                z=z_out)
    }
    
    #correct for multiple comparisons
    all_p_values <- sapply(results_pairs[[subfield]], function(x) x$p_value)
    all_p_values <- p.adjust(all_p_values, method="BH")
    
    results_pairs[[subfield]] <- Map(function(x, p_value) { x$p_value <- p_value; x }, 
                                     results_pairs[[subfield]], 
                                     all_p_values)
  }
}

