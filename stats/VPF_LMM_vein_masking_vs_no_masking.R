# init ------
pacman::p_load(pacman,R.matlab,car,reshape,dplyr,GGally,ggplot2,ggthemes,ggvis,
               httr,lubridate,plotly,rio,rmarkdown,shiny,stringr,tidyr,lme4,nlme,
               flexplot,foreach,doParallel,fitdistrplus,DescTools,buildmer,pbkrtest,lmerTest,
               rstatix,jtools) 

pre_vs_post <- TRUE



N_cores <- 5
cl <- makeCluster(N_cores,type = "FORK")
registerDoParallel(cores=N_cores)

masked_vs_no_masked_test <- function(inpdata){
  res <- pairwise_wilcox_test(inpdata,z~condition,comparisons=list(c("masked","not masked")),
                              paired=TRUE,exact=TRUE,p.adjust.method="none")
  
  eff_size <- wilcox_effsize(inpdata,z~condition,comparisons=list(c("masked","not masked")),
                             paired=TRUE,exact=TRUE,p.adjust.method="none")
  z_out <- qnorm(res$p/2)
  
  outlist <- list(p_value=res$p,effect_size=as.numeric(eff_size$effsize),z=z_out)
  return(outlist)
}
# read matrix and create data frame  ------------------------------------------
if (pre_vs_post==TRUE) {
  s_masked <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/pre_vs_post_aggregated_masked.mat")
  s_no_masked <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/pre_vs_post_aggregated_no_masked.mat")
  
} else {
  s_masked <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/memory_vs_math_aggregated_masked.mat")
  s_no_masked <- readMat("/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/memory_vs_math_aggregated_no_masked.mat")
  
}

subfields_levels <- c('Subiculum','CA1','CA2','CA3','CA4/DG')
depths <- c('SRLM','inner','midthickness','outer')
conditions <- c("masked","not masked")
N_subfields <- length(subfields_levels)
N_subjects <- dim(s_masked[[1]])[1]/N_subfields

df_masked    <- data.frame(matrix(s_masked[[1]], nrow = N_subjects*N_subfields, ncol=length(depths)))
df_no_masked <- data.frame(matrix(s_no_masked[[1]], nrow = N_subjects*N_subfields, ncol=length(depths)))
colnames(df_masked) <- depths
colnames(df_no_masked) <- depths

subjects <-c(1:N_subjects)
subjects <- rep(subjects,each=N_subfields)
subfields <- rep(subfields_levels,N_subjects)
condition<- rep(conditions,each=N_subjects*N_subfields)

df_masked$subfield <- subfields
df_masked$subject <- subjects

df_no_masked$subfield <- subfields
df_no_masked$subject <- subjects

df <- rbind(df_masked,df_no_masked)
df$condition <- condition

# transform to long ------
df_long <- melt(df,id.vars=c("subfield","subject","condition")) 
names(df_long)[names(df_long)=="variable"] <- "layer"
names(df_long)[names(df_long)=="value"] <- "z"

# create factors --------
df_long$subfield <- factor(df_long$subfield,levels=subfields_levels)
df_long$subject <- factor(df_long$subject,levels=subjects[seq(1, length(subjects), by = N_subfields)])
df_long$layer <- factor(df_long$layer,levels=depths,ordered = TRUE)
df_long$condition <-factor(df_long$condition,level=conditions)
#compare masked vs not masked only at layers where the hypothesis is that three i
#an effect based on the breathhold data------
#subiculum,DG/CA4 = inner
#CA1= SRLM,inner
#CA2,CA3 = outer
results_pairs <- vector("list",length=N_subfields)
names(results_pairs) <- subfields_levels
all_p_values <- vector("list",length=6)
count <- 1
for (subfield in unique(subfields)){
  inpdata_subfield = df_long[df_long$subfield==subfield,]
  
  if (subfield=="CA2" || subfield=="CA3"){
    layer_used <- "outer"
  } else {
    layer_used <- "inner"
  }
  
  inpdata = inpdata_subfield[inpdata_subfield$layer==layer_used,]
  
  outlist <- masked_vs_no_masked_test(inpdata)
  all_p_values[count] <- outlist$p_value
  
  cname <- paste("masked vs. not masked, layer: ",layer_used)
  results_pairs[[subfield]][[cname]] <- outlist
  count <- count+1
}

#SRLM of CA1 is missing. Look at it here
outlist <- masked_vs_no_masked_test(df_long[df_long$subfield=="CA1" & df_long$layer=="SRLM",])
cname <- "masked vs. not masked, layer: SRLM"
all_p_values[length(all_p_values)] <- outlist$p_value
results_pairs[["CA1"]][[cname]] <- outlist

all_p_values <- p.adjust(all_p_values, method="BH")

count <- 1
for (subfield in unique(subfields)){
  results_pairs[[subfield]][[1]]$p_value <- all_p_values[count]
  count <- count + 1
}
results_pairs[['CA1']][[2]]$p_value <- all_p_values[length(all_p_values)]

if (pre_vs_post==TRUE) {
  results_pre_vs_post <- results_pairs
} else {
  results_memory_vs_math <- results_pairs
}
