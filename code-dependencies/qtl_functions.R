### Functions for QTL analysis
library(lme4) # for modeling random effects
library(MESS) # for AUC calculation 
library(plyr)
library(shades)
library(scales)
library(extRemes)
library(RColorBrewer)
library(tidyverse)
library(snow)
library(pzfx)
source("code-dependencies/lmmultiresponse.R")

#-------------------------------Miscellaneous----------------------------------#

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# which_chr: given an R/qtl cross object and a marker
#            ID, determines which chromosome the marker
#            is on. 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
which_chr <- function(cross, marker) {
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o] 
  return(chr)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# has_xchr_perms: indicates whether a given permutation 
#                 object has X-chromosome-specific 
#                 permutations  
#+++++++++++++++++++++++++++++++++++++++++++++++++++
has_xchr_perms <- function(perm){
  return('xchr' %in% names(attributes(perm)))
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# add_attributes: adds an attribute to an object  
#+++++++++++++++++++++++++++++++++++++++++++++++++++
add_attribute <- function(obj, attr_name, attr_value){
  attr(obj, attr_name) <- attr_value
  return(obj)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# calc.num.markers: calculates number of markers in 
#                   a cross object
#+++++++++++++++++++++++++++++++++++++++++++++++++++
calc.num.markers <- function(cross){
  p = 0
  for (c in 1:length(cross$geno)){
    p = p + ncol(cross$geno[[c]]$data)
  }
  return(p)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# logit: performs logit transformation (for 
#        proportions/percentages)
#+++++++++++++++++++++++++++++++++++++++++++++++++++
logit <- function(p){log(p/(1-p))}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# trint: performs truncated-RINT transformation 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
trint <- function(y, theta=0.01){
  p <- theta + (rank(y)-1)*(1-2*theta)/(length(y)-1)
  qnorm(p)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# ensure_directory: creates directory if it doesn't 
#                   exist yet
#+++++++++++++++++++++++++++++++++++++++++++++++++++
ensure_directory <- function(directory){
  if(!dir.exists(directory)){
    dir.create(directory);
  }
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# make_logger: creates a logger function 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
make_logger <- function(filename, sep="\n"){
  if(file.exists(filename)){
    file.remove(filename);
  }
  function(...){
    text <- sprintf(...);
    cat(text, file=filename, sep=sep, append=T);
  }
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# calc.dosage: calculates dosage and adds to cross object 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
calc.dosage <- function(cross){
  if (summary(cross)$type == 'bc'){
    for (c in 1:length(cross$geno)){ # for each chr
      cross$geno[[c]]$dos <- 0*cross$geno[[c]]$prob[,,1] + 1*cross$geno[[c]]$prob[,,2]
    }
  } else if (summary(cross)$type == 'f2'){ ### X-chr needs to be validated? 
    for (c in 1:(length(cross$geno)-1)){ # for each autosome
      cross$geno[[c]]$dos <- 0*cross$geno[[c]]$prob[,,1] + 1*cross$geno[[c]]$prob[,,2] + 2*cross$geno[[c]]$prob[,,3]
    }
    cross$geno[[20]]$dos <- 0*cross$geno[[20]]$prob[,,1] + 1*cross$geno[[20]]$prob[,,2] # X-chr; g1 and g2
  }
  
  return(cross)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to add phenotype data to R/qtl data table 
# Input:
#   pheno.name = name of phenotype column in phenotype spreadsheet 
#   new.col = empty column 
#   col.pos = position of new column in R/qtl table 
# Output: dataframe updated with additional phenotype column 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
add_pheno <- function(pheno.name, new.col, col.pos){
  rqtl <- rqtl %>% add_column(placeholder.name = new.col, .before = col.pos)
  names(rqtl)[names(rqtl) == "placeholder.name"] <- pheno.name
  for (i in 3:nrow(rqtl)){
    rqtl[[pheno.name]][i] <- ifelse(rqtl$mouse_ID[i] %in% pheno$Geno_ID, 
                                    pheno[[pheno.name]][which(pheno$Geno_ID == rqtl$mouse_ID[i])], NA)
  }
  return(rqtl)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# readGRMBin: R script to read the GRM binary file
# Input:
#   prefix = file prefix 
# Output: GRM 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
readGRMBin <- function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N, grm=grm))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# fill_in_res:  when a phenotype is missing values, residuals(fit) will be 
#               shorter than the original phenotype vector/column. This 
#               function expans the residuals vector to the original length, 
#               filling in the missing spaces with NAs. 
# Input:
#     pheno: vector of phenotype values (original length)
#     fit: fit with residuals to use in new data vector 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
fill_in_res <- function(pheno, fit){
  has.data <- as.numeric(names(residuals(fit)))
  missing.data <- !seq(1,length(pheno)) %in% has.data # boolean vector (TRUE if missing)  
  
  tmp <- vector(length = length(pheno))
  tmp[missing.data] <- NA
  tmp[has.data] <- residuals(fit)
  
  return(tmp)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# load_themes: loads ggplot themes 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
load_themes <- function(){
  # To find default colors
  #show_col(hue_pal()(2))
  
  poster_theme <<- theme(axis.title.x = element_text(size = 20), 
                   axis.text.x = element_text(size = 18), 
                   axis.title.y = element_text(size = 20),
                   axis.text.y = element_text(size = 18),
                   legend.title = element_text(size = 20),
                   legend.text = element_text(size = 18),
                   plot.title = element_text(size = 22))
  
  poster_theme2 <<- theme(axis.title.x = element_text(size = 22), 
                    axis.text.x = element_text(size = 20), 
                    axis.title.y = element_text(size = 22),
                    axis.text.y = element_text(size = 20),
                    legend.title = element_text(size = 22),
                    legend.text = element_text(size = 20),
                    plot.title = element_text(size = 24))
  
  # For things needing big text  
  big_theme <<- theme(axis.title = element_text(size = 24), 
                      axis.text = element_text(size = 22),
                      plot.title = element_text(size = 24),
                      legend.title = element_text(size = 22),
                      legend.text = element_text(size = 18))
  
  # Same but white background
  bw_biggish_theme <<- theme(axis.title = element_text(size = 20), 
                             axis.text = element_text(size = 18),
                             plot.title = element_text(size = 22),
                             legend.title = element_text(size = 20),
                             legend.text = element_text(size = 18),
                             legend.key = element_rect(fill = "white"),
                             panel.background = element_rect(fill = "white",
                                                             colour = "white"),
                             panel.border = element_rect(colour = "lightgray", fill = NA, size = 1),
                             panel.grid.major = element_line(colour = "lightgray", size = 0.25),
                             panel.grid.minor = element_line(colour = "lightgray", size = 0.1))
  
  # Slightly bigger text
  bw_big_theme <<- theme(axis.title = element_text(size = 22), 
                         axis.text = element_text(size = 20),
                         plot.title = element_text(size = 24),
                         legend.title = element_text(size = 20),
                         legend.text = element_text(size = 18),
                         legend.key = element_rect(fill = "white"),
                         panel.background = element_rect(fill = "white",
                                                         colour = "white"),
                         panel.border = element_rect(colour = "lightgray", fill = NA, size = 1),
                         panel.grid.major = element_line(colour = "lightgray", size = 0.25),
                         panel.grid.minor = element_line(colour = "lightgray", size = 0.1))
  
  # For things needing big text but smaller x-axis labels   
  # (like HS and titer plots with x-axis = strain)
  big_theme2 <<- theme(axis.title = element_text(size = 24), 
                      axis.text.y = element_text(size = 22),
                      axis.text.x = element_text(size = 18),
                      plot.title = element_text(size = 24),
                      legend.title = element_text(size = 22),
                      legend.text = element_text(size = 20))
  
  # same but white bg
  bw_big_theme2 <<- theme(axis.title = element_text(size = 24), 
                       axis.text.y = element_text(size = 22),
                       axis.text.x = element_text(size = 18),
                       plot.title = element_text(size = 24),
                       legend.title = element_text(size = 22),
                       legend.text = element_text(size = 20),
                       legend.key = element_rect(fill = "white"),
                       panel.background = element_rect(fill = "white",
                                                       colour = "white"),
                       panel.border = element_rect(colour = "lightgray", fill = NA, size = 1),
                       panel.grid.major = element_line(colour = "lightgray", size = 0.25),
                       panel.grid.minor = element_line(colour = "lightgray", size = 0.1))
  
  # side-by-side PxG theme 
  sbs_pxg_theme <<- theme(plot.title = element_text(size = 26),
                          axis.text.y = element_text(size = 18),
                          axis.text.x = element_text(size = 17),
                          axis.title = element_text(size = 24))
  
  # For PxG plots needing big text  
  big_pxg_theme <<- theme(axis.title = element_text(size = 28), #26
                      axis.text = element_text(size = 26),
                      legend.title = element_text(size = 26),
                      legend.text = element_text(size = 24))
  
  # themes for Rmarkdown
  rmd_theme <<- theme(axis.title.x = element_text(size = 16), 
                    axis.text.x = element_text(size = 14), 
                    axis.title.y = element_text(size = 16),
                    axis.text.y = element_text(size = 14),
                    legend.title = element_text(size = 16),
                    legend.text = element_text(size = 14),
                    plot.title = element_text(size = 18))
  
  # theme for publication (multiple plots per figure)
  pub_theme <<- theme(legend.key.size = unit(0.03, 'npc'), # legend key is 0.03 of the plot size (basically)
                      axis.title = element_text(size = 8), 
                      axis.text = element_text(size = 6), 
                      legend.title = element_text(size = 8),
                      legend.text = element_text(size = 6),
                      legend.box.spacing = unit(0, 'pt')) # no space between plot and legend 
  
  # theme for publication (one plot)
  pub_theme2 <<- theme(axis.title = element_text(size = 22), 
                     axis.text = element_text(size = 20), 
                     plot.title = element_text(size = 24))
  
  # theme for publication PxG plots 
  pub_pxg_theme <<- theme(axis.title = element_text(size = 22), 
                          axis.text = element_text(size = 18))
  
}


#-----------------------------Phenotype plotting-------------------------------#

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# trajectory_plot: plots trajectory data from an R/QTL cross object
# Input:
#     cross.obj: R/qtl cross object 
#     steps: x-axis data 
#     phenos: y-axis data. Can either be the name of a column in cross.obj$pheno, 
#             or an integer which will be repeated as the phenotype. 
# Output: trajectory plot 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
cov_trajectory_plot <- function(cross, phenos, title = NULL, ylab = NULL, 
                                incl.parents = TRUE, parent.lty = 1, ylim = NULL,...){
  # load parent data 
  if (incl.parents == TRUE){
    CC_weight <- read_pzfx("source_data/screen/CC_mice.pzfx", "CC_weight")
    CC006_weight <- CC_weight[,c(10:13)] #16750
    CC044_weight <- CC_weight[,c(22:25)] #4410
    CC006dat <- data.frame(mouse_ID = 'CC006', 
                           dpi = c('0', '1', '2', '3', '4'), 
                           weight = rowMeans(CC006_weight, na.rm=T))
    CC044dat <- data.frame(mouse_ID = 'CC044', 
                           dpi = c('0', '1', '2', '3', '4'), 
                           weight = rowMeans(CC044_weight, na.rm=T))
  }
  
  # data setup for ggplot 
  dat <- cross$pheno$mouse_ID
  for (i in 1:length(phenos)){
    y <- cross$pheno[phenos[i]]
    dat <- cbind(dat, y)
  }
  dat <- as.data.frame(dat)
  colnames(dat) <- c('mouse_ID', '0', '1', '2', '3', '4')
  dat <- pivot_longer(dat, cols = c(2:6), names_to = "dpi", values_to = "weight")
  
  if (incl.parents == TRUE){
    wtloss <- ggplot(dat, aes(x = dpi, y = weight)) +
      geom_line(aes(group = mouse_ID, color = "F2", linetype = "F2", size = "F2")) +
      labs(x = "Days post-infection", y = ylab, title = title, color = "Legend", linetype = "Legend", size = "Legend") +
      geom_line(data = CC006dat, aes(x = dpi, y = weight, group = mouse_ID, color = "CC006", linetype = "CC006", size = "CC006")) + 
      geom_line(data = CC044dat, aes(x = dpi, y = weight, group = mouse_ID, color = "CC044", linetype = "CC044", size = "CC044")) +
      scale_linetype_manual(values = c("F2" = 1, "CC006" = parent.lty, "CC044" = parent.lty)) +
      scale_size_manual(values = c("F2" = 0.2, "CC006" = 1.3, "CC044" = 1.3)) + 
      scale_color_manual(values = c("F2" = "gray63", "CC006" = "#E41A1C", "CC044" = "#377EB8")) +
      {if(!is.null(ylim)) ylim(ylim) }+
      scale_x_discrete(expand = c(0.03,0.03))
  } else { # don't include parent trajectories 
    wtloss <- ggplot(dat, aes(x = dpi, y = weight)) +
      geom_line(aes(group = mouse_ID), color = "gray63", size = 0.2) +
      {if(!is.null(ylim)) ylim(ylim) }+
      labs(x = "Days post-infection", y = ylab, title = title) +
      scale_x_discrete(expand = c(0.03,0.03))
  }

  return(wtloss)
}


#--------------------------Calculate derived measures--------------------------#
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# calc_auc: this function calculates an "area under the curve" metric, which is 
# really the area above the curve and below the horizontal line at the first data 
# point. 
# Input:
#     cross: r/qtl cross object 
#     col.name: column name to be used for the new aggregate phenotype 
#     steps: x-axis data 
#     phenos: y-axis data - phenotypes to be used to calculate the new aggregate 
#             phenotype
# Output:
#     cross: r/qtl cross object updated with new aggregate phenotype 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
calc_auc <- function(cross, steps, col.name, phenos){
  cross$pheno[,col.name] <- rep(NA,nrow(cross$pheno))
  
  #avgs <- colMeans(cross$pheno[,phenos], na.rm = T)
  
  for (i in 1:nrow(cross$pheno)){
    line <- cross$pheno[i,phenos]
    auc <- auc(x = steps, y = line)
    aac <- line[1]*tail(steps, n=1) - auc
    cross$pheno[i,col.name] <- aac 
  }
  
  return(cross)
}


#----------------------------Working with QTL objects--------------------------#

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# create_models: this function creates all single-QTL models specified in the 
#                models table
# Input:
#     cross.obj: r/qtl cross object 
#     models: table containing information about phenotypes to create single-QTL 
#             models for 
#     covar: dataframe containing covariates 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
create_models <- function(cross.obj, models, covar = NULL, mod.dir){
  ensure_directory(mod.dir)
  
  for (i in 1:nrow(models)){
    filepath = paste(c(mod.dir, models[i,'obj'], ".Rdata"), collapse = "")
    
    if (file.exists(filepath)){ # model exists already, load .Rdata object 
      print(paste("Loading scanone object from", filepath))
      load(filepath, envir = .GlobalEnv)
    } else { # file doesn't exist, create model 
      print(paste("Creating scanone object for", models[i,'name'], "using", models[i,'type'], "model."))
      
      pheno.col <- models$colname[i]
      mod.type <- models$type[i]
      
      # define covariates, if specified
      covlist <- str_split(models$cov[i], ",")[[1]]
      addcovar <- NULL
      if(!is.na(models$cov[i])){
        addcovar <- covar[,covlist]
      }
      
      # create model 
      assign(models[i,'obj'], scanone(cross.obj, pheno.col = pheno.col, model = mod.type, 
                                      method = "hk", addcovar = addcovar), 
             envir = .GlobalEnv) 
      
      print(paste("Saving object to", filepath))
      save(list = models[i,'obj'], file = filepath)
    }
  }
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# create_perms: this function creates permutation objects for all single-QTL models 
#               specified in the models table 
# Input:
#     cross.obj: r/qtl cross object 
#     models: table containing information about phenotypes with single-QTL models 
#     covar: covariates 
#     perm.dir: directory with existing perm objects / to save perm objects to 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
create_perms <- function(cross.obj, models, perm.dir, n.perm = 1000,
                         perm.strata = NULL, covar = NULL, perm.Xsp = FALSE, ...){

  ensure_directory(perm.dir)
  
  #--------------------Create or load permutation object-----------------------#
  for (i in 1:nrow(models)){
    perm <- models$perm.obj[i]
    perm.fname <- gsub("\\.", "", perm) # remove period from permutation object name, if present
    permfile <- paste(perm.dir, perm.fname, ".Rdata", sep = "")
    
    # create or load permutation object
    if(file.exists(permfile)){ # permutation object already exists, load from Rdata file 
      load(permfile, envir = .GlobalEnv)
      print(paste("Loaded", perm, "from file."))
    } else {
      print(paste("Creating", perm, "..."))
      
      covlist <- str_split(models$cov[i], ",")[[1]]
      addcovar <- NULL
      if(!is.na(models$cov[i])){
        addcovar <- covar[,covlist]
      }
      
      set.seed(124)
      assign(perm, 
             scanone(cross.obj, pheno.col = models[i, 'colname'], method = "hk",
                     model = models[i,'type'], addcovar = addcovar, n.perm = n.perm,
                     perm.strata = perm.strata, perm.Xsp = perm.Xsp, n.cluster = 4), 
             envir = .GlobalEnv)

      save(list = perm, file = permfile)
    }
  }
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# doc_peaks: this function finds significant LOD peaks, calculates 95% Bayes credible 
#           intervals for those peaks, and creates a summary table with data for 
#           each one, including marker, chromosome, position, LOD and positions of 
#           Bayes credible intervals.
# Input: 
#     models: table containing information about phenotypes with single-QTL models 
#     sig.level: significance level (can be 0.05 or 0.10)
#+++++++++++++++++++++++++++++++++++++++++++++++++++
doc_peaks <- function(models, sig.level = 0.10){
  
  #--------------------Make sure all permutations exist------------------------#
  perms.exist = T
  for (m in 1:nrow(models)){
    if (!exists(models$perm.obj[m])){
      print(paste(models$perm.obj[m], "does not exist. All permutation objects must exist before attempting to identify significant peaks."))
      perms.exist = F
    }
  }
  if (perms.exist == F){stop("Missing permutation objects, execution halted.")}
  
  #-------------------------Count significant peaks----------------------------# 
  num.peaks = 0
  num.peaks.per.mod <- rep(NA, nrow(models))
  gevnames <- vector(mode = 'list', length = nrow(models))
  thresholds <- vector(mode = 'list', length = nrow(models))
  sig.ix <- ifelse(sig.level == 0.10, 2, ifelse(sig.level == 0.05, 1, NA)) ### Handle more than just 0.05 and 0.10
  for (m in 1:nrow(models)){
    mod <- models$obj[m]
    perm <- models$perm.obj[m]
    
    if ((has_xchr_perms(get(perm))) & (ncol(get(perm)$A) == 1)) { # X-chr specific thresholds, 1 LOD 
      num.peaks.per.mod[m] <- nrow(summary(get(mod), alpha = sig.level, perms = get(perm)))
    } else if ((has_xchr_perms(get(perm))) & (ncol(get(perm)$A) > 1)) { # X-chr specific thresholds >1 LOD (2part)
      num.peaks.per.mod[m] <- nrow(summary(get(mod), alpha = sig.level, perms = get(perm)))
    }
    
  }
  num.peaks <- sum(num.peaks.per.mod)
  
  
  #---------------------Get data for significant peaks-------------------------# 
  peak.data <- c('model', 'marker', 'chr', 'pos', 'lod', 'Bayes CI')
  peaks <- matrix(ncol = length(peak.data), nrow = num.peaks)
  colnames(peaks) <- peak.data
  
  row.ix = 0 # row in peaks table
  for (m in 1:nrow(models)){ # for each phenotype / model
    mod <- models$obj[m]
    perm <- models$perm.obj[m]
    
    if (num.peaks.per.mod[m] == 0){
      print(paste("For model", models[m,'obj'], ", no LOD peaks above threshold"))
    } else {
      for (p in 1:num.peaks.per.mod[m]){ # for each peak 
        row.ix = row.ix + 1 

        chrom = summary(get(mod), alpha = 0.10, perms = get(perm))[p,1]

        # bayesint() produces a table with 3 rows (lower interval, peak, upper interval), columns = chr, pos, lod 
        bayesCI <- bayesint(get(models[m,'obj']), chr = chrom, prob = 0.95)
        # Populate peaks dataframe 
        peaks[row.ix,] <- c(mod, rownames(bayesCI)[2], chrom, bayesCI[2,2], bayesCI[2,3], 
                            paste(bayesCI[1,2],"-",bayesCI[3,2]))
      } 
    }
  }
  
  peaks <- as.data.frame(peaks)
  
  # Remove leading/trailing white space 
  for (i in 1:ncol(peaks)){
    peaks[,i] <- trimws(peaks[,i], which = c("both"))
  }
  
  return(peaks)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# get_unadj_pval: this function gets the unadjusted p-value from the linear regression
#           between a particular phenotype and the genotype at a particular marker.
# Input: 
#     cross: R/qtl cross object 
#     pheno.col: phenotype column name 
#     marker: QTL marker 
#     covar_df: dataframe of covariates 
#     covar_names: list of covariates to use in fitting the model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
get_unadj_pval <- function(cross, pheno.col, marker, covar_df = NULL, covar_names = NULL){
  chr <- find.markerpos(cross, marker)[1,1]
  
  cross.imp <- fill.geno(cross) 
  marker.geno <- pull.geno(cross.imp, chr = chr)[,marker]
  
  if (!is.null(covar_names)){
    covlist <- str_split(covar_names, ",")[[1]]
  }
  
  formula <- paste0(pheno.col, " ~ ", 'marker.geno', 
                    ifelse(is.null(covar_names), "", " + "),
                    paste(covlist, collapse = " + "))
  
  if (is.null(covar_names)){
    lmdata <- data.frame(cross$pheno[,pheno.col], marker.geno)
    names(lmdata)[1] <- pheno.col # have to rename column 
  } else {
    lmdata <- data.frame(cross$pheno[,pheno.col], 
                         covar_df[,covlist], 
                         marker.geno)
    names(lmdata)[1] <- pheno.col
  }
  
  lmod <- lm(formula, lmdata)
  unadj_pval <- summary(lmod)$coefficients['marker.geno', 'Pr(>|t|)']
  return(unadj_pval)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# get_var_expl: this function gets the variance explained for a QTL.
# Input: 
#     cross: R/qtl cross object 
#     pheno.col: phenotype column name 
#     marker: QTL marker 
#     covar_df: dataframe of covariates 
#     covar_names: list of covariates to use in fitting the model 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
get_var_expl <- function(cross, pheno.col, marker, covar_df = NULL, covar_names = NULL){
  pheno <- pull.pheno(cross, pheno.col)
  chr <- find.markerpos(cross, marker)[1,1]
  
  cross.imp <- fill.geno(cross) 
  marker.geno <- pull.geno(cross.imp, chr = chr)[,marker]
  
  if (!is.null(covar_names)){
    covlist <- str_split(covar_names, ",")[[1]]
  }
  
  formula <- paste0(pheno.col, " ~ ", 'marker.geno', 
                    ifelse(is.null(covar_names), "", " + "),
                    paste(covlist, collapse = " + "))
  
  if (is.null(covar_names)){
    lmdata <- data.frame(cross$pheno[,pheno.col], 
                         marker.geno)
    names(lmdata)[1] <- pheno.col  
  } else {
    lmdata <- data.frame(cross$pheno[,pheno.col], 
                         covar_df[,covlist], 
                         marker.geno)
    names(lmdata)[1] <- pheno.col
  }
  lmdata$marker.geno <- as.factor(lmdata$marker.geno)
  
  lmod <- lm(formula, lmdata)
  
  sigma2 <- summary(lmod)$sigma^2

  if (summary(cross)$type == 'bc'){
    uAA <- mean(pheno[which(marker.geno == 1)], na.rm = TRUE) # 1 = AA
    uAB <- mean(pheno[which(marker.geno == 2)], na.rm = TRUE) # 2 = AB
    a <- uAB - uAA
    h2 <- a^2 / (a^2 + 4*sigma2)
  } else if (summary(cross)$type == 'f2'){
    uAA <- mean(pheno[which(marker.geno == 1)], na.rm = TRUE) # 1 = AA
    uAB <- mean(pheno[which(marker.geno == 2)], na.rm = TRUE) # 2 = AB
    uBB <- mean(pheno[which(marker.geno == 3)], na.rm = TRUE) # 3 = BB 
    a <- (uBB - uAA)/2
    d <- uAB - (uAA + uBB)/2
    h2 <- (2*a^2 + d^2)/(2*a^2 + d^2 + 4*sigma2)
  }

  return(h2)
}


#-----------------------------------Plotting-----------------------------------#

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# plot_scans: this function plots genome scans for all models specified in the 
#             models table. Adds a solid line for 5 % significance threshold and 
#             a dotted line for 10 % significance threshold, if the permutation
#             object exists. 
# Input:
#     models: table containing information about phenotypes with single-QTL models 
#     ylim: y limits, if they need to be consistent across plots. If not set, will 
#           use the default. 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
plot_scans <- function(models, save = FALSE, save.dir = NULL, ...){
  
  for (m in 1:nrow(models)){
    mod_name <- models$obj[m]
    mod <- get(mod_name)
    
    if (ncol(mod)>3) {lodcolumn=1:3} else {lodcolumn=1} # if ncol>3, multi-LOD model 
    
    if (save == TRUE){ # saving, so open connection to png and make axis labels bigger 
      fname <- paste0(save.dir, mod_name, '.png')
      png(fname, width = 750)
      plot(mod, ylab = "", xlab = "", main = models$name[m], lodcolumn=lodcolumn,
           alternate.chrid = T, bandcol = "gray90", cex.main = 2, cex.axis = 2, ...) # cex.main was 1.8, cex.axis was 1.5
      title(ylab = "LOD", line = 2.5, cex.lab = 2)
      title(xlab = "Chromosome", cex.lab = 2.3, line = 3) # cex.lab was 2
    } else { # not saving, just printing
      plot(mod, ylab = "LOD", main = models$name[m], lodcolumn=lodcolumn,
           alternate.chrid = T, bandcol = "gray90", ...)
    }
    
    # Add significance thresholds 
    perm <- models$perm.obj[m]

    if (exists(perm)){
      perm.sum = summary(get(perm))
      abline(h = perm.sum$A[,1], lty=1:2) # lod.p.mu
    }
    
    if (save == TRUE){ dev.off() } # close connection 
  }
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to plot phenotype x genotype plots for significant peaks
# Input:
#     cross.obj: r/qtl cross object 
#     cross.type: bc or f2
#     raw.data: raw phenotypes (cross.obj$pheno before transformations)
#     geno.map: list mapping A and B allele to more informative labels
#     peaks: dataframe of significant peaks 
#     plot.type: type of plot to produce for each marker (can be "effect" or "pxg")
#     qtl_map: dataframe mapping chromosome to QTL name (for plot titles)
#     theme: ggplot theme to use for PxG
#+++++++++++++++++++++++++++++++++++++++++++++++++++
plot_pxg <- function(cross.obj, cross.type, raw.data, peaks, plot.type, qtl_map, 
                     theme, type, save = FALSE, save.dir = NULL,...){
  for (k in 1:nrow(peaks)){
    
    if (save == TRUE){
      plotname <- paste0(peaks[k,'model'], '-chr', peaks[k,'chr'])
      fname <- paste0(save.dir, plotname, '.png')
      png(fname, width = 550)
    }
    
    if (plot.type == "effect"){
      effectplot(cross.obj, mname1 = peaks[k,'marker'], 
                 pheno.col = models[which(models$obj == peaks[k,'model']), 'colname'], 
                 ylab = models[which(models$obj == peaks[k,'model']), 'abbr'], 
                 main = paste(peaks[k,'marker'], " (chr", peaks[k,'chr'], ")", sep=""))
    } else if (plot.type == "pxg"){
      pheno.col = models[which(models$obj == peaks[k,'model']), 'colname']
      pheno = raw.data[,pheno.col]
      marker.name = peaks[k,'marker']
        
      p <- pxg(cross.obj, pheno = pheno, 
               marker = marker.name,  
               ylab = models[which(models$obj == peaks[k,'model']), 'abbr'], 
               theme = theme, ...)
      print(p)
    }
    
    if (save == TRUE){ dev.off() }
  }
}





#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to plot phenotype x genotype plots with ggplot. 
# Alternative to qtl::plotPXG, uses ggplot and more intuitive error bars  
# Input:
#     cross: r/qtl cross object 
#     pheno: vector of phenotype values for all mice 
#     marker: marker whose genotype to plot
#     geno.map: list mapping A and B allele to more informative labels 
#     qtl.map: list mapping chromosomes to QTL names for x-axis label
#     xlab: x axis label
#     ylab: y axis label 
#     title: title
#     theme: ggplot theme
#     type: scatter (default), boxplot or violin 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
pxg <- function(cross, pheno, marker, geno.map, qtl.map = NULL, xlab = NULL, ylab = NULL, ylim = NULL,
                title = NULL, bestfit = TRUE, theme = rmd_theme, type = 'scatter', ...){
  # cross type 
  cross.type <- class(cross)[1]
  
  # what chromosome is marker on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o]
  
  # linear regression (best-fit line only appropriate for autosomes) 
  if ((chr == 'X')&('M' %in% levels(cross$pheno$sex))){bestfit = FALSE} 
  
  if(is.null(xlab) & !is.null(qtl.map)){
    qtl.name <- qtl.map$qtl_name[qtl.map$chr == chr]
    xlab <- as.expression(bquote(italic(.(qtl.name))*" (chr"*.(chr)*") "*Genotype))
  } else if (is.null(xlab) & is.null(qtl.map)){
    xlab <- "Genotype"
  }
  
  cross.imp <- fill.geno(cross) 
  marker.genos <- cross.imp$geno[[chr]]$data[,marker]
  
  if (bestfit == TRUE){
    fit <- lm(pheno ~ marker.genos)
    int <- summary(fit)$coefficients[1,1]
    slope <- summary(fit)$coefficients[2,1]
    slope.pval <- summary(fit)$coefficients[2,4]
  }
  
  df <- data.frame(geno = marker.genos,
                   pheno = pheno)
  
  # create genotype labels from A/B mapping 
  if (cross.type == 'bc'){ # backcross, autosome or X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr != 'X')){ # F2 cross, autosome
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"))
  } else if ((cross.type == 'f2') & (chr == 'X')){ # F2 cross, X
    geno.labels = c(paste(geno.map$A, geno.map$A, sep = "/"),
                    paste(geno.map$A, geno.map$B, "f", sep = "/"),
                    paste(geno.map$A, geno.map$B, "r", sep = "/"),
                    paste(geno.map$B, geno.map$B, sep = "/"),
                    paste(geno.map$A, "Y", sep = "/"),
                    paste(geno.map$B, "Y", sep = "/"))
  }
  
  # map genotypes to genotype labels                 
  if (cross.type == 'bc'){ 
    df$geno <- mapvalues(df$geno, from = c(1,2), to = geno.labels)
  } else { # if cross.type = 'f2'
    if (chr != 'X'){
      df$geno <- mapvalues(df$geno, from = c(1,2,3), to = geno.labels) 
    } else { # X-chromosome marker 
      pgm <- pull.pheno(cross, "pgm")
      sex <- as.numeric(pull.pheno(cross, "sex") == "M")
      X.data <- data.frame(sex = sex, geno = df$geno, pgm = pgm)
      X.data$X.geno <- with(df, 
                            ifelse(sex==0 & geno==1 & pgm==0, 'AA', 
                            ifelse(sex==0 & geno==1 & pgm==1, 'BB', 
                            ifelse(sex==0 & geno==2 & pgm==0, 'ABf', 
                            ifelse(sex==0 & geno==2 & pgm==1, 'ABr', 
                            ifelse(sex==1 & geno==1, 'AY', 'BY'))))))
      
      df$geno <- mapvalues(X.data$X.geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                           to = geno.labels)
    }
  }
  
  df$geno <- factor(df$geno, levels = geno.labels)
  
  dfsum <- data.frame(geno = geno.labels,
                      mean = aggregate(df$pheno, by = list(df$geno), FUN = mean, na.rm = TRUE)$x, 
                      sd = aggregate(df$pheno, by = list(df$geno), FUN = sd, na.rm = TRUE)$x)
  
  if (type == 'scatter'){
    p <- ggplot(data = df, mapping = aes(x = geno, y = pheno)) + 
      geom_jitter(width = 0.1, height = 0) + 
      scale_y_continuous(limits = ylim) + 
      geom_errorbar(data = dfsum, mapping = aes(x = geno, y = mean, ymin = mean, ymax = mean), 
                    width = 0.2, size = 1, alpha = 0.5) + 
      geom_errorbar(data = dfsum, mapping = aes(x = geno, y = mean, ymin = mean-sd, ymax = mean+sd), 
                    width = 0.3, size = 1, alpha = 0.5) +
      labs(title = title, y = ylab, x = xlab) + 
      {if (bestfit == TRUE) geom_abline(slope = slope, intercept = int, color = 'red')} + 
      theme
    # bestfit line color was brewer.pal(n=5,name='Dark2')[4]
    # sd bar color was brewer.pal(n=5,name='Dark2')[5]
  } else if (type == 'boxplot'){
    p <- ggplot(data = df, mapping = aes(x = imp_geno, y = pheno)) +
      geom_boxplot(fill = "#1B9E77", notch = TRUE) + # add geom_jitter(width = 0.1, height = 0) to do a dotplot 
      labs(title = title, y = ylab, x = xlab) +
      theme
  } else if (type == 'violin'){
    p <- ggplot(data = df, mapping = aes(x = imp_geno, y = pheno)) +
      geom_violin(fill = "#1B9E77") + # add geom_jitter(width = 0.1, height = 0) to do a dotplot 
      labs(title = title, y = ylab, x = xlab) +
      theme
  }
  
  
  if ((type == 'scatter') & (bestfit == TRUE)){ # only do best fit line on scatter plots 
    print(paste("slope p-value: ", slope.pval))
  }
  
  return(p)
}
