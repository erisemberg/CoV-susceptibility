### Functions for QTL analysis (SARS and peanut allergy)
library(lme4) # for modeling random effects
library(MESS) # for AUC calculation 
library(plyr)
library(shades)
library(scales)
library(extRemes)
library(RColorBrewer)
library(tidyverse)
library(snow)
source("/code-dependencies/lmmultiresponse.R")


#-------------------------------Miscellaneous----------------------------------#
which_chr <- function(cross, marker) {
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o] 
  return(chr)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# has_xchr_perms: indicates whether permutatino object has x-chr-specific permutations  
#+++++++++++++++++++++++++++++++++++++++++++++++++++
has_xchr_perms <- function(perm){
  return('xchr' %in% names(attributes(perm)))
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# add_attributes: add an attribute to an object  
#+++++++++++++++++++++++++++++++++++++++++++++++++++
add_attribute <- function(obj, attr_name, attr_value){
  attr(obj, attr_name) <- attr_value
  return(obj)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# logit transformation (for proportions/percentages)
#+++++++++++++++++++++++++++++++++++++++++++++++++++
logit <- function(p){log(p/(1-p))}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Will's truncated-RINT transformation 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
trint <- function(y, theta=0.01){
  p <- theta + (rank(y)-1)*(1-2*theta)/(length(y)-1)
  qnorm(p)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to ensure directory exists before creating files in it 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
ensure_directory <- function(directory){
  if(!dir.exists(directory)){
    dir.create(directory);
  }
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to create logger function 
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
# Function to add phenotype data to R/qtl data table 
# Input:
#   pheno.name = name of phenotype column in phenotype spreadsheet 
#   new.col = empty column 
#   col.pos = position of new column in R/qtl table 
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
# Output:
#     cross: r/qtl cross object updated with new aggregate phenotype 
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
# fill_in_res:  when a phenotype is missing values (e.g. flow data in SARS 
#               experiments), residuals(fit) will be shorter than the original 
#               phenotype vector/column. This function expans the residuals 
#               vector to the original length, filling in the missing spaces 
#               with NAs. 
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
                             panel.border = element_rect(colour = "lightgray", fill = NA, linewidth = 1),
                             panel.grid.major = element_line(colour = "lightgray", linewidth=0.25),
                             panel.grid.minor = element_line(colour = "lightgray", linewidth=0.1))
  
  # Slightly bigger text
  bw_big_theme <<- theme(axis.title = element_text(size = 22), 
                         axis.text = element_text(size = 20),
                         plot.title = element_text(size = 24),
                         legend.title = element_text(size = 20),
                         legend.text = element_text(size = 18),
                         legend.key = element_rect(fill = "white"),
                         panel.background = element_rect(fill = "white",
                                                         colour = "white"),
                         panel.border = element_rect(colour = "lightgray", fill = NA, linewidth = 1),
                         panel.grid.major = element_line(colour = "lightgray", linewidth=0.25),
                         panel.grid.minor = element_line(colour = "lightgray", linewidth=0.1))
  
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
                       panel.border = element_rect(colour = "lightgray", fill = NA, linewidth = 1),
                       panel.grid.major = element_line(colour = "lightgray", linewidth=0.25),
                       panel.grid.minor = element_line(colour = "lightgray", linewidth=0.1))
  
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
  # add axis.title.y = element_text(hjust = 0) or something to move y-axis away from plot
  
}



#-----------------------------Phenotype plotting-------------------------------#

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# trajectory_plot: plots trajectory data from a r/QTL cross object
# Input:
#     cross.obj: r/qtl cross object 
#     steps: x-axis data 
#     phenos: y-axis data - can either be the name of a column in cross.obj$pheno, 
#             or an integer which will be repeated as the phenotype. 
### to-do: 
###     1. Convert to ggplot  
#+++++++++++++++++++++++++++++++++++++++++++++++++++
cov_trajectory_plot <- function(cross, phenos, title = NULL, ylab = NULL, 
                                incl.parents = TRUE, parent.lty = 1, ylim = NULL,...){
  # load parent data 
  if (incl.parents == TRUE){
    CC_weight <- read_pzfx("source_data/CC_mice.pzfx", "CC_weight")
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
#     method: whether to use existing R/qtl functions ("rqtl", default);
#             linear mixed modeling code ("lmm"); or
#             gene-by-treatment effect modeling code ("gxt")
#+++++++++++++++++++++++++++++++++++++++++++++++++++
create_models <- function(cross.obj, models, covar = NULL, use = "rqtl", 
                          trt = NULL, mod.dir){
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
      
      if ((use == "lmm") & (mod.type == "normal")){ 
        assign(models[i,'obj'], hk(cross.obj, pheno.col = pheno.col, envir = .GlobalEnv))
      } else if ((use == "gxt") & (mod.type == "normal")){  
        assign(models[i,'obj'], gxt(cross.obj, pheno.col = pheno.col, 
                                    covar = covar, trt = trt), envir = .GlobalEnv)
      } else { # use rqtl::scanone
        
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
      }
      
      print(paste("Saving object to", filepath))
      save(list = models[i,'obj'], file = filepath)
    }
  }
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++
# summary.permgev: summary() function for permutation object
#+++++++++++++++++++++++++++++++++++++++++++++++++++
summary.permgev <- function(object, ...){
  # Calculate 5% and 10% significance thresholds 
  thresholds <- matrix(qevd(p = c(0.95, 0.90), 
                            loc = object$results$par[1], 
                            scale = object$results$par[2], 
                            shape = object$results$par[3]))
  colnames(thresholds) <- c('lod')
  rownames(thresholds) <- c('5%', '10%')
  
  return(thresholds)
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
    
    
    #------------------------Fit GEV distribution(s)---------------------------#  
    if (!has_xchr_perms(get(perm))){ # no X-chr-specific thresholds 
      # 1 GEV distribution
      gevname <- paste0(perm, '.gev')
      print(paste("Creating GEV fit", gevname, "..."))
      fitgev <- fevd(as.data.frame(get(perm))$lod, type = "GEV")
      assign(gevname, fitgev, envir = .GlobalEnv)
      # add custom 'permgev' class so that we can use custom summary function to get GEV-defined thresholds  
      assign(gevname,
             add_attribute(get(gevname), 'class', c('permgev', 'fevd')),
             envir = .GlobalEnv)
    } else if ((has_xchr_perms(get(perm))) & (ncol(get(perm)$A) == 1)) { # X-chr thresholds, 1 lod column
      # 2 GEV distributions - 1 for autosomes, 1 for X-chr
      for (chr in c('A', 'X')){
        gevname <- paste0(perm, chr, '.gev')
        print(paste("Creating GEV fit", gevname, "..."))
        fitgev <- fevd(as.data.frame(get(perm)[[chr]])$lod, type = "GEV")
        assign(gevname, fitgev, envir = .GlobalEnv)
        # add custom 'permgev' class for use of custom summary function to get GEV-defined thresholds  
        assign(gevname,
               add_attribute(get(gevname), 'class', c('permgev', 'fevd')),
               envir = .GlobalEnv)
      }
    } else if ((has_xchr_perms(get(perm))) & (ncol(get(perm)$A) > 1)){ # X-chr specific thresholds, >1 LOD column (e.g. 2part model)
      # 6 GEV distributions - 3 for autosomes (1 for each LOD score); 3 for X-chr
      for (chr in c('A', 'X')){
        for (lod.type in colnames(get(perm)[[chr]])){
          gevname <- paste0(perm, chr, '.', lod.type, '.gev')
          print(paste("Creating GEV fit", gevname, "..."))
          fitgev <- fevd(as.data.frame(get(perm)[[chr]])[[lod.type]], type = "GEV")
          assign(gevname, fitgev, envir = .GlobalEnv)
          # add custom 'permgev' class for use of custom summary function to get GEV-defined thresholds  
          assign(gevname,
                 add_attribute(get(gevname), 'class', c('permgev', 'fevd')),
                 envir = .GlobalEnv)
        }
      }
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
    
    if (!has_xchr_perms(get(perm))){ # no X-chr specific thresholds (e.g. backcross)
      
      # not using GEV 
      num.peaks.per.mod[m] <- nrow(summary(get(mod), alpha = sig.level, perms = get(perm)))
      
      ### Using GEV - copy from code on peanut GitHub 
      # gevnames[[m]] <- paste0(perm, '.gev') # GEV dist fit from perm object
      ### this method of defining length 20 thresholds won't work after Karl fixes bug 
      #thresholds[[m]] <- summary(get(gevnames[[m]]))[sig.ix]
      #num.peaks.per.mod[m] <- nrow(summary(get(mod), threshold = thresholds[[m]]))
      
      
    } else if ((has_xchr_perms(get(perm))) & (ncol(get(perm)$A) == 1)) { # X-chr specific thresholds, 1 LOD 
      
      # not using GEV 
      # num.peaks.per.mod[m] <- nrow(summary(get(mod), alpha = sig.level, perms = get(perm)))
      
      # using GEV but not handling X-chr specific thresholds (fine for virus research paper bc no X-chr QTL)
      gevnames[[m]] <- c(paste0(perm, 'A.gev'), paste0(perm, 'X.gev'))
      thresholds[[m]] <- summary(get(gevnames[[m]][1]))[sig.ix]
      num.peaks.per.mod[m] <- nrow(summary(get(mod), threshold = thresholds[[m]]))
      
      ### Handle X-chr specific thresholds 
      # this method of defining thresohlds won't work after Karl fixes bug 
      # thresholds[[m]] <- c(rep(summary(get(gevnames[[m]][1]))[sig.ix], 19), 
      #                    summary(HSpermX.gev)[sig.ix])
      # num.peaks.per.mod[m] <- nrow(summary(get(mod), threshold = thresholds[[m]]))
      
    } else if ((has_xchr_perms(get(perm))) & (ncol(get(perm)$A) > 1)) { # X-chr specific thresholds >1 LOD (2part)
      
      # not using GEV 
      num.peaks.per.mod[m] <- nrow(summary(get(mod), alpha = sig.level, perms = get(perm)))
      
      # gevnames[[m]] <- c(paste0('Titerperm', c(rep('A.',3), rep('X.', 3)), colnames(Titerperm$A), '.gev'))
      ### Need to define thresholds[[m]] in order for later code in doc_peaks to work
      ### This threshold code doesn't work though 
      # thresholds[[m]] <- data.frame(lod.p.mu = c(rep(summary(get(gevnames[[m]][1]))[1], 19), summary(get(gevnames[[m]][4]))[1]),
      #                   lod.p = c(rep(summary(get(gevnames[[m]][2]))[1], 19), summary(get(gevnames[[m]][5]))[1]),
      #                   lod.mu = c(rep(summary(get(gevnames[[m]][3]))[1], 19), summary(get(gevnames[[m]][6]))[1]))
      # num.peaks.per.mod[m] <- nrow(summary(get(mod), threshold = thresholds[[m]]))
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
    #gevs <- gevnames[[m]]
    
    if (num.peaks.per.mod[m] == 0){
      print(paste("For model", models[m,'obj'], ", no LOD peaks above threshold"))
    } else {
      for (p in 1:num.peaks.per.mod[m]){ # for each peak 
        row.ix = row.ix + 1 
        
        # not using GEV 
        #chrom = summary(get(mod), alpha = 0.10, perms = get(perm))[p,1]
        
        # using GEV - one threshold 
        chrom = summary(get(mod), threshold = thresholds[[m]])[p,1]
        
        # using GEV - multiple thresholds 
        #chrom = summary(get(mod), threshold = summary(get(gevs[[1]]))[sig.ix])[p,1]
        
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
#     cross: 
#     pheno.col: 
#     marker:
#     covariates 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
get_unadj_pval <- function(cross, pheno.col, marker, covar = NULL){
  pheno <- pull.pheno(cross, pheno.col) ### COMMENT THIS OUT??
  chr <- find.markerpos(cross, marker)[1,1]
  
  cross.imp <- fill.geno(cross) 
  marker.geno <- pull.geno(cross.imp, chr = chr)[,marker]
  
  if (!is.null(covar)){
    covlist <- str_split(covar, ",")[[1]]
  }
  
  formula <- paste0(pheno.col, " ~ ", 'marker.geno', 
                    ifelse(is.null(covar), "", " + "),
                    paste(covlist, collapse = " + "))
  
  if (is.null(covar)){
    lmdata <- data.frame(cross$pheno[,pheno.col], marker.geno)
    names(lmdata)[1] <- pheno.col # have to rename column 
  } else {
    lmdata <- data.frame(cross$pheno[,c(pheno.col,covlist)], marker.geno)
  }
  
  lmod <- lm(formula, lmdata)
  unadj_pval <- summary(lmod)$coefficients['marker.geno', 'Pr(>|t|)']
  return(unadj_pval)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# get_unadj_pval: this function gets the unadjusted p-value from the linear regression
#           between a particular phenotype and the genotype at a particular marker.
# Input: 
#     cross: 
#     pheno.col: 
#     marker:
#     covariates 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
get_var_expl <- function(cross, pheno.col, marker, covar = NULL){
  pheno <- pull.pheno(cross, pheno.col)
  chr <- find.markerpos(cross, marker)[1,1]
  
  cross.imp <- fill.geno(cross) 
  marker.geno <- pull.geno(cross.imp, chr = chr)[,marker]
  
  if (!is.null(covar)){
    covlist <- str_split(covar, ",")[[1]]
  }
  
  formula <- paste0(pheno.col, " ~ ", 'marker.geno', 
                    ifelse(is.null(covar), "", " + "),
                    paste(covlist, collapse = " + "))
  
  if (is.null(covar)){
    lmdata <- data.frame(cross$pheno[,pheno.col], marker.geno)
    names(lmdata)[1] <- pheno.col # have to rename column 
  } else {
    lmdata <- data.frame(cross$pheno[,c(pheno.col,covlist)], marker.geno)
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
#
# ### NEED TO UPDATE THIS CODE TO USE GEV AND PLOT SEPARATE LINES FOR X-CHR
# ### NOT IMPORTANT FOR PNUT PAPER BC GEV THRESHOLDS ARE NOT SIGNIFICANTLY 
# ### DIFFERENT EXCEPT IN THE CASE OF X-CHR-SPECIFIC PERMUTATIONS  
# 
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
      if (!has_xchr_perms(get(perm)) & ncol(mod)==3){ # no X-chr-specific perms, 1 lod column
        permgev <- paste0(perm, '.gev')
        perm.sum <- summary(get(permgev))
        abline(h = c(perm.sum[1], perm.sum[2]), lty=1:2)
      } else if (!has_xchr_perms(get(perm)) & ncol(mod)>3){ # no X-chr-specific perms, >1 lod column 
        perm.sum = summary(get(perm))
        abline(h = perm.sum[,1], lty=1:2) # lod.p.mu
        abline(h = perm.sum[,2], lty=1:2, col = 'blue') # lod.p 
        abline(h = perm.sum[,3], lty=1:2, col = 'red') # lod.mu
      } else if (has_xchr_perms(get(perm)) & ncol(mod)==3){ # X-chr-specific perms, 1 LOD column
        ### Should handle X-chr specific thresholds, but not necessary for virus research paper bc no X-chr QTL
        # perm.sum = summary(get(perm))
        # abline(h = c(perm.sum$A[1], perm.sum$A[2]), lty = 1:2)
        permgevA <- paste0(perm, 'A.gev')
        perm.sumA <- summary(get(permgevA))
        abline(h = c(perm.sumA[1], perm.sumA[2]), lty=1:2)
        # permgevX <- paste0(perm, 'X.gev')
      } else if (has_xchr_perms(get(perm)) & ncol(mod)>3){ # X-chr-specific perms, >1 LOD column 
        perm.sum = summary(get(perm))
        abline(h = perm.sum$A[,1], lty=1:2) # lod.p.mu
        abline(h = perm.sum$A[,2], lty=1:2, col = 'blue') # lod.p 
        abline(h = perm.sum$A[,3], lty=1:2, col = 'red') # lod.mu
      }
    }
    
    if (save == TRUE){ dev.off() } # close connection 
  }
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++
# plot_scans: use different thresholds for X-chr 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# plot_scans <- function(models, ylim = NULL){
#   for (i in 1:nrow(models)){
#     mod <- get(models[i, 'obj'])
#     
#     if (ncol(mod)>3) {lodcolumn=1:3} else {lodcolumn=1} # if ncol>3, multi-LOD model 
#     
#     if (is.null(ylim)){
#       plot(mod, ylab = "LOD", main = models[i,'name'], alternate.chrid = T, 
#            bandcol = "gray90", lodcolumn=lodcolumn)
#     } else {
#       plot(mod, ylab = "LOD", ylim = ylim, main = models[i,'name'], 
#            alternate.chrid = T, bandcol = "gray90", lodcolumn=lodcolumn)
#     }
#     
#     if (exists(models[i,'perm.obj'])){
#       perm.sum <- summary(get(models[i,'perm.obj']))
#       
#       if (length(perm.sum) == 10){ # separate thresholds for X chromosome 
#         ### CHANGE TO if (length(perm.sum == 2))
#         # Determine x-axis segments representing autosomes and X
#         temp <- mod
#         begend <- matrix(unlist(tapply(temp[,2],temp[,1],range)),ncol=2,byrow=TRUE)
#         rownames(begend) <- unique(mod[,1])
#         chr <- unique(as.character(mod[,1]))
#         begend <- begend[as.character(chr),,drop=FALSE]
#         len <- begend[,2]-begend[,1]
#         start <- c(0,cumsum(len))-c(begend[,1],0)
#         maxx <- sum(len)
#         
#         # end of chr1
#         x1=tail(mod[mod[,1] == chr[m],2], n=1)
#         
#         # beginning and end of chr2
#         
#         
#         # beginning and end of all chromosomes 
#         i <- c(1:20)
#         
#         segments(x0 = 0, x1 = 2600, y0 = perm.sum[[1]][1], lty=1)
#         segments(x0 = 0, x1 = 2600, y0 = perm.sum[[1]][2], lty=2)
#         segments(x0 = 2600, x1 = 3000, y0 = perm.sum[[2]][1], lty=1)
#         segments(x0 = 2600, x1 = 3000, y0 = perm.sum[[2]][2], lty=2)
#       } else { # only one set of thresholds 
#         ### REMOVE [[1]] from perm.sum calls 
#         abline(h = c(perm.sum[[1]][1], perm.sum[[1]][2]), lty=1:2)
#       }
#     }
#   }
# }


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
      
      # using built-in r/qtl function
      # plotPXG(cross.obj, marker = peaks[k,'marker'], 
      #         pheno.col = models[which(models$obj == peaks[k,'model']), 'colname'], 
      #         ylab = models[which(models$obj == peaks[k,'model']), 'abbr'], 
      #         main = paste(peaks[k,'marker'], " (chr", peaks[k,'chr'], ")", sep=""))
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





#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Adding support for subsetting 
# Removed code that colored imputed genotypes differently 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
pxg_subset <- function(cross, pheno, marker, cross.type, geno.map, xlab = 'Genotype', ylab, 
                title = NULL, bestfit = TRUE, groupcol, groupselect, theme = rmd_theme){
  # what chromosome is marker on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o]
  
  # linear regression (best-fit line only appropriate for autosomes) 
  if (chr == 'X'){bestfit = FALSE} 
  
  # plot title 
  if (is.null(title)){
    title = paste('Genotype vs. phenotype at chr', chr, 'marker', marker)
  }
  
  cross <- fill.geno(cross) # calculate imputed genotypes, just for plotting
  marker.genos <- cross$geno[[chr]]$data[,marker]
  marker.genos <- marker.genos[pull.pheno(cross, groupcol)==groupselect]
  
  if (bestfit == TRUE){
    fit <- lm(pheno ~ marker.genos)
    int <- summary(fit)$coefficients[1,1]
    slope <- summary(fit)$coefficients[2,1]
    slope.pval <- summary(fit)$coefficients[2,4]
  }
  
  df <- data.frame(geno = marker.genos, pheno = pheno)
  
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
  
  if (cross.type == 'bc'){ 
    df$geno <- mapvalues(df$geno, from = c(1,2), to = geno.labels)
  } else { # if cross.type = 'f2'
    if (chr != 'X'){
      df$geno <- mapvalues(df$geno, from = c(1,2,3), to = geno.labels) 
    } else { # X-chromosome marker 
      pgm <- pull.pheno(cross, "pgm")
      sex <- as.numeric(pull.pheno(cross, "sex") == "M")
      X.data <- data.frame(sex = sex, geno = geno, pgm = pgm)
      X.data$X.geno <- ifelse(sex==0 & geno==1 & pgm==0, 'AA', 
                              ifelse(sex==0 & geno==1 & pgm==1, 'BB', 
                                     ifelse(sex==0 & geno==2 & pgm==0, 'ABf', 
                                            ifelse(sex==0 & geno==2 & pgm==1, 'ABr', 
                                                   ifelse(sex==1 & geno==1, 'AY', 'BY')))))
      
      temp <- X.data$X.geno
      temp[which(is.na(df$geno))] <- NA
      df$geno <- temp
      
      df$geno <- mapvalues(df$geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                           to = geno.labels)
      }
  }
  
  df$geno <- factor(df$geno, levels = geno.labels)
  
  dfsum <- data.frame(geno = geno.labels,
                      mean = aggregate(df$pheno, by = list(df$imp_geno), FUN = mean, na.rm = TRUE)$x, 
                      sd = aggregate(df$pheno, by = list(df$imp_geno), FUN = sd, na.rm = TRUE)$x)
  
  p <- ggplot(data = df, mapping = aes(x = geno, y = pheno)) + 
      geom_jitter(width = 0.1, height = 0) + 
      geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean, ymax = mean), width = 0.2) + 
      geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean-sd, ymax = mean+sd), width = 0.3) +
      labs(title = title, y = ylab, x = xlab) + 
      {if (bestfit == TRUE) geom_abline(slope = slope, intercept = int, color = "red")} + 
      theme
  
  if (bestfit == TRUE){
    print(paste("slope p-value: ", slope.pval))
  }
  
  return(p)
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Adding support for multiple groups (e.g. infection, control) to pxg()
# Trying to make it so it works with and without groups, so there's only one function
# Maybe only have one pxg function, and that function calls pxg_group if there's a group 
# and pxg_ind if not? 
### 1. Works for markers with no missing genotypes, need to validate with markers that 
###    have missing genotypes. 
### 2. Need to validate with no grouping (groupcolumn = NULL); I know that df doesn't 
###    group column when groupcolumn=NULL
### 3. Need to validate with bestfit = FALSE (validate on X chr)
###    groupcolumn = 'infection'/NULL with bestfit= TRUE/FALSE 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
pxg_group <- function(cross, pheno, marker, cross.type, geno.map, xlab = 'Genotype', ylab, 
                plot.title = NULL, groupcolumn = NULL, groupcolors = NULL, bestfit = TRUE, 
                theme = rmd_theme){
  # what chromosome is marker on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o]
  
  # linear regression (best-fit line only appropriate for autosomes) 
  if (chr == 'X'){bestfit = FALSE} 
  
  if (is.null(plot.title)){
    plot.title = paste('Genotype vs. phenotype at chr', chr, 'marker', marker)
  }
  
  marker.genos <- cross$geno[[chr]]$data[,marker] # get genotypes at marker 
  which.missing <- which(is.na(marker.genos)) # are any genotypes missing at this marker? 
  cross.imp <- fill.geno(cross) # calculate imputed genotypes (for plotting)
  
  group <- cross$pheno[,groupcolumn] # category for grouping 
  grouplevels <- levels(group) # group names 
  num.groups <- length(grouplevels) # how many groups? 
  
  if (bestfit == TRUE){
    # allocate vectors for storing linear regression info 
    int <- vector(length = num.groups)
    slope <- vector(length = num.groups)
    slope.pval <- vector(length = num.groups)
    
    # compute intercept and slope for each group with linear reg 
    for (i in 1:num.groups){
      group.pheno <- pheno[which(group == grouplevels[i])]
      group.marker.genos <- marker.genos[which(group == grouplevels[i])]
      
      fit <- lm(group.pheno ~ group.marker.genos)
      int[i] <- summary(fit)$coefficients['(Intercept)', 'Estimate']
      slope[i] <- summary(fit)$coefficients['group.marker.genos', 'Estimate']
      slope.pval[i] <- summary(fit)$coefficients['group.marker.genos', 'Pr(>|t|)'] 
    }
  }
  
  df <- data.frame(geno = marker.genos, pheno = pheno, group = group,
                   imp_geno = cross.imp$geno[[chr]]$data[,marker])
  
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
  
  if (cross.type == 'bc'){ 
    df$geno <- mapvalues(df$geno, from = c(1,2), to = geno.labels)
    df$imp_geno <- mapvalues(df$imp_geno, from = c(1,2), to = geno.labels) 
  } else { # if cross.type = 'f2'
    if (chr != 'X'){
      df$geno <- mapvalues(df$geno, from = c(1,2,3), to = geno.labels)
      df$imp_geno <- mapvalues(df$imp_geno, from = c(1,2,3), to = geno.labels) 
    } else { # X-chromosome marker 
      pgm <- pull.pheno(cross, "pgm")
      sex <- as.numeric(pull.pheno(cross, "sex") == "M")
      geno <- df$imp_geno
      X.data <- data.frame(sex = sex, geno = geno, pgm = pgm)
      X.data$X.geno <- ifelse(sex==0 & geno==1 & pgm==0, 'AA', 
                       ifelse(sex==0 & geno==1 & pgm==1, 'BB', 
                       ifelse(sex==0 & geno==2 & pgm==0, 'ABf', 
                       ifelse(sex==0 & geno==2 & pgm==1, 'ABr', 
                       ifelse(sex==1 & geno==1, 'AY', 'BY')))))
      
      temp <- X.data$X.geno
      temp[which(is.na(df$geno))] <- NA
      df$geno <- temp
      
      df$imp_geno <- X.data$X.geno
      
      df$geno <- mapvalues(df$geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                           to = geno.labels)
      df$imp_geno <- mapvalues(df$imp_geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                               to = geno.labels) 
    }
  }
  
  df$geno <- factor(df$geno, levels = geno.labels)
  df$imp_geno <- factor(df$imp_geno, levels = geno.labels)
  
  dfsum <- data.frame(imp_geno = if(is.null(groupcolumn)){geno.labels} else {rep(geno.labels, num.groups)},
                      group = rep(grouplevels, each=length(levels(df$imp_geno))),
                      mean = aggregate(df$pheno, by = list(df$imp_geno, df$group), FUN = mean, na.rm = TRUE)$x, 
                      sd = aggregate(df$pheno, by = list(df$imp_geno, df$group), FUN = sd, na.rm = TRUE)$x)
  
  if (length(which.missing) > 0){ ### no idea if this is going to work
    p <- ggplot(data = df[-which.missing,], mapping = aes(x = geno, y = pheno, color = group)) +
      geom_jitter(width = 0.1, height = 0, shape = 20) + 
      geom_jitter(data = df[which.missing,], mapping = aes(x = imp_geno, y = pheno, color = scale_color_manual(values = saturation(groupcolors, 0.5))), width = 0.1) +
      {if (is.null(groupcolumn)) geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean, ymax = mean), width = 0.2)} +
      {if (is.null(groupcolumn)) geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean-sd, ymax = mean+sd), width = 0.3)} + 
      labs(title = plot.title, y = ylab, x = xlab) + 
      {if (!is.null(groupcolumn)) scale_color_manual(values = groupcolors)} +
      {if (bestfit == TRUE) geom_abline(slope = slope, intercept = int, color = groupcolors)} + 
      theme
  } else { 
    p <- ggplot(data = df, mapping = aes(x = geno, y = pheno, color = group)) + 
      geom_point(width = 0.1, height = 0, shape=20, position = position_jitterdodge(dodge.width = 0.3)) +
      #geom_line(aes(x = imp_geno, y = mean, group = group), data=dfsum) +
      #{if (is.null(groupcolumn)) geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean, ymax = mean), width = 0.2)} +
      {if (!is.null(groupcolumn)) geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean-sd, ymax = mean+sd), position = position_dodge(width=0.3), width = 0.3)} +
      labs(title = plot.title, y = ylab, x = xlab, color = groupcolumn) + 
      {if (!is.null(groupcolumn)) scale_color_manual(values = groupcolors)} +
      {if (bestfit == TRUE) geom_abline(slope = slope, intercept = int, color = groupcolors)} + 
      theme
  }
  
  if (bestfit == TRUE){
    print(paste("slope p-value: ", slope.pval))
  }
  
  return(p)
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Adding support for coloring dots by 2nd marker genotype
# Remove color-by-imputed-or-not because we're coloring by the genotype at a 2nd marker 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
pxg_colbym <- function(cross, pheno, marker, geno.map, xlab = 'Genotype', ylab, 
                title = NULL, bestfit = FALSE, theme = rmd_theme, marker2,
                type = 'scatter'){
  # cross type 
  cross.type <- class(cross)[1]
  
  # what chromosome is marker on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o]
  
  # what chromosome is marker2 on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker2)
  chr2 <- names(cross$geno)[o]
  
  cross <- fill.geno(cross)
  marker.genos <- cross$geno[[chr]]$data[,marker]
  marker2.genos <- cross$geno[[chr2]]$data[,marker2]
  
  df <- data.frame(geno = marker.genos,
                   geno2 = marker2.genos,
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
    df$geno2 <- mapvalues(df$geno2, from = c(1,2), to = geno.labels)
  } else { # if cross.type = 'f2'
    
    ### NEED TO VALIDATE THIS
    
    if (chr != 'X'){
      df$geno <- mapvalues(df$geno, from = c(1,2,3), to = geno.labels)
    } else { # X-chromosome marker 
      pgm <- pull.pheno(cross, "pgm")
      sex <- as.numeric(pull.pheno(cross, "sex") == "M")
      geno <- df$geno
      X.data <- data.frame(sex = sex, geno = geno, pgm = pgm)
      X.data$X.geno <- ifelse(sex==0 & geno==1 & pgm==0, 'AA', 
                              ifelse(sex==0 & geno==1 & pgm==1, 'BB', 
                                     ifelse(sex==0 & geno==2 & pgm==0, 'ABf', 
                                            ifelse(sex==0 & geno==2 & pgm==1, 'ABr', 
                                                   ifelse(sex==1 & geno==1, 'AY', 'BY')))))
      
      temp <- X.data$X.geno
      temp[which(is.na(df$geno))] <- NA
      df$geno <- temp
      
      df$geno <- mapvalues(df$geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                           to = geno.labels)
    }
  }
  
  df$geno <- factor(df$geno, levels = geno.labels)
  df$geno2 <- factor(df$geno2, levels = geno.labels)
  
  dfsum <- data.frame(geno = geno.labels,
                      mean = aggregate(df$pheno, by = list(df$geno), FUN = mean, na.rm = TRUE)$x, 
                      sd = aggregate(df$pheno, by = list(df$geno), FUN = sd, na.rm = TRUE)$x)
  
  if (type == 'scatter'){
    p <- ggplot(data = df, mapping = aes(x = geno, y = pheno, col = geno2)) + 
      geom_jitter(width = 0.1, height = 0) + 
      labs(title = title, y = ylab, x = xlab) + 
      scale_color_brewer(palette = 'Dark2', name = bquote(italic(Qpa2)~Genotype)) +
      theme
  } else if (type == 'boxplot'){
    p <- ggplot(data = df, mapping = aes(x = geno, y = pheno, fill = geno2)) + 
      geom_boxplot() +
      scale_fill_brewer(palette = 'Dark2', name = bquote(italic(Qpa2)~Genotype)) + 
      labs(title = title, y = ylab, x = xlab) +
      theme
  }
  
  
  if ((bestfit == TRUE) & (type == 'scatter')){ # only do best fit line on scatter plots 
    print(paste("slope p-value: ", slope.pval))
  }
  
  return(p)
}








#-------------------------------Statistical model------------------------------#
# Calculates dosage and adds to cross object 
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

# Calculates number of markers in a cross object 
calc.num.markers <- function(cross){
  p = 0
  for (c in 1:length(cross$geno)){
    p = p + ncol(cross$geno[[c]]$data)
  }
  return(p)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
### This function does mixed modeling (using Haley-Knott regression)
### To-do:
###   1. Figure out how to dynamically create random effect terms based on length(randoms).
###   For now, assume there are two (cage and batch) - line 44
###   2. Handle covariates in addition to random effects (also dynamically)
###   3. X-chr approach validated for BC but needs to be validated for F2 (i.e. compare 
###   with results from R/qtl) - line 24
###   4. Permutations 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
hk <- function(cross, pheno.col, covar = NULL){
  # Calculate genotype probabilities if necessary 
  if (!with(cross$geno[[1]], exists('prob'))){  #calc.genoprob hasn't been run
    cross <- calc.genoprob(cross)
    print("Running calc.genoprob...")
  }
  
  cross <- calc.dosage(cross) # calculate dosage (of B allele)
  
  p <- calc.num.markers(cross) # number of markers
  n = nrow(cross$pheno) # number of mice / observations
  
  # Pre-allocate vectors for results 
  markers <- vector(length = p)
  chrs <- vector(length = p)
  positions <- vector(length = p)
  lods <- vector(length = p)
  
  
  df <- covar # data frame for lm()
  df$pheno <- pull.pheno(cross, pheno.col) # add phenotype (outcome variable)
  
  # Create null model: just covariates 
  # formula0 <- paste0("pheno ~ ", paste(colnames(covar), collapse=" + "))
  # fit0 <- lm(formula0, df)
  
  # pheno <- cross$pheno[,pheno.col]
  # cage <- cross$pheno[,"Cage"]
  # batch <- cross$pheno[,"Batch"]
  
  marker.pos = 0
  for (c in 1:length(cross$geno)){ # for each chr
    for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker
      marker.pos = marker.pos + 1
      geno <- cross$geno[[c]]$dos[,k]
      
      #-----------------------------Create models------------------------------#
      ### For dynamically creating random effect terms
      #fit <- lmer(pheno ~ cross$geno[[c]]$dos[,k] +
      #              (1|cross$pheno[,randoms[1]]) + (1|cross$pheno[,randoms[2]]))
      
      fit <- lmer(pheno ~ geno + (1|batch/cage))
      
      # fit0 <- lm(pheno ~ sex + batch)
      # fit1 <- lm(pheno ~ sex + batch + geno)
      
      #----------------------------Calculate LODs------------------------------#
      # aov <- anova(fit1, fit0)
      # Fval <- aov$F[2] 
      # df <- abs(aov$Df[2])
      # lod <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      
      t <- summary(fit)$coefficients['geno','t value']
      Fval <- t^2
      df <- 1 # (in backcross, 2 in f2)
      lod <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      
      #-------------------------------Save LODs--------------------------------#
      marker.name <- colnames(cross$geno[[c]]$data)[k]
      markers[marker.pos] <- marker.name
      chrs[marker.pos] = c
      positions[marker.pos] = cross$geno[[c]]$map[marker.name]
      lods[marker.pos] = lod
    }
  }
  
  # Create object 
  model.df <- data.frame(chrs, positions, lods) 
  names(model.df) <- c("chr", "pos", "lod") # same column names as scanone object
  rownames(model.df) <- markers 
  class(model.df) <- c("scanone", "data.frame") # so we can use R/qtl functions (see Fig 5.4 in guide)
  
  return(model.df)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# This function does GxT modeling
# Input: 
#   cross*: an rqtl cross object with two infection groups (mock + infection) and 
#          a treatment covariate 
#   pheno.col*: name of phenotype
#   covar*: dataframe of covariates 
#   trt*: treatment to use in GxT term - must be a column in covar 
# Output: scanone object with three LOD columns  
### TO-DO: 
###   1. Handle variance heterogeneity 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
gxt <- function(cross, pheno.col, covar, trt){
  # Calculate genotype probabilities if necessary
  if (!with(cross$geno[[1]], exists('prob'))){  # calc.genoprob hasn't been run
    cross <- calc.genoprob(cross)
    print("Running calc.genoprob...")
  }

  cross <- calc.dosage(cross) # calculate dosage (of B allele)

  p <- calc.num.markers(cross) # number of markers
  n <- nrow(cross$pheno) # number of mice / observations

  # Pre-allocate vectors for results
  markers <- vector(length = p)
  chrs <- vector(length = p)
  positions <- vector(length = p)
  lods1 <- vector(length = p)
  lods2 <- vector(length = p)
  lods3 <- vector(length = p)

  lmdata <- covar # data frame for lm()
  lmdata$pheno <- pull.pheno(cross, pheno.col) # add phenotype (outcome variable)
  lmdata$geno <- rep(NA, nrow(covar)) # to hold genotypes

  # Create null model: just covariates
  formula0 <- paste0("pheno ~ ", paste(colnames(covar), collapse=" + "))
  fit0 <- lm(formula0, lmdata)

  marker.pos = 0
  for (c in 1:length(cross$geno)){ # for each chr
    for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker
      lmdata$geno <- NA # clear geno col
      marker.pos = marker.pos + 1
      geno <- cross$geno[[c]]$dos[,k]
      lmdata$geno <- geno # update geno col with marker genotype data

      #-----------------------------Create models------------------------------#
      formula1 <- paste0(formula0, " + geno")
      fit1 <- lm(formula1, lmdata)

      formula2 <- paste0(formula1, " + ", paste0("geno*", trt))
      fit2 <- lm(formula2, lmdata)

      #----------------------------Calculate LODs------------------------------#
      # Marginal effect QTL (fit1 vs fit0)
      aov <- anova(fit0, fit1)
      Fval <- aov$F[2]
      df <- abs(aov$Df[2]) # abs() so that this works regardless of direction models are specified in
      lod1 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)

      # GxT QTL (fit2 vs fit1)
      aov <- anova(fit1, fit2)
      Fval <- aov$F[2]
      df <- abs(aov$Df[2])
      lod2 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)

      # Marginal, GxT or some combo of both (fit2 vs fit0). Not as powerful at
      # detecting purely marginal or purely interactive QTL, but will be more powerful
      # for QTL having both a marginal and an interactive component.
      aov <- anova(fit0, fit2)
      Fval <- aov$F[2]
      df <- abs(aov$Df[2])
      lod3 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)

      #------------------------------Save LODs---------------------------------#
      marker.name <- colnames(cross$geno[[c]]$data)[k]
      markers[marker.pos] <- marker.name
      chrs[marker.pos] = c
      positions[marker.pos] = cross$geno[[c]]$map[marker.name]

      lods1[marker.pos] = lod1
      lods2[marker.pos] = lod2
      lods3[marker.pos] = lod3
    }
  }

  # Create object
  model.df <- data.frame(chrs, positions, lods1, lods2, lods3)

  # Names for 2part model: lod.p.mu, lod.p, lod.mu
  names(model.df) <- c("chr", "pos", "lod1", "lod2", "lod3")

  rownames(model.df) <- markers
  class(model.df) <- c("scanone", "data.frame") # so we can use R/qtl functions (see Fig 5.4 in guide)

  return(model.df)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++
# This function does GxT modeling handling variance heterogeneity with a dglm 
# Input: 
#   cross*: an rqtl cross object with two infection groups (mock + infection) and 
#          a treatment covariate 
#   pheno.col*: name of phenotype
#   covar*: dataframe of covariates 
#   trt*: treatment to use in GxT term - must be a column in covar 
#   ### Note - handles variance heterogeneity between groups 
# Output: scanone object with three LOD columns  
#+++++++++++++++++++++++++++++++++++++++++++++++++++
gxt.dglm <- function(cross, pheno.col, covar, trt){
  # Calculate genotype probabilities if necessary
  # if (!with(cross$geno[[1]], exists('prob'))){  # calc.genoprob hasn't been run
  #   cross <- calc.genoprob(cross)
  #   print("Running calc.genoprob...")
  # }
  #cross <- calc.dosage(cross) # calculate dosage (of B allele)
  
  p <- calc.num.markers(cross) # number of markers
  n <- nrow(cross$pheno) # number of mice / observations
  
  # Pre-allocate vectors for results
  markers <- vector(length = p)
  chrs <- vector(length = p)
  positions <- vector(length = p)
  lods1 <- vector(length = p) # marginal QTL
  lods2 <- vector(length = p) # GxT QTL
  lods3 <- vector(length = p) # marginal + GxT QTL 
  lods4 <- vector(length = p) # vQTL
  
  lmdata <- covar # data frame for lm()
  lmdata$pheno <- pull.pheno(cross, pheno.col) # add phenotype (outcome variable)
  lmdata$geno <- rep(NA, nrow(covar)) # to hold genotypes
  
  # Create null model: just covariates
  formula0str <- paste0("pheno ~ ", paste(colnames(covar), collapse=" + "))
  formula0 <- as.formula(formula0str)
  dformulastr <- paste0("~", trt)
  dformula <- as.formula(dformulastr)
  fit0 <- dglm(formula0, as.formula(dformula), data = lmdata, method = 'ml')
  
  marker.pos = 0
  for (c in 1:length(cross$geno)){ # for each chr
    for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker
      lmdata$geno <- NA # clear geno col
      marker.pos = marker.pos + 1
      
      #geno <- cross$geno[[c]]$dos[,k]
      geno <- as.factor(cross$geno[[c]]$data[,k])
      lmdata$geno <- geno # update geno col with marker genotype data
      
      #-----------------------------Create models------------------------------#
      formula1str <- paste0(formula0str, " + geno")
      formula1 <- as.formula(formula1str)
      fit1 <- dglm(formula1, dformula, data = lmdata, method = 'ml')
      
      formula2str <- paste0(formula1str, " + ", paste0("geno*", trt))
      formula2 <- as.formula(formula2str)
      fit2 <- dglm(formula2, dformula, data = lmdata, method = 'ml')
      
      dformulavstr <- paste0(dformulastr, " + geno")
      dformulav <- as.formula(dformulastr)
      fitv <- dglm(formula1, dformulav, data = lmdata, method = 'ml')
      
      #----------------------------Calculate -logP-----------------------------#
      # Marginal effect QTL (fit1 vs fit0)
      aov <- anova(fit0, fit1)
      # Fval <- aov$F[2]
      # df <- abs(aov$Df[2]) # abs() so that this works regardless of direction models are specified in
      # lod1 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      pval1 <- -log10(aov$Seq.P[1])
      
      # GxT QTL (fit2 vs fit1)
      aov <- anova(fit1, fit2)
      # Fval <- aov$F[2]
      # df <- abs(aov$Df[2])
      # lod2 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      pval2 <- -log10(aov$Seq.P[1])
      
      # Marginal, GxT or some combo of both (fit2 vs fit0). Not as powerful at
      # detecting purely marginal or purely interactive QTL, but will be more powerful
      # for QTL having both a marginal and an interactive component.
      aov <- anova(fit0, fit2)
      # Fval <- aov$F[2]
      # df <- abs(aov$Df[2])
      # lod3 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      pval3 <- -log10(aov$Seq.P[1])
      
      # vQTL 
      aov <- anova(fit1, fitv)
      pval4 <- -log10(aov$Seq.P[1])
      
      #------------------------------Save LODs---------------------------------#
      marker.name <- colnames(cross$geno[[c]]$data)[k]
      markers[marker.pos] <- marker.name
      chrs[marker.pos] = c
      positions[marker.pos] = cross$geno[[c]]$map[marker.name]
      
      lods1[marker.pos] = pval1
      lods2[marker.pos] = pval2
      lods3[marker.pos] = pval3
      lods4[marker.pos] = pval4
    }
  }
  
  # Create object
  model.df <- data.frame(chrs, positions, lods1, lods2, lods3, lods4)
  
  # Names for 2part model: lod.p.mu, lod.p, lod.mu
  names(model.df) <- c("chr", "pos", "lod1", "lod2", "lod3", "lod4")
  
  rownames(model.df) <- markers
  class(model.df) <- c("scanone", "data.frame") # so we can use R/qtl functions (see Fig 5.4 in guide)
  
  return(model.df)
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++
# This function does GxT modeling with random effect modeling 
# Input: 
#   cross*: an rqtl cross object with two infection groups (mock + infection) and 
#          a treatment covariate 
#   pheno.col*: name of phenotype
#   covar*: dataframe of covariates 
#   trt*: treatment to use in GxT term - must be a column in covar 
#   randoms: dataframe of covariates to model as random effects 
# Output: list of three scanone objects 
### TO-DO: 
###   1. Handle variance heterogeneity 
###   2. Handle random effect covariates? 
### NOT SURE HOW TO COMPARE FIT1 TO FIT0 WHEN BOTH HAVE RANDOM EFFECTS
#+++++++++++++++++++++++++++++++++++++++++++++++++++
gxt_rand <- function(cross, pheno.col, covar, trt, randoms = NULL){
  # Calculate genotype probabilities if necessary 
  if (!with(cross$geno[[1]], exists('prob'))){  # calc.genoprob hasn't been run
    cross <- calc.genoprob(cross)
    print("Running calc.genoprob...")
  }
  
  cross <- calc.dosage(cross) # calculate dosage (of B allele)
  
  p <- calc.num.markers(cross) # number of markers
  n = nrow(cross$pheno) # Number of mice / observations
  
  # Pre-allocate vectors for results 
  markers <- vector(length = p)
  chrs <- vector(length = p)
  positions <- vector(length = p)
  lods1 <- vector(length = p)
  lods2 <- vector(length = p)
  lods3 <- vector(length = p)
  
  lmdata <- cbind(covar,randoms) # data frame for lm()
  lmdata$pheno <- pull.pheno(cross, pheno.col) # add phenotype (outcome variable)
  lmdata$geno <- rep(NA, nrow(covar)) # to hold genotypes 

  # Create null model: just covariates 
  fixedcovar <- paste(colnames(covar), collapse=" + ")
  randcovar <- paste0("(1|", colnames(randoms), ")", collapse=" + ")
  formula0 <- paste0("pheno ~ ", fixedcovar, " + ", randcovar)
  fit0 <- lmer(formula0, lmdata)
  
  marker.pos = 0
  for (c in 1:length(cross$geno)){ # for each chr 
    for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker 
      lmdata$geno <- NA # clear geno col 
      marker.pos = marker.pos + 1
      geno <- cross$geno[[c]]$dos[,k]
      lmdata$geno <- geno # update geno col with marker genotype data 
      
      #-----------------------------Create models------------------------------#
      formula1 <- paste0(formula0, " + geno")
      fit1 <- lmer(formula1, lmdata)
      
      formula2 <- paste0(formula1, " + ", paste0("geno*", trt))
      fit2 <- lmer(formula2, lmdata)
      
      #----------------------------Calculate LODs------------------------------#
      ### Need to use MLE instead of LMER with ANOVA on LMMs 
      ### ANOVA gives you a Chisq value (not an F val) and p-val
      
      # Marginal effect QTL (fit1 vs fit0)
      aov <- anova(fit0, fit1)
      Fval <- aov$F[2] 
      df <- abs(aov$Df[2]) # abs() so that this works regardless of direction models are specified in 
      lod1 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      
      # GxT QTL (fit2 vs fit1)
      aov <- anova(fit1, fit2)
      Fval <- aov$F[2] 
      df <- abs(aov$Df[2])
      lod2 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      
      # Marginal, GxT or some combo of both (fit2 vs fit0). Not as powerful at 
      # detecting purely marginal or purely interactive QTL, but will be more powerful
      # for QTL having both a marginal and an interactive component. 
      aov <- anova(fit0, fit2)
      Fval <- aov$F[2] 
      df <- abs(aov$Df[2])
      lod3 <- (n/2)*log10(Fval*(df/(n-df-1)) + 1)
      
      #------------------------------Save LODs---------------------------------#
      marker.name <- colnames(cross$geno[[c]]$data)[k]
      markers[marker.pos] <- marker.name
      chrs[marker.pos] = c
      positions[marker.pos] = cross$geno[[c]]$map[marker.name]
      
      lods1[marker.pos] = lod1
      lods2[marker.pos] = lod2
      lods3[marker.pos] = lod3
    }
  }
  
  # Create object 
  model.df <- data.frame(chrs, positions, lods1, lods2, lods3) 
  
  # Names for 2part model: lod.p.mu, lod.p, lod.mu
  names(model.df) <- c("chr", "pos", "lod1", "lod2", "lod3") 
  
  rownames(model.df) <- markers 
  class(model.df) <- c("scanone", "data.frame") # so we can use R/qtl functions (see Fig 5.4 in guide)
  
  return(model.df)
}


#-------------------------------Permutation tests------------------------------#
#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to perform permutation tests more efficiently than R/qtl. Instead 
# of performing n.perm tests from start to finish and logging the highest LOD, 
# perform n.perm tests at each marker. 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
scanone_perm <- function(cross, pheno.col, addcovar, n.perm){
  # create 1000 permutations of the phenotype, add to matrix 
  
  # lm multi response, add to matrix of results 
  
  # find highest LOD from each genome scan and create distribution 
  
  # return - vector of 1000 LOD scores with some attributes 
  # method = "em"
  # model = "normal" 
  # type = "bc"
  # class = scanoneperm, matrix 
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Permutation tests for GxT 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
gxt_perm <- function(cross, pheno.col, perm.strata=NULL, covar, n.perm, trt){
  
  pheno.og <- pull.pheno(combined, pheno.col)
  # remove NAs bc don't want to permute those 
  # use !is.na(pheno.og) to filter genotypes later on 
  pheno <- pheno.og[!is.na(pheno.og)]
  
  # create 1000 permutations of phenotype 
  pheno.perms <- matrix(nrow = length(pheno), ncol = n.perm)
  for (c in 1:n.perm){
    pheno.perms[,c] <- pheno[sample(length(pheno))]
  }
  
  # calculate dosage for hk regression 
  cross <- calc.dosage(cross)
  
  # at each marker, run lm.multiresponse(), add to matrix of results 
  # create empty matrix for storing results 
  permLODs <- matrix(nrow = n.perm, ncol = calc.num.markers(cross))
  
  # For reference: 
  # lm.multiresponse(formula, response.matrix, data, null.formula  = NULL,
  #                  null.fit = NULL, rsquared = FALSE, pvalue = FALSE,
  #                  logP = FALSE, LOD = FALSE, total.ss = FALSE,
  #                  weights = rep(1, nrow(data)), model.args = list(),
  #                  verbose.at.every = 0)
  # 
  ### Do we need this dataframe? I think this should be passed in as 'data' 
  lmdata <- covar[!is.na(pheno.og),] # data frame for lm()
  lmdata$geno <- rep(NA, nrow(lmdata)) # to hold genotypes
  
  # Create null model: just covariates
  formula0 <- paste0("pheno ~ ", paste(colnames(covar), collapse=" + "))
  fit0 <- lm(formula0, lmdata)
  
  marker.pos = 0
  for (c in 1:length(cross$geno)){ # for each chr
    for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker
      lmdata$geno <- NA # clear geno col
      marker.pos = marker.pos + 1
      
      geno <- cross$geno[[c]]$dos[,k]
      geno <- geno[!is.na(pheno.og)] # remove genotypes for mice that don't have the phenotype 
      
      lmdata$geno <- geno # update geno col with marker genotype data
      
      #-----------------------------Create models------------------------------#
      formula1 <- paste0(formula0, " + geno")
      fit1 <- lm(formula1, lmdata)
      
      formula2 <- paste0(formula1, " + ", paste0("geno*", trt))
      fit2 <- lm(formula2, lmdata)
      
      #----------------------------Calculate LODs------------------------------#
      res <- lm.multiresponse(formula = formula2, 
                       response.matrix = pheno.perms, 
                       data = lmdata, 
                       null.formula = formula1, 
                       null.fit = fit1, 
                       LOD = TRUE)
      
      #------------------------------Save LODs---------------------------------#
      resLODs <- res$LOD
      permLODs[,marker.pos] <- resLODs
    }
  }
  
  # find highest LOD from each genome scan and create distribution 
  maxLODs <- apply(permLODs, 1, max)
  
  # return - vector of 1000 LOD scores with some attributes 
  attributes(maxLODs) <- list(method = 'hk', model = 'normal',
                              type = 'f2', class = c('scanoneperm', 'matrix'))
  
  return(maxLODs) 
}










#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Function to plot phenotype x genotype plots with ggplot
### THIS VERSION OF PXG() PLOTS IMPUTED GENOTYPES AS GRAY AND DIRECTLY GENOTYPED
### GENOTYPES AS BLACK. NO LONGER USED 
# Alternative to qtl::plotPXG, uses ggplot and more intuitive error bars  
# Input:
#     cross: r/qtl cross object 
#     pheno: vector of phenotype values for all mice 
#     marker: marker whose genotype to plot
#     geno.map: list mapping A and B allele to more informative labels 
#     xlab: x axis label
#     ylab: y axis label 
#     title: title
#     theme: ggplot theme
#     type: scatter (default), boxplot or violin 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
pxg_old <- function(cross, pheno, marker, geno.map, xlab = NULL, ylab, 
                title = NULL, bestfit = TRUE, theme = rmd_theme, type = 'scatter', ...){
  # cross type 
  cross.type <- class(cross)[1]
  
  # what chromosome is marker on?
  o <- sapply(cross$geno, function(a, b) b %in% colnames(a$data), marker)
  chr <- names(cross$geno)[o]
  
  # linear regression (best-fit line only appropriate for autosomes) 
  if ((chr == 'X')&('M' %in% levels(cross$pheno$sex))){bestfit = FALSE} 
  
  # plot title 
  # if (is.null(title)){
  #   title = paste('Genotype vs. phenotype at chr', chr, 'marker', marker)
  # }
  
  if(is.null(xlab)){
    qtl_name <- qtl_map$qtl_name[qtl_map$chr == chr]
    xlab <- as.expression(bquote(italic(.(qtl_name))*" (chr"*.(chr)*") "*Genotype))
  }
  
  marker.genos <- cross$geno[[chr]]$data[,marker]
  which.missing <- which(is.na(marker.genos))
  # calculate imputed genotypes, for plotting and regression. Need to keep separate
  # from observed genotypes for plotting (imputed genotypes gray, observed genotypes black)
  cross.imp <- fill.geno(cross) 
  imp.marker.genos <- cross.imp$geno[[chr]]$data[,marker]
  
  if (bestfit == TRUE){
    fit <- lm(pheno ~ imp.marker.genos)
    int <- summary(fit)$coefficients[1,1]
    slope <- summary(fit)$coefficients[2,1]
    slope.pval <- summary(fit)$coefficients[2,4]
  }
  
  df <- data.frame(geno = marker.genos,
                   imp_geno = imp.marker.genos,
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
    df$imp_geno <- mapvalues(df$imp_geno, from = c(1,2), to = geno.labels) 
  } else { # if cross.type = 'f2'
    if (chr != 'X'){
      df$geno <- mapvalues(df$geno, from = c(1,2,3), to = geno.labels)
      df$imp_geno <- mapvalues(df$imp_geno, from = c(1,2,3), to = geno.labels) 
    } else { # X-chromosome marker 
      pgm <- pull.pheno(cross, "pgm")
      sex <- as.numeric(pull.pheno(cross, "sex") == "M")
      geno <- df$imp_geno
      X.data <- data.frame(sex = sex, geno = geno, pgm = pgm)
      X.data$X.geno <- ifelse(sex==0 & geno==1 & pgm==0, 'AA', 
                              ifelse(sex==0 & geno==1 & pgm==1, 'BB', 
                                     ifelse(sex==0 & geno==2 & pgm==0, 'ABf', 
                                            ifelse(sex==0 & geno==2 & pgm==1, 'ABr', 
                                                   ifelse(sex==1 & geno==1, 'AY', 'BY')))))
      
      temp <- X.data$X.geno
      temp[which(is.na(df$geno))] <- NA
      df$geno <- temp
      
      df$imp_geno <- X.data$X.geno
      
      df$geno <- mapvalues(df$geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                           to = geno.labels)
      df$imp_geno <- mapvalues(df$imp_geno, from = c('AA', 'ABf', 'ABr', 'BB', 'AY', 'BY'), 
                               to = geno.labels) 
    }
  }
  
  df$geno <- factor(df$geno, levels = geno.labels)
  df$imp_geno <- factor(df$imp_geno, levels = geno.labels)
  
  dfsum <- data.frame(imp_geno = geno.labels,
                      mean = aggregate(df$pheno, by = list(df$imp_geno), FUN = mean, na.rm = TRUE)$x, 
                      sd = aggregate(df$pheno, by = list(df$imp_geno), FUN = sd, na.rm = TRUE)$x)
  
  if (type == 'scatter'){
    if (length(which.missing) > 0){
      p <- ggplot(data = df[-which.missing,], mapping = aes(x = geno, y = pheno)) +
        geom_jitter(width = 0.1, height = 0) + scale_x_discrete(drop = FALSE) + 
        geom_jitter(data = df[which.missing,], mapping = aes(x = imp_geno, y = pheno), color = "gray", width = 0.1) +
        geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean, ymax = mean), width = 0.2) +
        geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean-sd, ymax = mean+sd), width = 0.3) + 
        labs(title = title, y = ylab, x = xlab) + 
        {if (bestfit == TRUE) geom_abline(slope = slope, intercept = int, color = "red")} + 
        theme
    } else { 
      p <- ggplot(data = df, mapping = aes(x = geno, y = pheno)) + 
        geom_jitter(width = 0.1, height = 0) + 
        geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean, ymax = mean), width = 0.2) + 
        geom_errorbar(data = dfsum, mapping = aes(x = imp_geno, y = mean, ymin = mean-sd, ymax = mean+sd), width = 0.3) +
        labs(title = title, y = ylab, x = xlab) + 
        {if (bestfit == TRUE) geom_abline(slope = slope, intercept = int, color = "red")} + 
        theme
    }
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



