#!/usr/bin/env RScript
# Title: [ENTER TITLE]
# Authors: [ENTER AUTHORS]
# Brief description: [ENTER DESCRIPTION]

# Initialize workspace -----------------------------------------------------------------------
rm(list=ls())
library(data.table)
library(ggplot2)
source('~/theme_alex.R')
library(stringr)

# Attempt to identify contaminants based on some criteria and make a report
# Report as a TAB delimited file where eah row is a read
# Colums:   [1] read_id, 
#           [2] likely species, 
#           [3] taxid, 
#           [4] molecule length,
#           [5] mismatched Ts, 
#           [6] matched Cs,
#           [7] belonging to known contaminant list, 
#           [8] is this read belonging to a species with high Cs,
#           [9] probability


# Command line arguments ---------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

pe_tblat = args[[1]] 
output = args[[2]]
LUT = args[[3]]

#pe_tblat = '/workdir/apc88/WGBS_pipeline/sample_output/filter_unconverted/grammy/MET-1-01.tblat.1'
#LUT = '/workdir/apc88/WGBS_pipeline/LUTGrammy/taxids_names_lengths_tax.tab'

# Functions ----------------------------------------------------------------------------------
quantify_chars <- function(mismatch_string, n){
  return(str_count(mismatch_string, n))
}

assign_unknown_read <- function(val, pass_prob){
  if( is.na(val) ){
    val <- sample(x=c(TRUE, FALSE), size = 1, prob = c(pass_prob, 1-pass_prob))
  }
  return(val)
}

evaluate_species <- function(df1){
  df1$pass <- NA
  df1$pass[df1$C_count>5*df1$T_count]<-FALSE
  df1$pass[df1$T_count>5*df1$C_count]<-TRUE
  
  if ( !all(is.na(df1$pass)) ){ #randomly assign reads probabilistically
    pass_prob <- sum(df1$pass[!is.na(df1$pass)])/length(df1$pass[!is.na(df1$pass)])
    df1$pass <- sapply(X=df1$pass, FUN=assign_unknown_read, pass_prob)
  }
  
  # species with uninformative C regions are "passed" but can fail later
  df1$pass[is.na(df1$pass)]<-TRUE
  return(df1)
}

# Process files ------------------------------------------------------------------------------
pe_tblat = fread(pe_tblat)
colnames(pe_tblat)[[1]]<-'Taxid'
LUT = fread(LUT, fill=TRUE)

df <- merge(pe_tblat, LUT[, c('Taxid', 'Name', 'superkingdom', 'Length', 'species')])
df$C_count <- sapply(X = df$mismatch_string, FUN = quantify_chars, n='C')
df$T_count <- sapply(X = df$mismatch_string, FUN = quantify_chars, n='T')

#handle each species seperately
df_list <- split(df, df$species)

ark <- rbindlist(lapply(df_list, evaluate_species))
ark$pass <- NULL #for now
fwrite(x = ark, file = output, sep = '\t', quote = FALSE, col.names = TRUE)



