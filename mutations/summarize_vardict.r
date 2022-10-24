#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author  : zhangsiwen
# @ChangeLog
#     20221018, first version
#     20221024, change StatsByKey(), add GetKeyContent() and CountByKey(),  add --mode, add fixed parameters


usage <- "
  Usage:
    Rscript --vanilla summarize_vardict.r --input vardict_output [--out outprefix --dir ./ --mode tsv|xlsx|both]
  
  Arguments:
    --input     vardict output
    --out       output name without path, default is same as input name
    --dir       output path, default is same as input file
    --mode      output mode, default is 'both' 
                    if choose 'both', both tsv and xlsx files will be generated,
                    if choose 'tsv', 18 files will be generated,
                    if choose 'xlsx', only 1 xlsx file will be generated.
    
  Vardict output:
    https://github.com/AstraZeneca-NGS/VarDictJava#output-columns
  
"

# test if there is at least one argument: if not, return an error
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop(usage, call.=FALSE)
} else if (args[1]=='-h'||args[1]=='--help') {
  stop(usage, call.=FALSE)
}


library(optparse)
library(data.table)
library(dplyr)
library(xlsx)
library(gtools)
#library(ggplot2)
#theme_set(theme_classic())
#library(XLConnect)


# fixed parameters
if ( TRUE ) {
  keys <- c('MutationType', 'VAFRange', 'ContextType')
  columns <- c('PMean', 'QMean', 'NM')
  cutoffs <- c(20, 20, 3)
  bases <- c('A','C','G','T')
  bases_perm <- permutations(4, 2, bases)
  mutation_types <- paste0(bases_perm[,1], '>', bases_perm[,2])
  context_df <- expand.grid(permutations(4, 1, bases), mutation_types, permutations(4, 1, bases))
  context_types <- sort(paste0(context_df$Var1, '(', context_df$Var2, ')', context_df$Var3))
  vafrange_types <- c('[0,0.0001)', '[0.0001,0.0005)', '[0.0005,0.001)', '[0.001,0.01)')
}



# read file and process
ReadVardict <- function ( file ) {
  
  # read and check header
  df <- fread(file)
  if (colnames(df)[1] != "Sample") {
    if (dim(df)[2]==36) {
      colnames(df) <- c('Sample','Gene','Chr','Start','End','Ref','Alt','Depth','AltDepth','RefFwdReads','RefRevReads','AltFwdReads','AltRevReads','Genotype','AF','Bias','PMean','PStd','QMean','QStd','MAPQ','QRATIO','HIFREQ','EXTRAFR','SHIFT3','MSI','MSILEN','NM','HICNT','HICOV','5pFlankSeq','3pFlankSeq','SEGMENT:CHR_START_END','VARTYPE','DUPRATE','SV')
    } else if (dim(df)[2]==37) {
      colnames(df) <- c('Sample','Gene','Chr','Start','End','Ref','Alt','Depth','AltDepth','RefFwdReads','RefRevReads','AltFwdReads','AltRevReads','Genotype','AF','Bias','PMean','PStd','QMean','QStd','MAPQ','QRATIO','HIFREQ','EXTRAFR','SHIFT3','MSI','MSILEN','NM','HICNT','HICOV','5pFlankSeq','3pFlankSeq','SEGMENT:CHR_START_END','VARTYPE','DUPRATE','SV', 'CRISPR')
    } else {
      stop("Only accept vardict result, see more information at:
           https://github.com/AstraZeneca-NGS/VarDictJava#output-columns")
    }
  }
  
  # filter mutation and frequency
  df <- df[df$Ref!=df$Alt & df$AF<0.01 & df$VARTYPE=='SNV', ]
  df$AF <- df$AltDepth/df$Depth
  
  # test strand bias
  # for (i in 1:dim(df)[1]) {
  #   mat <- matrix(unlist(data.frame(c(df[i,'RefFwdReads'],df[i,'AltFwdReads'],df[i,'RefRevReads'],df[i,'AltRevReads']))),2)
  #   df[i,'FisherTestStrandBias'] <- fisher.test(mat)$p.value
  # }
  
  # mutation type
  df$MutationType <- paste0(df$Ref,'>',df$Alt)
  
  # context type
  df$Upstream <- substr(df$`5pFlankSeq`, nchar(df$`5pFlankSeq`), nchar(df$`5pFlankSeq`))
  df$Downstream <- substr(df$`3pFlankSeq`, 1, 1)
  df$ContextType <- paste0(df$Upstream,'(',df$MutationType,')',df$Downstream)
  
  # VAF ranges
  df$VAFRange <- ifelse(df$AF<0.0001, '[0,0.0001)',
                      ifelse(df$AF<0.0005, '[0.0001,0.0005)', 
                             ifelse(df$AF<0.001, '[0.0005,0.001)', '[0.001,0.01)')))
  
  return(data.frame(df))
  
}



# get key content from key name
GetKeyContent <- function ( key ) {
  if ( key == 'MutationType' ) {
    key_types <- mutation_types
  } else if ( key == 'ContextType' ) {
    key_types <- context_types
  } else if ( key == 'VAFRange' ) {
    key_types <- vafrange_types
  }
  return(key_types)
}



# count number of key content
CountByKey <- function (df, key='MutationType') {
  
  # count for mutation types
  key_types <- GetKeyContent(key)
  key_count <- data.frame(key_types)
  tmp <- data.frame(table(df[,key]))
  key_count <- merge(key_count, tmp, by.x = "key_types", by.y = "Var1", all.x=TRUE)
  key_count[is.na(key_count )] <- 0
  colnames(key_count) <- c(key, 'Count')
  rownames(key_count) <- key_types
  
  # ratio and rank
  total <- sum(key_count$Count)
  key_count$Ratio <- paste0(round(100*key_count$Count/total,4), '%')
  
  return(key_count)
  
}



# test strand bias
StrandBiasSummary <- function ( df, key='MutationType' ) {
  
  # calculate p.value
  if (! 'FisherTestStrandBias' %in% colnames(df)) {
    for (i in 1:dim(df)[1]) {
      mat <- matrix(unlist(data.frame(c(df[i,'RefFwdReads'],df[i,'AltFwdReads'],df[i,'RefRevReads'],df[i,'AltRevReads']))),2)
      df[i,'FisherTestStrandBias'] <- fisher.test(mat)$p.value
    }    
  }
  
  # single key type
  key_count <- CountByKey(df, key=key) 
  bias_count <- CountByKey(df[df$FisherTestStrandBias<=0.05,], key=key)
  colnames(bias_count) <- c(key, 'BiasCount')
  
  bias_table <- merge(key_count[,-3], bias_count[,-3], by = key)
  bias_table$BiasCountRatio <- paste0(round(100*bias_table$BiasCount/bias_table$Count,4), '%')
  
  key_total <- sum(bias_table$Count)
  bias_table$BiasCountRatioInTotalCount <- paste0(round(100*bias_table$BiasCount/key_total,4), '%')
  
  return(bias_table)
  
}



# summary stats
SummarizeStats <- function ( df, id, vars ) {
  
  mat <- melt(df, id.vars=id, measure.vars=vars)
  stats <- mat[,c(1,3)] %>% 
    group_by(across(all_of(id))) %>% 
    summarise(Min = min(value),
              Q1 = quantile(value,probs=0.25),
              Median = median(value),
              Mean = mean(value),
              Q3 = quantile(value,probs=0.75),
              Max = max(value),
              LowerBoundary = quantile(value,probs=0.25) - 1.5*(quantile(value,probs=0.75)-quantile(value,probs=0.25)),
              UpperBoundary = quantile(value,probs=0.75) + 1.5*(quantile(value,probs=0.75)-quantile(value,probs=0.25)),
              )
  
  stats <- data.frame(stats)
  rownames(stats) <- stats[,1]
  stats[,1] <- NULL
  
  if (max(stats$Max)>=1) {
    stats <- stats %>% mutate_if(is.numeric, round)
  } else {
    stats <- stats %>% mutate_if(is.numeric, round, 8)
  }
  
  return(stats)
  
}



# draw plot, not used
DrawPlot <- function ( df, x, y ) {

  g <- ggplot(df, aes_string(x=reorder(x,-y,na.rm = TRUE), y=y, color=x, shape=x)) +
    geom_boxplot(alpha=0.5, outlier.color=NA) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    ylab(y) + xlab(x) +
    theme_classic(base_size = 14) +
    theme(axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          strip.text=element_text(size = 14))

  return(g)
  
}



# summarize stats by specific key
StatsByKey <- function ( df, key='MutationType', condition='NA' ) {
  
  # count for mutation types
  key_count <- CountByKey(df, key=key) 

  # quantiles
  key_stats <- SummarizeStats(df, key, 'AF')
  key_df <- merge(as.data.frame(key_count), as.data.frame(key_stats), by="row.names", all.x=TRUE)[,-1]
  
  key_df$Condition <- condition
  return(key_df)
  
}



# wilcox test for multiple dataframes
TestDataframes <- function ( df_list, col='Count', method='wilcox') {
  
  combinations <- combn(names(df_list), 2)
  result <- data.frame(matrix(NA, nrow=dim(combinations)[2], ncol=3))
  colnames(result) <- c('Compares', 'Pvalue', 'Result')
  
  for ( i in 1:dim(combinations)[2] ) {
    
    df1 <- df_list[combinations[,i][1]][[1]]
    df2 <- df_list[combinations[,i][2]][[1]]
    # if not enough observations
    if ( length(df1[,col])>1 & length(df2[,col])>1 ) {
      if ( method == 'wilcox' ) {
        test <- wilcox.test(df1[,col], df2[,col])
      } else if (method == 'student' ) {
        test <- t.test(df1[,col], df2[,col])
      }
    } else {
      test <- list()
      test$p.value <- NA
    }

    # result 
    result[i,1] <- paste(combinations[,i][1], 'VS', combinations[,i][2])
    result[i,2] <- test$p.value
    if ( is.na(test$p.value) ) {
      result[i,3] <- 'Not enough observations'
    } else if ( test$p.value <= 0.05 ) {
      result[i,3] <- 'Significantly different'
    } else {
      result[i,3] <- 'Consistency'
    }
  }
  
  return(result)
  
}



# summarize stats of key (like PMean) values in single condition (like MutationType)
SummarizeKeyStatsByColumn <- function ( df, alldata_stats, key, col, cutoff=20 ) {
  
  # labels
  higher_label <- paste0(col, '>', cutoff)
  lower_label <- paste0(col, '<', cutoff)
  
  # column cutoff stats
  higherDf <- df[df[,col]>=cutoff,]
  lowerDf <- df[df[,col]<cutoff,]
  higher_cutoff_stats <- StatsByKey(higherDf, key=key, higher_label)
  lower_cutoff_stats <- StatsByKey(lowerDf, key=key, lower_label)
  
  # alldata comparing column cutoff
  alldata_compare_list <- list(alldata_stats, higher_cutoff_stats, lower_cutoff_stats)
  names(alldata_compare_list) <- c('All', higher_label, lower_label)
  all_test <- TestDataframes(alldata_compare_list, col='Count', method='wilcox')
  
  # key values comparing column cutoff
  key_types <- sort(unique(df[,key]))
  key_test <- list()
  for ( k in key_types ) {
    compare_list <- list(higherDf[higherDf[,key]==k,], lowerDf[lowerDf[,key]==k,])
    names(compare_list) <- c(higher_label, lower_label)
    tmp <- TestDataframes(compare_list, col='AF', method='student')
    tmp$Group <- k
    key_test[[k]] <- tmp
  }
  key_compare <- do.call(rbind, key_test)
  
  # return list
  results <- list(higher_cutoff_stats, lower_cutoff_stats, all_test, key_compare)
  names(results) <- c(higher_label, lower_label, 'WilcoxTestAlldata', 'WilcoxTestKeys')
  return(results)
  
}



# write to sheet
WriteToSheet <- function ( data, excelbook, sheetname ) {
  
  # check sheet if exist
  if ( ! sheetname %in% names(getSheets(excelbook)) ) {
    sheet <- createSheet(excelbook, sheetName=sheetname)
  } else {
    sheet <- getSheets(excelbook)[[sheetname]]
  }
  
  # write to sheet according to data type
  if ( class(data) == 'list' ) {
    writerow <- 1
    for ( df in names(data) ) {
      addDataFrame(data[[df]], sheet=sheet, col.names = TRUE, row.names = FALSE, startRow = writerow, startColumn = 1)
      writerow <- writerow + dim(data[[df]])[1] + 2
    }
    
  } else if ( class(data) == 'data.frame' ) {
    for ( element in names(data) ) {
      writerow <- 1
      addDataFrame(data, sheet=sheet, col.names = TRUE, row.names = FALSE, startRow = writerow, startColumn = 1)
    }
  }

}



# main 

Main <- function ( file, output='output', dir='./', mode='both' ) {
  
  # output files
  if ( is.null(dir) ) {
    dir <- dirname(file)
  } 
  dir.create(dir, recursive=T, showWarnings=F)
  if ( is.null(output) ) {
    output <- basename(file)
  } 
  output <- paste(dir, output, sep='/')
  
  # read input file
  df <- ReadVardict( file )
  
  # all data stats list
  alldata_stats_list <- list()
  alldata_strandbias_list <- list()
  excelbook <- createWorkbook()
  
  # summarize by each categories
  for ( k in keys ) {
    
    # all data stats
    mutation_AllData <- StatsByKey( df, key=k, 'AllData' )
    alldata_stats_list[[k]] <- mutation_AllData
    alldata_strandbias_list[[k]] <- StrandBiasSummary( df, key=k )

    # for each condition
    for ( c in 1:length(columns) ) {
      label <- paste0(k, 'By', columns[c], cutoffs[c])
      mutation_condition <- SummarizeKeyStatsByColumn( df, mutation_AllData, key=k, col=columns[c], cutoff=cutoffs[c] )
      
      # output format
      if ( mode=='both' | mode=='xlsx' ) {
        # to xlsx
        WriteToSheet(mutation_condition, excelbook = excelbook, sheetname = label)
      }
      if ( mode=='both' | mode=='tsv' ) {
        # to tsv
        write.table(mutation_condition[[1]], paste0(output,'.',label,'higher.tsv'), sep="\t", row.names=F, quote=F)
        write.table(mutation_condition[[2]], paste0(output,'.',label,'lower.tsv'), sep="\t", row.names=F, quote=F)
      }
      
    }
    
  }
  
  if ( mode=='both' | mode=='xlsx' ) {
    WriteToSheet(alldata_stats_list, excelbook = excelbook, sheetname = 'AlldataStats')
    WriteToSheet(alldata_strandbias_list, excelbook = excelbook, sheetname = 'AlldataStrandBias')
    saveWorkbook(excelbook, file=paste0(output,'.xlsx'))
  }
}





if (FALSE) {
  # comments
  excelbook <- createWorkbook()
  sheets <- getSheets(excelbook)
  sheet <- createSheet(excelbook, sheetName = 'test')
  
  addDataFrame(df, sheet=sheet, col.names = TRUE, row.names = FALSE, startRow = 1, startColumn = 1)
  # Add the plot created previously
  addPicture("boxplot.png", sheet, scale = 1, startRow = 4, startColumn = 1)
  # Remove the plot from the disk
  res<-file.remove("boxplot.png")
  saveWorkbook(excelbook, file=outfile)
#write.xlsx(alldata_vafstats$MutationDF, outfile, sheetName = "MutationStats", 
#           col.names = TRUE, row.names = FALSE, append = FALSE)
}



# define arguments
option_list = list(
  make_option("--input", type="character", default=NULL, help="vardict output", metavar="character"),
  make_option("--out", type="character", default=NULL, help="output file name (without path), default is input file name", metavar="character"),
  make_option("--dir", type="character", default=NULL, help="output directory name, default same as input file", metavar="character"),
  make_option("--mode", type="character", default='both', help="output file format, 'tsv' OR 'xlsx' OR 'both'", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

Main(opt$input, opt$out, opt$dir, opt$mode)
#Main('F:\\项目\\0930-突变模式\\test.alt_mut', output='test.alt_mut', dir='F:\\项目\\0930-突变模式\\test', mode='both')
