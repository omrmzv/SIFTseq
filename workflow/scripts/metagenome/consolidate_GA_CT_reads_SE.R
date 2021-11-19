library(data.table)

args= commandArgs(trailingOnly = TRUE)

GA_file = args[[1]] 
CT_file = args[[2]]

tblat1_file = args[[3]]

GA = fread(GA_file, col.names = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'strand', 'taxid'))
CT = fread(CT_file, col.names = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'strand', 'taxid'))

GA = GA[GA$length/GA$qlen >=0.9, ]
CT = CT[CT$length/CT$qlen >=0.9, ]

GA$taxid = as.character(GA$taxid)
CT$taxid = as.character(CT$taxid)

tblat1 = rbind(GA, CT)

set.seed(42)

rows <- sample(nrow(tblat1))

tblat1 <- data.frame(tblat1)[rows, ]

to_keep = !duplicated(tblat1$qseqid)

tblat1 = tblat1[to_keep, ]

tblat1$qseqid <- paste0(tblat1$qseqid, '-1')

fwrite(x = tblat1, file = tblat1_file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)