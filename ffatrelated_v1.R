####### ####### ####### ####### FFAT-RELATED ####### ####### ####### ####### ####### ####### 
## RUN AS

#' R CMD BATCH --vanilla --slave --args --i=cert.txt --o=o.txt --all ffatrelated.R
#' --o outputfile
#' --i input file fasta format .txt is also fine
#' --all add if you want to screen each residue - in case you have multiple motifs

## Cite Cabukusta et al. 2020 Cell Reports


checkOptions <- function( args )
{
  i = grep("--o",args); 
  outputFile = ifelse( !any(i), NA, strsplit(args[i],"=")[[1]][2])
  i = grep("--i",args); 
  inputFile = ifelse( !any(i), NA, strsplit(args[i],"=")[[1]][2])
  i = grepl("--all",args); 
  all.or.not = as.logical(sum(i))
  bb = list(outputFile=outputFile,inputFile=inputFile, all.or.not=all.or.not);
  return (bb);
}
options(echo=F);
opts = checkOptions(commandArgs())

library(dnar)
x <- read.fa(opts$inputFile)
allornot <- opts$all.or.not
print(opts)

#SCORES OF PWM
scores <- 
  data.frame(
    Ideal = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","X"),
    p1=c(1,  1,  0, .5,  1,  1,  1,  1, 1.5, 1,  1, .5,  1,  1,  2,  0, .5,  1,  1,  1,  2),
    p2=c(4,  4,  4,  4,  0,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  2,  0,  4), 
    p3=c(1,  0,  1,  1,  0,  1,  1,  1,  1, .5,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  2), 
    p4=c(.5, 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  0,  2,  2,  0, .5,  2,  2,  2,  2), 
    p5=c(.5, 1,  4,  4,  2,  1,  2,  2, 2.5, 2,  2,  0,  2,  2,  3,  0 , 0  ,2,  2,  2,  2),
    p6=c(.5,.5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5,  0, .5, .5, .5, .5, .5, .5, .5,  2),
    p7=c(0,  1,  4,  4,  1,  0,  1,  1, 1.5, 1,  1, .5,  1,  1, 1.5, 0, .5,  1,  1,  1,  2) 
  )


#SCORING FUNCTION
ffatty5.related <- function(x, all.scores=F) {
  x.scores <- data.frame(name = NA, score = NA, pos = NA, seq = NA)
  for(i.q in 1:nrow(x)) {
    temp.seq <- x[i.q,2]
    temp.seq <- unlist(strsplit(temp.seq, ""))
    temp.scores <- data.frame(name = NA, score = NA, pos = NA, seq = NA)
    
    
    #N-term append
    temp.seq = append(c("X","X","X","X","X","X","X","X"), values = temp.seq)
    temp.seq = append(temp.seq, values = c("X","X","X","X","X","X","X","X"))
    
    #without N terminal acidic tract
    all.bits <- c()
    for (ii in 1:(length(temp.seq)-14)) {
      query.N <-temp.seq[(ii):(ii+6)]
      query.core <- temp.seq[(ii+7):(ii+13)]
      bit.core <- 
        sum(
          scores[grep(query.core[1],scores$Ideal),'p1'],
          scores[grep(query.core[2],scores$Ideal),'p2'],
          scores[grep(query.core[3],scores$Ideal),'p3'],
          scores[grep(query.core[4],scores$Ideal),'p4'],
          scores[grep(query.core[5],scores$Ideal),'p5'],
          scores[grep(query.core[6],scores$Ideal),'p6'],
          scores[grep(query.core[7],scores$Ideal),'p7']
        )
      all.bits <- c(all.bits, bit.core)
      seq <- c()
      for(i.allbits in 1:length(all.bits)) {
        seq <- c(seq, paste0(temp.seq[(i.allbits):(i.allbits+16)], collapse = "") )
      }
    }
    if (all.scores == F) {
      pos <-  match(min(all.bits, na.rm = T), all.bits)
      temp.scores$name <- x[i.q,1]
      temp.scores$score <- all.bits[pos]
      temp.scores$pos <- pos
      temp.scores$seq <- paste0(temp.seq[(pos):(pos+16)], collapse = "")
      x.scores <- rbind(x.scores, temp.scores)
    } else {
      temp.scores = data.frame(name = x[i.q,1],
                               pos = 1:length(all.bits),
                               score = all.bits,
                               seq = seq)
      x.scores <- rbind(x.scores, temp.scores)
    }
    
  }
  return(x.scores)
}

  

## TEST
#ffatty5.related(dnar::read.fa('EMD.txt'), all.scores = T)
#ffatty5.related(dnar::read.fa('CERT.txt'), all.scores = F)
#ffatty5.related(dnar::read.fa('test.txt'), all.scores = F)
#ffatty5.related(dnar::read.fa('cert+emd.txt'), all.scores = F)
#test <- ffatty5.related(dnar::read.fa('cert+emd.txt'), all.scores = T)

scores <- ffatty5.related(x, all.scores = allornot)
scores <- scores[complete.cases(scores),]

#output result
if(is.na(opts$outputFile))
{
  write.table(scores, 'output.txt',quote = F, sep = '\t',row.names = F)
}else
{
  write.table(scores, opts$outputFile,quote = F, sep = '\t',row.names = F)
  
}


#' Cabukusta 2020
#' CELL REPORTS 2020