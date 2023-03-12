library(hash)
library(stringr)
library(rstudioapi)

createFolder <- function(path){
  dir.create(path, showWarnings = TRUE)
}


#Symmetric KLD

My_KLD <- function(P,Q){
  p <- P$y
  q <- Q$y
  if (all(q==0)==TRUE){eps <- min(p[which(p!=0)])*0.1}
  else{eps <-  min(min(p[which(p!=0)]),min(q[which(q!=0)]))*0.1}
  w_p <- p < eps
  if (any(w_p)) p[w_p] <- eps # to avoid numerical instabilities
  w_q <- q < eps  
  if (any(w_q)) q[w_q] <- eps # to avoid numerical instabilities
  f <- function(i){p[i]*log(p[i]/q[i])}
  dx <- P$x[2]-P$x[1]
  KLD_value <- dx*0.5*(f(1)+2*sum(f(2:(length(p)-1) ))+f(length(p)))    
  return(KLD_value)
}



kld <- function(models, algorithms, thinnings,  pairs, samplesFolder, resultFolder, chunkCut){
  
  createFolder(resultFolder)
  
  for (model in models){
    createFolder(paste(resultFolder , "/", model, sep = ""))
      for (algorithm in algorithms){
        if(algorithm == "chrr"){
          rowNames = "Row"
        }else{
          rowNames = "X"
        }
        for (thinning in thinnings[[algorithm]]){
          if(algorithm  == "cbs3"){
            createFolder(paste(resultFolder , model, "/" ,algorithm, "groupedBy", thinning , sep=""))
            createFolder(paste(resultFolder , model, "/" ,algorithm, "groupedBy", thinning , "/chunks/", sep=""))
          }else{
            createFolder(paste(resultFolder , model, "/" , algorithm, "Thinning", thinning , sep=""))
            createFolder(paste(resultFolder , model, "/" ,algorithm, "Thinning", thinning , "/chunks/", sep=""))
          }
          
          
          
          if(algorithm == "cbs3"){
            
            dfTest = read.csv(file = paste(samplesFolder ,model ,"/" , algorithm , "groupedBy" , thinning , "/" ,
                                           pairs[[algorithm]][[1]][1], ".csv" ,sep = ""),header=TRUE,  row.names=rowNames)
            
          }else{
            
          
            
              dfTest = read.csv(file = paste(samplesFolder ,model ,"/" , algorithm , "Thinning" , thinning , "/" ,
                                             pairs[[algorithm]][[1]][1], ".csv" ,sep = ""),header=TRUE, row.names=rowNames)
           
            
            
        
          }
          
         
         
          
          nReactions = ncol(dfTest)
          chunk = 0
          result_df = data.frame(matrix(ncol = nReactions + 1 , nrow = 0))
          colnames(result_df) = append(list("pair") , colnames(dfTest))
          
          
          for (pair in pairs[[algorithm]]){
            
            if(algorithm == "cbs3"){
              
              
              file1 = paste(samplesFolder ,model ,"/" , algorithm , "groupedBy" , thinning , "/" ,
                            pair[[1]] , ".csv", sep="")
              file2 = paste(samplesFolder ,model ,"/" , algorithm , "groupedBy" , thinning , "/" ,
                            pair[[2]] , ".csv", sep="")
            }else{
              

              
              file1 = paste(samplesFolder ,model ,"/" , algorithm , "Thinning" , thinning , "/" ,
                            pair[[1]] , ".csv", sep="")
              file2 = paste(samplesFolder ,model ,"/" , algorithm , "Thinning" , thinning , "/" ,
                            pair[[2]] , ".csv", sep="")
            }
            
            
            df1 = read.csv(file = file1, header=TRUE, row.names=rowNames)
            df2 = read.csv(file = file2, header=TRUE,  row.names=rowNames)
            
            resList = list(paste(pair[[1]] ,pair[[2]], sep = "_"))
            
            for(h in seq(1, ncol(df1), by = 1)){
              X = df1[, h]
              Y = df2[, h]
              res1 = My_KLD(density(X, n=5000,from =min(X),to=max(X),cut=10),
                            density(Y, n=5000,from =min(X) ,to=max(X), cut=10 ))
              res2 = My_KLD(density(Y, n=5000,from =min(Y),to=max(Y),cut=10),
                            density(X, n=5000,from =min(Y) ,to=max(Y), cut=10 ))
              resList = append(resList, (res1+res2)/2)
            }
            
            result_df[nrow(result_df) + 1,] <- resList
            
            if(nrow(result_df)==chunkCut){
              if(algorithm  == "cbs3"){
                write.csv(result_df, paste( resultFolder , model, "/" ,algorithm, "groupedBy", thinning , "/chunks/" , "chunk_", 
                                            chunk, ".csv", sep = ""), row.names = TRUE)
              }else{
                write.csv(result_df, paste( resultFolder , model, "/" ,algorithm, "Thinning", thinning , "/chunks/" , "chunk_", 
                                            chunk, ".csv", sep = ""), row.names = TRUE)
              }
             
              result_df = data.frame(matrix(ncol = nReactions + 1, nrow = 0))
              colnames(result_df) = append(list("pair") ,colnames(dfTest))
              chunk = chunk + 1
            }
            
          }
          
          if(algorithm  == "cbs3"){
            files = list.files(path=paste(resultFolder , model, "/" ,algorithm, "groupedBy", thinning , "/chunks/", sep=""), 
                               all.files = TRUE, 
                               full.names = TRUE, 
                               pattern="chunk_")
          }else{
            files = list.files(path=paste(resultFolder , model, "/" ,algorithm, "Thinning", thinning , "/chunks/", sep=""), 
                               all.files = TRUE, 
                               full.names = TRUE, 
                               pattern="chunk_")
          }
          
          
          
          files = str_sort(files, numeric = TRUE)
          isFirst = TRUE
          
          for(file in files){
            if(isFirst){
              df_result = read.csv(file = file, header=TRUE, row.names="X")
              isFirst = FALSE
            }else{
              df_result = rbind(df_result,read.csv(file = file, header=TRUE, row.names="X"))
            }
          }
          
          if(algorithm  == "cbs3"){
            write.csv(df_result, paste(resultFolder , model, "/" ,algorithm, "groupedBy", thinning , "/kld.csv", sep = ""), row.names = TRUE)
          }else{
            write.csv(df_result, paste(resultFolder , model, "/" ,algorithm, "Thinning", thinning , "/kld.csv", sep = ""), row.names = TRUE)
          }
          
    
        }
      }
   
    
  }
  
}



models <- c('ENGRO 1', 'ENGRO 2')

algorithms <- c('achr', 'optgp', 'chrr', 'cbs3')


thinnings <- hash()
thinnings[["achr"]] <- list(1, 2, 3)
thinnings[["optgp"]] <- list(1, 2, 3)
thinnings[["chrr"]] <- list(1, 2, 3)
thinnings[["cbs3"]] <- list(1000)


samples <- hash()
samples[["achr"]] <- seq(1000, 4000, 1000)
samples[["optgp"]] <- seq(1000, 4000, 1000)
samples[["chrr"]] <- seq(1000, 4000, 1000)
samples[["cbs3"]] <- c(1000)

executionsPerSamples = 3
executionsPerSamplesCbs3 = 3

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

samplesFolder = "../../samples/"

resultFolder = "../../results/KLD/"

chunkCut = 3

pairs <- hash()

for (algorithm in algorithms){
  pairlist = list()
  for(nsample in samples[[algorithm]]){
    if(algorithm == "cbs3"){
      executionsPerSamplesUse = executionsPerSamplesCbs3
    }else{
      executionsPerSamplesUse = executionsPerSamples
    }
    for(nexec1 in seq(0, executionsPerSamplesUse-2, by=1)){
      for(nexec2 in seq(nexec1 + 1, executionsPerSamplesUse-1, by = 1)){
        if(nexec1 != nexec2){
          if(algorithm != "cbs3"){
            el1 = paste(nsample, nexec1, algorithm ,sep = "_") #1000_0_achr
            el2 = paste(nsample, nexec2, algorithm ,sep = "_")
            pair = list(el1, el2)
            pairlist = append(pairlist, list(pair))
          }else{
            el1 = paste(nexec1, 0, algorithm ,sep = "_")
            el2 = paste(nexec2, 0, algorithm ,sep = "_")
            pair = list(el1, el2)
            pairlist = append(pairlist, list(pair))
          }
          
        }
      }
    }
  }
  pairs[[algorithm]] = pairlist
}


kld(models, algorithms, thinnings,  pairs, samplesFolder, resultFolder , chunkCut)
