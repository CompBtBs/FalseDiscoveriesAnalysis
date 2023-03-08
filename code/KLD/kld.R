library(hash)
library(stringr)

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


pairs <- hash() 

for (algo in c('achr', 'optgp', 'chrr')){
  pairlist = list()
  #for(nsample in seq(1000, 30000, by=1000)){
  for(nsample in c(1000, 5000, 10000, 30000)){
    for(nexec1 in seq(0, 18, by=1)){
      for(nexec2 in seq(nexec1 + 1, 19, by = 1)){
        if(nexec1 != nexec2){
          el1 = paste(nsample, nexec1, algo ,sep = "_") #1000_0_achr
          el2 = paste(nsample, nexec2, algo ,sep = "_")
          pair = list(el1, el2)
          pairlist = append(pairlist, list(pair))
        }
      }
    }
  }
  pairs[[algo]] = pairlist
}


for (model in c('ENGRO 2')){
  for (algo in c( 'chrr')){
    for (thinning in c(1, 10, 100)){
      print(paste(model, algo, thinning, sep = " "))
      chunk = 0
      result_df = data.frame(matrix(ncol = 496, nrow = 0))
      colnames(result_df) = append(list("pair") ,colnames(read.csv(file = paste("D:/flux_sampling_analysis/sampling/" , model , "/" , toupper(algo), 
                                                                                "Thinning"  , thinning, "/", "1000_0_",algo ,".csv" ,sep = ""), header=TRUE, row.names="X")))
      for(pair in pairs[[algo]]){
        file1 = paste("D:/flux_sampling_analysis/sampling/" , model , "/" , toupper(algo), 
                      "Thinning" , thinning, "/", pair[[1]], ".csv" ,sep = "")
        file2 = paste("D:/flux_sampling_analysis/sampling/" , model , "/" , toupper(algo), 
                      "Thinning"  , thinning, "/", pair[[2]], ".csv" ,sep = "")
        
        df1 = read.csv(file = file1, header=TRUE, row.names="X")
        df2 = read.csv(file = file2, header=TRUE, row.names="X")
        
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
        if(nrow(result_df)==50){
          write.csv(result_df, paste("D:/flux_sampling_analysis/result/", model, "/kld/symmetricChunks/", toupper(algo),"Thinning",thinning,"chunk_", chunk,
                                     ".csv", sep = ""), row.names = TRUE)
          result_df = data.frame(matrix(ncol = 496, nrow = 0))
          colnames(result_df) = append(list("pair") ,colnames(read.csv(file = paste("D:/flux_sampling_analysis/sampling/" , model , "/" , toupper(algo), 
                                                                                    "Thinning"  , thinning, "/", "1000_0_",algo ,".csv" ,sep = ""), header=TRUE, row.names="X")))
          chunk = chunk + 1
        }
        
      }
      write.csv(result_df, paste("D:/flux_sampling_analysis/result/", model, "/kld/symmetricChunks/", toupper(algo),"Thinning",thinning,"chunk_", chunk,
                                 ".csv", sep = ""), row.names = TRUE)
    }
  }
}

#Merge chunks
for (model in c('ENGRO 2')){
  for (algo in c( 'chrr')){
    for (thinning in c(1, 10, 100)){
      files = list.files(path=paste("D:/flux_sampling_analysis/result/", model, "/kld/symmetricChunks/", sep=""), 
                   all.files = TRUE, 
                   full.names = TRUE, 
                   pattern=paste(toupper(algo), "Thinning", thinning,"chunk", sep=""))
      
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
      print(nrow(df_result))
      write.csv(df_result, paste("D:/flux_sampling_analysis/result/", model, "/kld/", toupper(algo),"Thinning",thinning,
                                 ".csv", sep = ""), row.names = TRUE)
      
    }
  }
}





