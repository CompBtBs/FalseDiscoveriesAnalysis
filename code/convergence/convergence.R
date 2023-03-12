require(coda)
library(data.table)
library(hash)
library(rstudioapi)

createFolder <- function(path){
  dir.create(path, showWarnings = TRUE)
}

#3746 error for raftery Lewis--> if chain length not > 3746 --> rl = Nan
#geweke.diag returns Nan if there is perfect convergence and -Inf with total absence of convergence

convergence <- function(models, algorithms, thinnings, executionsPerSamples, executionsPerSamplesCbs3,  tests, samplesFolder, resultFolder){

  for (model in models){
    createFolder(paste(resultFolder , model, sep = ""))
    for(test in tests){
      if(test == "geweke" || test == "raftery-lewis"){
        createFolder(paste(resultFolder,   sep=""))
        createFolder(paste(resultFolder, model,    sep=""))
        createFolder(paste(resultFolder , "/", model, "/" , test, sep=""))
      }else{
        return ("ERROR - test not supported")
      }
      for (algorithm in algorithms){
        for (thinning in thinnings[[algorithm]]){
          
          if(algorithm  == "cbs3"){
            createFolder(paste(resultFolder , model , "/" , test, "/", algorithm, "groupedBy", thinning ,sep=""))
          }else{
            createFolder(paste(resultFolder , model , "/" , test, "/", algorithm, "Thinning", thinning ,sep=""))
          }
          
          for (nsample in samples[[algorithm]]){
            if(algorithm == "cbs3"){
              executionsPerSamplesToUse = executionsPerSamplesCbs3 - 1
            }else{
              executionsPerSamplesToUse = executionsPerSamples - 1
            }
            for(nexec in seq(0, executionsPerSamplesToUse, by = 1)){
              
              if(algorithm  == "cbs3"){
                file = paste(samplesFolder , model ,"/" , algorithm , "groupedBy" , thinning , "/" ,
                             nexec , "_" , 0 , "_" , algorithm, ".csv", sep="")
              }else{
                file = paste(samplesFolder ,model ,"/" , algorithm , "Thinning" , thinning , "/" ,
                             nsample , "_" , nexec , "_" , algorithm, ".csv", sep="")
              }
              df=fread(file, sep=",",header = TRUE)
              if(algorithm == "chrr"){
                rownames(df) = df$Row
                df$Row=NULL
              }else{
                
                rownames(df) = df$V1
                df$V1=NULL
              }
              
              df_matrix=as.matrix(df)
              result_list_geweke=c()
              result_list_rl=c()
              nreactions=dim(df_matrix)[2]
              for (i in 1:nreactions){
                b=as.mcmc(df_matrix[,i])
                if(test == "geweke"){
                  result_list_geweke[i]=as.numeric(geweke.diag(b[1:nsample])$z)
                }else{
                  result_list_rl[i]= as.numeric(raftery.diag(b[1:nsample])$resmatrix[4])
                }
                
              }
              if(test == "geweke"){
                
                data_final_geweke <- data.frame(
                  reactions = colnames(df),
                  result = result_list_geweke
                )
                
                if(algorithm  == "cbs3"){
                  write.csv(data_final_geweke,paste(resultFolder , model , "/" , test, "/", algorithm, "groupedBy", thinning , 
                                                    "/" , nexec , "_" , 0 , '.csv', sep=""))
                }else{
                  write.csv(data_final_geweke,paste(resultFolder , model , "/" , test, "/", algorithm, "Thinning", thinning  ,
                                                    "/" , nsample , "_" , nexec , '.csv', sep=""))
                }
                
              }else{
                
          
                
                data_final_rl <- data.frame(
                  reactions = colnames(df),
                  result = result_list_rl
                )
                
                if(algorithm  == "cbs3"){
                  write.csv(data_final_rl,paste(resultFolder , model , "/" , test, "/", algorithm, "groupedBy", thinning , 
                                                    "/" , nexec , "_" , 0 , '.csv', sep=""))
                }else{
                  write.csv(data_final_rl,paste(resultFolder , model , "/" , test, "/", algorithm, "Thinning", thinning  ,
                                                    "/" , nsample , "_" , nexec , '.csv', sep=""))
                }
                
                
              }
             
             
            }
          }
        }
      }
    }
    
  }
  
}

##########
models <- c('ENGRO 1', 'ENGRO 2')

algorithms <- c('achr', 'optgp', 'chrr', 'cbs3')

thinnings <- hash()
thinnings[["achr"]] <- list(1, 2, 3)
thinnings[["optgp"]] <- list(1, 2, 3)
thinnings[["chrr"]] <- list(1, 2, 3)
thinnings[["cbs3"]] <- list(1000)

tests <- c('geweke')

samples <- hash()
samples[["achr"]] <- seq(1000, 4000, 1000)
samples[["optgp"]] <- seq(1000, 4000, 1000)
samples[["chrr"]] <- seq(1000, 4000, 1000)
samples[["cbs3"]] <- c(1000)

executionsPerSamples = 3
executionsPerSamplesCbs3 = 3

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

samplesFolder = "../../samples/"

resultFolder = "../../results/convergence/"

convergence(models, algorithms, thinnings, executionsPerSamples, executionsPerSamplesCbs3,  tests, samplesFolder, resultFolder)

########


models <- c('ENGRO 1', 'ENGRO 2')

algorithms <- c('achr', 'optgp', 'chrr', 'cbs3')

thinnings <- hash()
thinnings[["achr"]] <- list(1, 10, 100)
thinnings[["optgp"]] <- list(1, 10, 100)
thinnings[["chrr"]] <- list(1, 10, 100)
thinnings[["cbs3"]] <- list(1000)

tests <- c('geweke', 'raftery-lewis')

samples <- hash()
samples[["achr"]] <- seq(1000, 30000, 1000)
samples[["optgp"]] <- seq(1000, 30000, 1000)
samples[["chrr"]] <- seq(1000, 30000, 1000)
samples[["cbs3"]] <- list(-1)

executionsPerSamples = 20
executionsPerSamplesCbs3 = 20

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

samplesFolder = "../../samples/"

resultFolder = "../../results/convergence/"

convergence(models, algorithms, thinnings, executionsPerSamples, tests, samplesFolder, resultFolder)



df=fread("C:/Users/LM856702/Documents/FalseDiscoveriesAnalysis/samples/ENGRO 1/achrThinning1/1000_0_achr.csv" , sep=",",header = TRUE)
df_chrr=fread("C:/Users/LM856702/Documents/FalseDiscoveriesAnalysis/samples/ENGRO 1/chrrThinning1/1000_0_chrr.csv" , sep=",",header = TRUE)


