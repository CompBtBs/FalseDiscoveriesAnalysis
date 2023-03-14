###############################################
#
# convergence.R  Geweke and Raftery-lewis tests for convergence
#
# Version 1.0 February 2023
#
# Authors: 
#    - Bruno G. Galuzzi <bruno.galuzzi@unimib.it> (Department of Biotechnology and Biosciences, University of Milano-Bicocca)
#    - Luca Milazzo <l.milazzo1@campus.unimib.it> (Department of Informatics, Systems, and Communications, University of Milano-Bicocca)
#    - Chiara Damiani <chiara.damiani@unimib.it> (Department of Biotechnology and Biosciences, University of Milano-Bicocca)
#
###############################################

require(coda)
library(data.table)
library(hash)
library(rstudioapi)

# Function used to create a new folder by the path parameter
createFolder <- function(path){
  dir.create(path, showWarnings = TRUE)
}



# Computes the convergence coefficients of the sampled fluxes
# with the Geweke and Raftery-Lewis diagnostics
# Parameters
# - models --> models names
# - algorithms --> list of algorithms
# - thinnings --> list of thinnings
# - executionsPerSamples --> executions per samples
# - executionsPerSamplesCbs3 --> executions per samples for CBS3
# - tests --> list of diagnostic to compute
# - samplesFolder --> target folder for all the samples
# - resultFolder --> target folder for the results
# Notes:
# - It is not possibile to compute the raftery Lewis coefficients when the chain
# has a length smaller than 3746 samples (it will return Nan values).
# - The geweke.diag function returns Nan if there is perfect convergence and
# -Inf with a total absence of convergence.
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
samples[["cbs3"]] <- c(1000)

executionsPerSamples = 20
executionsPerSamplesCbs3 = 20

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

samplesFolder = "../../samples/"

resultFolder = "../../results/convergence/"

convergence(models, algorithms, thinnings, executionsPerSamples, executionsPerSamplesCbs3,  tests, samplesFolder, resultFolder)
