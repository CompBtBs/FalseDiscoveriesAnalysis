require(coda)
library(data.table)

nsample = 1000

# geweke.diag returns Nan if ther's perfect convergence and -Inf with total absence of convergence

for (model in c('ENGRO 1', 'ENGRO 2')){
     for (nexec1 in seq(0, 19, by=1)){
        print(paste("**********", model , nexec1, sep = " "))
        nexec2 = 0
        file=paste("D:/flux_sampling_analysis/sampling/" ,model ,"/randomOF/"  ,
                   nexec1 , "_" , nexec2 ,  ".csv", sep="")
        df=fread(file, sep=",",header = TRUE)
        rownames(df) = df$V1
        df$V1=NULL
        df_matrix=as.matrix(df)
        result_list_geweke=c()
        result_list_rl=c()
        nreactions=dim(df_matrix)[2]
        for (i in 1:nreactions){
            b=as.mcmc(round(df_matrix[,i], 6))
            result_list_geweke[i]=as.numeric(geweke.diag(b[1:nsample])$z)
        }
        data_final_geweke <- data.frame(
            reactions = colnames(df),
            result = result_list_geweke
          )
        write.csv(data_final_geweke,paste("D:/flux_sampling_analysis/result/", model, "/convergence/geweke/randomOF_"  , nexec1 , "_" , nexec2 , '.csv', sep=""))
        
    }
  }
      
    



