require(coda)
library(data.table)

# geweke.diag returns Nan if ther's perfect convergence and -Inf with total absence of convergence

for (model in c('ENGRO 2')){
    for (algo in c('chrr')){
        for (thinning in c(1,10,100)){
            for (nsample in seq(1000, 30000, by=1000)){
                print(paste("**********", model , toupper(algo), "Thinning" , thinning   , nsample, sep = " "))
                for(nexec in seq(0, 19, by = 1)){
                    file=paste("D:/flux_sampling_analysis/sampling/" ,model ,"/" , toupper(algo) , "Thinning" , thinning , "/" ,
                                nsample , "_" , nexec , "_" , algo , ".csv", sep="")
                    df=fread(file, sep=",",header = TRUE)
                    rownames(df) = df$V1
                    df$V1=NULL
                    df_matrix=as.matrix(df)
                    result_list_geweke=c()
                    result_list_rl=c()
                    nreactions=dim(df_matrix)[2]
                    for (i in 1:nreactions){
                        b=as.mcmc(df_matrix[,i])
                        result_list_geweke[i]=as.numeric(geweke.diag(b[1:nsample])$z)
                        result_list_rl[i]= as.numeric(raftery.diag(b[1:nsample])$resmatrix[4])
                        #3746 error --> if chain length not > 3746 --> rl = Nan
                    }
                    data_final_geweke <- data.frame(
                        reactions = colnames(df),
                        result = result_list_geweke
                      )
                    write.csv(data_final_geweke,paste("D:/flux_sampling_analysis/result/", model, "/convergence/geweke/"  , toupper(algo) , "Thinning" , thinning ,
                                             "_" , nsample , "_" , nexec , '.csv', sep=""))
                    data_final_rl <- data.frame(
                        reactions = colnames(df),
                        result = result_list_rl
                      )
                      write.csv(data_final_rl,paste("D:/flux_sampling_analysis/result/", model, "/convergence/rafteryLewis/"  , toupper(algo) , "Thinning" , thinning ,
                                                    "_" , nsample , "_" , nexec , '.csv', sep=""))
                }
        }
      }
    }
}


