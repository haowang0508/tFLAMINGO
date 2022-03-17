# just the draft, need to add parameters for chromosome id and assembly
args <- commandArgs(T)
library(GenomicFeatures)
library(Matrix)


chr_name <- args[1]
resolution <- as.numeric(args[2])
input_path <- args[3]
opt_path <- args[4]
assembly <- args[5]
chr_size <- getChromInfoFromUCSC(assembly)
load("./data/Linear_model_for_transformation.Rdata")
setwd(paste0(input_path,chr_name))
file <- dir()
#file <- file[grep('RawCount',file)]
order_id <- as.numeric(sapply(strsplit(file,'_'),function(x){x[2]}))
file <- file[order(order_id,decreasing=F)]
sc_mat <- list()
for(i in 1:length(file)){
  sc_mat[[i]] <- as.matrix(as.data.frame(data.table::fread(file[i])))
}

# normalize the val towards the averaged count
averaged_count <- mean(sapply(sc_mat,sum))

for(i in 1:length(sc_mat)){
  sc_mat[[i]] <- sc_mat[[i]]/sum(sc_mat[[i]])*averaged_count
}


# missing rate
missing_rate <- c()
for(i in sc_mat){
  missing_rate <- c(missing_rate,1-mean(i>0))
}

# given a matrix, get the values, row index and column index
sc_frame <- list()
for(i in 1:length(sc_mat)){
  tmp_frame <-  summary(Matrix(sc_mat[[i]],sparse=T))
  id_x <- tmp_frame[,1]
  id_y <- tmp_frame[,2]
  sc_frame[[i]] <- data.frame(x=tmp_frame[,1],y=tmp_frame[,2],dist = abs(tmp_frame[,2]-tmp_frame[,1])*resolution,IF=tmp_frame[,3],
                              missing_rate=missing_rate[i],if_x = diag(as.matrix(sc_mat[[i]]))[id_x],if_y=diag(as.matrix(sc_mat[[i]]))[id_y])
}

pred_val <- function(x,lm_list){
  df_list <- split(x,floor(exp(x[,4])/5e6)+1)
  res <- list()
  for(i in names(df_list)){
    i <- as.numeric(i)
    lm_model <- lm_list[[ min(i,length(lm_list)) ]]
    tmp_res <- predict(lm_model,df_list[[as.character(i)]][,-c(1,2)])
    res[[i]] <- data.frame(df_list[[as.character(i)]][,c(1,2)],tmp_res)
  }
  return(res)
}

x <- sc_frame[[1]]
x <- x[,c(1,2,4,3,5)]
x <- as.matrix(x)
x <- subset(x,x[,4]>0)
x[,3] <- log(x[,3])
x[,4] <- log(x[,4])
x <- as.data.frame(x)
x$'ifxy' <- diag(as.matrix(sc_mat[[1]] ))[x[,1]]*diag(as.matrix(sc_mat[[1]] ))[x[,2]]
transfered_if <- pred_val(x,lm_list)
transfered_if <- do.call(rbind,transfered_if)
N <- floor(chr_size[match(chr_name,chr_size[,1]),2]/resolution)+1
tf_mat <- sparseMatrix(
  i = transfered_if[,1], 
  j = transfered_if[,2], 
  x = transfered_if[,3],
  dims = c(N, N), 
)
#tf_mat <- as.matrix(tf_mat)
y <- exp(tf_mat)-1

y <- y^(-0.25)
rm_id <- which(apply(y,1,min) == Inf)
y[which(y==Inf)] <- NA


pred_sc_if <- list()
for(i in 1:length(sc_frame)){
x <- sc_frame[[i]]
x <- x[,c(1,2,4,3,5)]
x <- as.matrix(x)
x <- subset(x,x[,4]>0)
x[,3] <- log(x[,3])
x[,4] <- log(x[,4])
x <- as.data.frame(x)
x$'ifxy' <- diag(as.matrix(sc_mat[[1]] ))[x[,1]]*diag(as.matrix(sc_mat[[1]] ))[x[,2]]
transfered_if <- pred_val(x,lm_list)
transfered_if <- do.call(rbind,transfered_if)
N <- floor(chr_size[match(chr_name,chr_size[,1]),2]/resolution)+1
tf_mat <- sparseMatrix(
  i = transfered_if[,1], 
  j = transfered_if[,2], 
  x = transfered_if[,3],
  dims = c(N, N), 
)
pred_sc_if[[i]] <- as.matrix(exp(tf_mat)-1)
}

for(i in 1:length(pred_sc_if) ){
  data.table::fwrite(pred_sc_if[[i]],paste0(opt_path,'/RawCount_Cell_',i,".txt"),col.names=F,row.names=F,sep='\t',quote=F,na='NA')
}


