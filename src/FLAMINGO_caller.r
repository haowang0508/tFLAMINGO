args <- commandArgs(T)
input_pd <- read.table(args[1])
input_if <- read.table(args[2])
cell_id <- args[3]
output_file_path <- args[4]
backbone_res <- as.numeric(args[5])
final_res <- as.numeric(args[6])
N <- dim(input_pd)[1]
m <- backbone_res/final_res
bin_count <- 1+floor(N/m)
for(i in 1:bin_count){
    start <- (i-1)*m+1
    end <- min(i*m,N)
    tmp_bin_pd <- input_pd[start:end,start:end]
    tmp_bin_if <- input_if[start:end,start:end]
    pd_fn <- paste0(output_file_path,"/PD_Cell_",cell_id,"_bin_",i,".txt")
    if_fn <- paste0(output_file_path,"/IF_Cell_",cell_id,"_bin_",i,".txt")
    write.table(tmp_bin_pd,pd_fn,sep='\t',quote=F,col.names=F,row.names=F)
    write.table(tmp_bin_if,if_fn,sep='\t',quote=F,col.names=F,row.names=F)
}
