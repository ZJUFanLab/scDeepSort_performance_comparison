# SCINA function
SCINA=function(exp,signatures,max_iter=100,convergence_n=10,convergence_rate=0.99,
               sensitivity_cutoff=1,rm_overlap=1,allow_unknown=1,log_file='SCINA.log'){
  check.inputs=function(exp, signatures, max_iter, convergence_n, convergence_rate, sensitivity_cutoff, rm_overlap, log_file){
    # Initialize parameters.
    quality=1
    def_max_iter=1000
    def_conv_n=10
    def_conv_rate=0.99
    def_dummycut=0.33
    allgenes=row.names(exp)
    # Check sequence matrices.
    if (any(is.na(exp))){
      cat('NA exists in expression matrix.',file=log_file,append=T)
      cat('\n',file=log_file,append=T)
      quality=0
    }
    # Check signatures.
    if (any(is.na(signatures))){
      cat('Null cell type signature genes.',file=log_file,append=T)
      cat('\n',file=log_file,append=T)
      quality=0
    }else{
      signatures=sapply(signatures,function(x) unique(x[(!is.na(x)) & (x %in% allgenes)]),simplify=F)
      # Remove duplicate genes.
      if(rm_overlap==1){
        tmp=table(unlist(signatures))
        signatures=sapply(signatures,function(x) x[x %in% names(tmp[tmp==1])],simplify=F)
      }
      # Check if any genes have all 0 counts
      signatures=sapply(signatures,function(x) x[apply(exp[x,,drop=F],1,sd)>0],simplify=F)
    }
    # Clean other parameters.
    if (is.na(convergence_n)){
      cat('Using convergence_n=default',file=log_file,append=T)
      cat('\n',file=log_file,append=T)
      convergence_n=def_conv_n
    }
    if (is.na(max_iter)){
      cat('Using max_iter=default',file=log_file,append=T)
      cat('\n',file=log_file,append=T)
      max_iter=def_max_iter
    }else{
      if (max_iter<convergence_n){
        cat('Using max_iter=default due to smaller than convergence_n.',file=log_file,append=T)
        cat('\n',file=log_file,append=T)
        max_iter=convergence_n 
      }
    }
    if (is.na(convergence_rate)){
      cat('Using convergence_rate=default.',file=log_file,append=T)
      cat('\n',file=log_file,append=T)
      convergence_rate=def_conv_rate
    }
    if (is.na(sensitivity_cutoff)){
      cat('Using sensitivity_cutoff=default.',file=log_file,append=T)
      cat('\n',file=log_file,append=T)
      sensitivity_cutoff=def_dummycut
    }
    # Return cleaned parameters.
    return(list(qual=quality,sig=signatures,
                para=c(max_iter,convergence_n,convergence_rate,sensitivity_cutoff)))
  }
  density_ratio=function(e,mu1,mu2,inverse_sigma1,inverse_sigma2){
    tmp1=colSums((e-mu1)*(inverse_sigma1%*%(e-mu1)))
    tmp2=colSums((e-mu2)*(inverse_sigma2%*%(e-mu2)))
    tmp=exp(-1/2*(tmp1+log(1/det(inverse_sigma1))-tmp2-log(1/det(inverse_sigma2))))
    tmp[tmp>1e200]=1e200
    tmp[tmp<1e-200]=1e-200
    return(tmp)
  }
  cat('Start running SCINA.',file=log_file,append=F)
  cat('\n',file=log_file,append=T)
  #Create a status file for the webserver
  status_file=paste(log_file,'status',sep='.')
  all_sig=unique(unlist(signatures))
  # Create low-expression signatures.
  invert_sigs=grep('^low_',all_sig,value=T)
  if(!identical(invert_sigs, character(0))){
    cat('Converting expression matrix for low_genes.',file=log_file,append=T)
    cat('\n',file=log_file,append=T)
    invert_sigs_2add=unlist(lapply(invert_sigs,function(x) strsplit(x,'_')[[1]][2]))
    invert_sigs=invert_sigs[invert_sigs_2add%in%row.names(exp)]
    invert_sigs_2add=invert_sigs_2add[invert_sigs_2add%in%row.names(exp)]
    sub_exp=-exp[invert_sigs_2add,,drop=F]
    row.names(sub_exp)=invert_sigs
    exp=rbind(exp,sub_exp)
    rm(sub_exp,all_sig,invert_sigs,invert_sigs_2add)
  }
  # Check input parameters.
  quality=check.inputs(exp,signatures,max_iter,convergence_n,convergence_rate,sensitivity_cutoff,rm_overlap,log_file)
  if(quality$qual==0){
    cat('EXITING due to invalid parameters.',file=log_file,append=T)
    cat('\n',file=log_file,append=T)
    cat('0',file=status_file,append=F)
    stop('SCINA stopped.')
  }
  signatures=quality$sig
  max_iter=quality$para[1]
  convergence_n=quality$para[2]
  convergence_rate=quality$para[3]
  sensitivity_cutoff=quality$para[4]
  # Initialize variables.
  exp=as.matrix(exp)
  exp=exp[unlist(signatures),,drop=F]
  labels=matrix(0,ncol=convergence_n, nrow=dim(exp)[2])
  unsatisfied=1
  if(allow_unknown==1){
    tao=rep(1/(length(signatures)+1),length(signatures))
  }else{tao=rep(1/(length(signatures)),length(signatures))}
  theta=list()
  for(i in 1:length(signatures)){
    theta[[i]]=list()
    theta[[i]]$mean=t(apply(exp[signatures[[i]],,drop=F],1,function(x) quantile(x,c(0.7,0.3))))
    tmp=apply(exp[signatures[[i]],,drop=F],1,var)
    theta[[i]]$sigma1=diag(tmp,ncol = length(tmp))
    theta[[i]]$sigma2=theta[[i]]$sigma1
  }
  theta1<- NULL
  for(marker_set in 1:length(theta)){
    if(rlang::is_empty(theta[[marker_set]]$sigma1) == TRUE){
      theta1 <- c(theta1,marker_set)
    }
  }
  if (!is.null(theta1)) {
    theta<- theta[-theta1]
    signatures<- signatures[-theta1]
    tao<- tao[-theta1]
  }
  sigma_min=min(sapply(theta,function(x) min(c(diag(x$sigma1),diag(x$sigma2)))))/100
  remove_times=0
  # Run SCINA algorithm.
  while(unsatisfied==1){
    prob_mat=matrix(tao,ncol=dim(exp)[2],nrow=length(tao))
    row.names(prob_mat)=names(signatures)
    iter=0
    labels_i=1
    remove_times=remove_times+1
    while(iter<max_iter){
      iter=iter+1
      # E step: estimate variables.
      for(i in 1:length(signatures)){
        theta[[i]]$inverse_sigma1=theta[[i]]$inverse_sigma2=chol2inv(chol(theta[[i]]$sigma1))
      }
      for (r in 1:dim(prob_mat)[1]){
        prob_mat[r,]=tao[r]*density_ratio(exp[signatures[[r]],,drop=F],theta[[r]]$mean[,1],
                                          theta[[r]]$mean[,2],theta[[r]]$inverse_sigma1,theta[[r]]$inverse_sigma2)
      }
      prob_mat1<- which(is.nan(prob_mat[,1]))
      if (length(prob_mat1) > 0) {
        prob_mat<- prob_mat[-prob_mat1,]
        theta<- theta[-prob_mat1]
        signatures<- signatures[-prob_mat1]
      }
      prob_mat=t(t(prob_mat)/(1-sum(tao)+colSums(prob_mat)))
      # M step: update sample distributions.
      tao=rowMeans(prob_mat)
      for(i in 1:length(signatures)){
        theta[[i]]$mean[,1]=(exp[signatures[[i]],]%*%prob_mat[i,])/sum(prob_mat[i,])
        theta[[i]]$mean[,2]=(exp[signatures[[i]],]%*%(1-prob_mat[i,]))/sum(1-prob_mat[i,])
        keep=theta[[i]]$mean[,1]<=theta[[i]]$mean[,2]
        if(any(keep)){
          theta[[i]]$mean[keep,1]=rowMeans(exp[signatures[[i]][keep],,drop=F])
          theta[[i]]$mean[keep,2]=theta[[i]]$mean[keep,1]
        }
        tmp1=t((exp[signatures[[i]],,drop=F]-theta[[i]]$mean[,1])^2)
        tmp2=t((exp[signatures[[i]],,drop=F]-theta[[i]]$mean[,2])^2)
        diag(theta[[i]]$sigma1)=diag(theta[[i]]$sigma2)=
          colSums(tmp1*prob_mat[i,]+tmp2*(1-prob_mat[i,]))/dim(prob_mat)[2]
        diag(theta[[i]]$sigma1)[diag(theta[[i]]$sigma1)<sigma_min]=sigma_min
        diag(theta[[i]]$sigma2)[diag(theta[[i]]$sigma2)<sigma_min]=sigma_min
      }
      labels[,labels_i]=apply(rbind(1-colSums(prob_mat),prob_mat),2,which.max)-1
      # Compare estimations with stop rules.
      if(mean(apply(labels,1,function(x) length(unique(x))==1))>=convergence_rate){
        cat('Job finished successfully.',file=log_file,append=T)
        cat('\n',file=log_file,append=T)
        cat('1',file=status_file,append=F)
        break
      }
      labels_i=labels_i+1
      if(labels_i==convergence_n+1){
        labels_i=1
      }
      if(iter==max_iter){
        cat('Maximum iterations, breaking out.',file=log_file,append=T)
        cat('\n',file=log_file,append=T)
      }
    }
    #Build result matrices.
    colnames(prob_mat)=colnames(exp)
    row.names(prob_mat)=names(signatures)
    row.names(labels)=colnames(exp)
    # Attempt to remove unused signatures. 
    dummytest=sapply(1:length(signatures),function(i) mean(theta[[i]]$mean[,1]-theta[[i]]$mean[,2]==0))
    if(all(dummytest<=sensitivity_cutoff)){
      unsatisfied=0
    }else{
      rev=which(dummytest>sensitivity_cutoff);
      cat(paste('Remove dummy signatures:',rev,sep=' '),file=log_file,append=T)
      cat('\n',file=log_file,append=T)
      signatures=signatures[-rev]
      tmp=1-sum(tao)
      tao=tao[-rev]
      tao=tao/(tmp+sum(tao))
      theta=theta[-rev]
    }
  }
  return(list(cell_labels=c("unknown",names(signatures))[1+labels[,labels_i]],probabilities=prob_mat))
}


# CellMatch is used to make marker file for each dataset. (https://github.com/ZJUFanLab/scCATCH)
# select a specific species and tissue type for each dataset.
# species ('Human' or 'Mouse'); tissue (e.g.,'Blood','Brain',etc.)

CellMatch<- CellMatch[CellMatch$speciesType == species & CellMatch$tissueType %in% tissue & CellMatch$cancerType == 'Normal',]
cellname<- unique(CellMatch$cellName)
cell_marker<- list()
for (j in 1:length(cellname)) {
  cell_marker[[j]]<- CellMatch[CellMatch$cellName == cellname[j],]$geneSymbol
  names(cell_marker)[j]<- cellname[j]
}

# Single-cell transcriptomics and cell type information are both curated from literature and pre-processed generating ndata and the corresponding celltype objects.
# ndata represents the normolized dgCMatrix object using Seurat. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

test_res<- SCINA(exp = as.matrix(ndata),signatures = cell_marker)
celltype$SCINA<- test_res$cell_labels
