############ ALL FUNCTIONS DEFINED HERE ################################################
########################################################################################
W2E <-function(x) cbind(which(x!=0,arr.ind=TRUE),x[x!=0])
########################################################################################
########################################################################################
########################################################################################
recoderFunc <- function(data, 
                        oldvalue, 
                        newvalue) {
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  newvec <- data
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  newvec
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
setup <- function (data,
                   sep,
                   header,
                   out,
                   plot,
                   ar) {
  
  files <- list.files(data, full.names=TRUE)
  subjects <- length(files)
  trackparts <- matrix(0,nrow=subjects,ncol=2)
  trackparts[,2] <- sapply(strsplit(basename(files),"\\."), 
                           function(x) paste(x[1:(length(x)-1)], collapse=".")) 
  for (p in 1:subjects) 
  {
    trackparts[p,1] <- p; file <- files[p]
    num.series = as.matrix(read.table(file,sep=sep,header=header))
    first <- as.matrix(num.series[1:(nrow(num.series)-1),])
    second <- as.matrix(num.series[2:nrow(num.series),])
    full <- cbind(first,second)
    resout <- paste(out,"/","data","_",trackparts[p,2],".txt",sep="")
    write.table(full, resout, row.names=FALSE, col.names=FALSE,quote=FALSE)
  }
  rois <- ncol(num.series)
  vars <- rois*2
  varnames <- colnames(num.series)
  varnames <- rep(varnames,2)
  cutoffind <- qchisq(.99,1)
  for (j in 1:(rois)) varnames[j] <- paste(varnames[j],"lag",sep="")
  out.files <- list.files(out,full.names=TRUE,pattern=".txt")
  lvarnames <- character()
  for (j in 1:rois) lvarnames[j] <- paste("VAR",j,"lag",sep="")
  for (j in (rois+1):(rois*2)) lvarnames[j] <- paste("VAR",j-rois,sep="")
  x <- seq(1:vars)
  y <- substring(lvarnames,4)
  betas <- paste(out,"/betas",sep="")
  dir.create(betas, showWarnings=FALSE)
  SEs <- paste(out,"/SEs",sep="")
  dir.create(SEs, showWarnings=FALSE)
  fitted <- paste(out,"/fitted",sep="")
  dir.create(fitted, showWarnings=FALSE)
  if (plot==TRUE) {plots <- paste(out,"/plots",sep="")
                   dir.create(plots,showWarnings=FALSE)
                   plot.names <-  varnames[(rois+1):(rois*2)]
  } else {plots<- ""; plot.names <-""}
  line1 <- paste(capture.output(for (i in 1:(rois*2)){ 
    cat(lvarnames[i],"=~","1*",varnames[i],sep="","\n")}),collapse="\n")
  line2 <- paste(capture.output(for (i in 1:rois) {for (j in 1:i){ 
    cat(lvarnames[i],"~~",lvarnames[j],sep="","\n")}}),collapse="\n")
  line3 <- paste(capture.output(for (i in (rois+1):(rois*2)) {
    cat (lvarnames[i],"~~",lvarnames[i],sep="","\n")}),collapse="\n")
  if (ar==TRUE) {line4 <- paste(capture.output(for (j in 1:rois) {
    cat (lvarnames[j+rois],"~",lvarnames[j],sep="","\n")}),collapse="\n")}
  if (ar==FALSE) {line4 <- paste(capture.output(for (j in 1:rois) {
    cat (lvarnames[j],"~0*",lvarnames[j+rois],sep="","\n")}),collapse="\n")}
  syntax <- paste(line1,line2,line3,line4,sep="\n")
  candidate.paths <- capture.output(
    for (i in (rois+1):(rois*2)) 
    {for (j in 1:(rois*2)) 
    {cat(lvarnames[i],"~",lvarnames[j],sep="","\n")}})
  list <- list("subjects"=subjects,
               "rois"=rois, 
               "x"=x,
               "y"=y,
               "varnames"=varnames,
               "trackparts"=trackparts,
               "out.files"=out.files,
               "vars"=vars,
               "lvarnames"=lvarnames, 
               "cutoffind"=cutoffind,
               "betas"=betas,
               "SEs"=SEs,
               "fitted"=fitted,
               "plots"=plots,
               "plot.names"=plot.names,
               "syntax"=syntax,
               "candidate.paths"=candidate.paths)
  return(list)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
fit.model <- function (varnames, 
                       syntax, 
                       data.file) {
  
  fit <- try(lavaan(syntax, 
                    data = data.file, 
                    model.type = "sem", 
                    sample.nobs = time,
                    missing = "fiml", 
                    estimator = "ml",
                    int.ov.free = TRUE,
                    int.lv.free = FALSE,
                    auto.fix.first = TRUE,
                    auto.var = TRUE,
                    auto.cov.lv.x = TRUE, 
                    auto.th = TRUE,
                    auto.delta = TRUE,
                    auto.cov.y = FALSE,
                    auto.fix.single = TRUE))
  return(fit)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
miSEM <- function (vars,
                   continue = 1,
                   varnames,
                   syntax,
                   out.files,
                   candidate.paths,
                   subjects,
                   count.group.paths) {
  param=NULL
  Freq=NULL
  mi.index <- matrix(1:((vars-1)*vars),nrow=((vars-1)*vars),ncol=1)
  cutoff <- qchisq(.95,subjects)
  cutoffgroup <- qchisq(1-.05/subjects,1)
  while (continue == 1) {
    mi.list <- matrix(0, nrow=vars*(vars-1)*subjects, ncol=6)
    colnames(mi.list) <- c("subject", "index", "lhs", "op", "rhs", "mi")
    for (k in 1:subjects)
    {
      data.file <- read.table(out.files[k],sep=" ")
      time <- nrow(data.file)
      colnames(data.file) <- c(varnames)
      fit <- fit.model(varnames = varnames,  
                       syntax = syntax,
                       data.file = data.file)
      singular <- tryCatch(modindices(fit),error=function(e) e)
      check.singular <- any(grepl("singular",singular)==TRUE)
      converge <- lavInspect(fit, "converged")
      if (check.singular==FALSE & converge==TRUE) {
        print(k)
        mi <- as.matrix(singular[singular$op == "~",])[,1:4]
      }
      # if it doesn't converge or is computationally singular
      if (converge==FALSE | check.singular==TRUE) {
        mi <- matrix(NA, nrow=((vars-1)*vars), ncol=4)}
      #stacking matrices
      mi.subject <- matrix(k,nrow=((vars-1)*vars),ncol=1)
      mi.list[(((nrow(mi)*k)-nrow(mi))+1):(nrow(mi)*k),(1:6)] <- 
        as.matrix(cbind(mi.subject,mi.index,mi),rownames.force=FALSE) 
    }  
    mi.all <- as.data.frame(mi.list[complete.cases(mi.list),])
    mi.all$param <- paste(mi.all$lhs, mi.all$op, mi.all$rhs, sep="")
    mi.all <- mi.all[-c(3:5)]
    mi.all <- subset(mi.all,param %in% candidate.paths) 
    mi.all[,3] <- as.numeric(as.character(mi.all[,3]))
    mi.all[,1] <- as.numeric(as.character(mi.all[,1]))
    mi.all <- transform(mi.all, sum = ave(mi, param, FUN=sum))
    mi.high <- subset(mi.all, mi>cutoffgroup)
    mi.count <- subset(as.data.frame(table(mi.high$param)),Freq>0)
    mi.high.count <- subset(mi.high, !duplicated(param))
    mi.merge <- merge(x=mi.high.count, y=mi.count, by.x="param", by.y="Var1")
    paramadd <- mi.merge[which.max(mi.merge$Freq),1]
    paramsum <- as.numeric(mi.merge[which.max(mi.merge$Freq),5])
    y <- subset(mi.all, param==paramadd, select=mi)$mi
    x <- seq(0, max(y), len=(subjects)) 
    drop.element <- ifelse(max(y)<cutoffgroup,TRUE,FALSE)
    n_elements <- append(0,hist(y,x,plot=FALSE)$counts)
    p <- rev(coef(lm.fit(outer(x,0:4,'^'),n_elements)))
    integratelow <- as.numeric((p[1]/(5))*cutoffgroup^(4+1)+(p[2]/4)*cutoffgroup^(3+1)+
                                 (p[3]/3)*cutoffgroup^(2+1)+(p[4]/2)*cutoffgroup^(1+1)+
                                 (p[5]/1)*cutoffgroup^(1))
    integratehigh <- as.numeric(((p[1]/(5))*max(y)^(4+1)+(p[2]/4)*max(y)^(3+1)+
                                   (p[3]/3)*max(y)^(2+1)+(p[4]/2)*max(y)^(1+1)+
                                   (p[5]/1*max(y)^(1)))-integratelow)
    
    if (integratelow >= integratehigh | drop.element==TRUE) {continue <-0 
    } else if (paramsum >= cutoff) {
      continue <- 1
      syntax <- paste(syntax,as.name(paramadd),sep="\n")
      count.group.paths <- count.group.paths +1          
    }
  }  #end of while continue>0
  list <- list("syntax"=syntax,"count.group.paths"=count.group.paths)
  return(list)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
evalbetas <- function (bad,
                       count.group.paths,
                       syntax,
                       out.files,
                       subjects,
                       rois,
                       varnames) {
  op=NULL
  sig=NULL
  param=NULL
  z=NULL
  #getting coefficients for final model
  while (bad == 1) {
    list.all <- matrix(NA, nrow=(rois+count.group.paths)*subjects, ncol=2)
    colnames(list.all) <- c("param", "z")
    for (k in 1:subjects)
    {
      data.file <- read.table(out.files[k],sep=" ")
      time <- nrow(data.file)
      colnames(data.file) <- c(varnames)
      fit <- fit.model(varnames = varnames, 
                       syntax = syntax,
                       data.file = data.file)
      singular <- tryCatch(modindices(fit),error=function(e) e)
      check.singular <- any(grepl("singular",singular)==TRUE)
      converge <- lavInspect(fit, "converged")
      if (check.singular==FALSE & converge==TRUE) {
        print(k)
        # new code to store coefficients and t values
        z.list <- subset(standardizedSolution(fit),op=="~")
        z.list$param <- paste(z.list$lhs, z.list$op, z.list$rhs, sep="")
        z.list <- z.list[c("param","z")]
        z.list <- z.list[complete.cases(z.list),]
      }
      # if it doesn't converge
      if (converge==FALSE | check.singular==TRUE) {
        z.list <- matrix(NA, nrow=(rois+count.group.paths), ncol=2) #if AR
        print("nonconvergence")
      }  
      list.all[(((nrow(z.list)*k)-nrow(z.list)+1):(nrow(z.list)*k)),] <- as.matrix(z.list)
    }
    list.all <- as.data.frame(list.all)
    list.all$z <- as.numeric(as.character(list.all$z))
    list.all$sig <- ifelse(abs(list.all$z)>1.96, 1, 0)
    list.all <- transform(list.all, sum=ave(sig, param, FUN=sum))
    if ((nrow(list.all)==0)==TRUE) bad <- 0
    paramdrop <- as.character(list.all[which.min(list.all$sum),1])
    #create histogram of t/z values
    y2 <- abs(subset(list.all, param==paramdrop, select=z)$z)  
    x2 <-  seq(0, max(y2), len=(subjects)) 
    drop.element <- ifelse(max(y2)<1.96,TRUE,FALSE)
    n_elements2 <- append(0,hist(y2,x2,plot=FALSE)$counts)
    p2 <- rev(coef(lm.fit(outer(x2,0:4,'^'),n_elements2)))
    integratelow2 <- as.numeric((p2[1]/(5))*1.96^(4+1)+(p2[2]/4)*1.96^(3+1)+
                                  (p2[3]/3)*1.96^(2+1)+(p2[4]/2)*1.96^(1+1)+
                                  (p2[5]/1)*1.96^(1))
    integratehigh2 <- as.numeric(((p2[1]/(5))*max(y2)^(4+1)+(p2[2]/4)*max(y2)^(3+1)+
                                    (p2[3]/3)*max(y2)^(2+1)+(p2[4]/2)*max(y2)^(1+1)+
                                    (p2[5]/1*max(y2)^(1)))-integratelow2)
    
    if (integratelow2 >= integratehigh2 | drop.element==TRUE) {
      syntax <- unlist(strsplit(syntax, "[\n]"))
      syntax <- syntax[!syntax %in% paramdrop]
      syntax <- paste(syntax,sep="",collapse="\n")
      bad <- 1
      count.group.paths <- count.group.paths - 1
      if ((count.group.paths==0)==TRUE) bad <- 0
    } else {bad <- 0
    }
    
  }
  list <- list("syntax"=syntax,
               "count.group.paths"=count.group.paths)
  return(list)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
setupsemi <- function(paths,
                      data,
                      sep,
                      header,
                      out,
                      plot,
                      ar){
  
  if (header==TRUE){
    for (v in 1:length(varnames)){
      paths <- gsub(varnames[v],lvarnames[v],paths)
    }
  }
  start <- paste(capture.output(cat(paths,sep="\n")),collapse="\n")
  count <- length(paths)
  files <- list.files(data, full.names=TRUE)
  subjects <- length(files)
  trackparts <- matrix(0,nrow=subjects,ncol=2)
  trackparts[,2] <- sapply(strsplit(basename(files),"\\."), 
                           function(x) paste(x[1:(length(x)-1)], collapse=".")) 
  for (p in 1:subjects) 
  {
    trackparts[p,1] <- p; file <- files[p]
    num.series = as.matrix(read.table(file,sep=sep,header=header))
    first <- as.matrix(num.series[1:(nrow(num.series)-1),])
    second <- as.matrix(num.series[2:nrow(num.series),])
    full <- cbind(first,second)
    resout <- paste(out,"/","data","_",trackparts[p,2],".txt",sep="")
    write.table(full, resout, row.names=FALSE, col.names=FALSE,quote=FALSE)
  }
  rois <- ncol(num.series)
  vars <- rois*2
  varnames <- colnames(num.series)
  varnames <- rep(varnames,2)
  cutoffind <- qchisq(.99,1)
  for (j in 1:(rois)) varnames[j] <- paste(varnames[j],"lag",sep="")
  out.files <- list.files(out,full.names=TRUE,pattern=".txt")
  lvarnames <- character()
  for (j in 1:rois) lvarnames[j] <- paste("VAR",j,"lag",sep="")
  for (j in (rois+1):(rois*2)) lvarnames[j] <- paste("VAR",j-rois,sep="")
  x <- seq(1:vars)
  y <- substring(lvarnames,4)
  betas <- paste(out,"/betas",sep=""); dir.create(betas, showWarnings=FALSE)
  SEs <- paste(out,"/SEs",sep=""); dir.create(SEs, showWarnings=FALSE)
  fitted <- paste(out,"/fitted",sep=""); dir.create(fitted, showWarnings=FALSE)
  if (plot==TRUE) {plots <- paste(out,"/plots",sep="")
                   dir.create(plots,showWarnings=FALSE)
                   plot.names <-  varnames[(rois+1):(rois*2)]
  } else {plots<- ""; plot.names <-""}
  line1 <- paste(capture.output(for (i in 1:(rois*2)){ 
    cat(lvarnames[i],"=~","1*",varnames[i],sep="","\n")}),collapse="\n")
  line2 <- paste(capture.output(for (i in 1:rois) {for (j in 1:i){ 
    cat(lvarnames[i],"~~",lvarnames[j],sep="","\n")}}),collapse="\n")
  line3 <- paste(capture.output(for (i in (rois+1):(rois*2)) {
    cat (lvarnames[i],"~~",lvarnames[i],sep="","\n")}),collapse="\n")
  if (ar==TRUE) {line4 <- paste(capture.output(for (j in 1:rois) {
    cat (lvarnames[j+rois],"~",lvarnames[j],sep="","\n")}),collapse="\n")}
  if (ar==FALSE) {line4 <- paste(capture.output(for (j in 1:rois) {
    cat (lvarnames[j],"~0*",lvarnames[j+rois],sep="","\n")}),collapse="\n")}
  syntax <- paste(line1,line2,line3,line4,start,sep="\n")
  candidate.paths <- capture.output(
    for (i in (rois+1):(rois*2)) 
    {for (j in 1:(rois*2)) 
    {cat(lvarnames[i],"~",lvarnames[j],sep="","\n")}})
  list <- list("subjects"=subjects,"rois"=rois, "x"=x,"y"=y,
               "varnames"=varnames,"trackparts"=trackparts,"out.files"=out.files,
               "vars"=vars,"lvarnames"=lvarnames, "cutoffind"=cutoffind,
               "betas"=betas,"SEs"=SEs,"fitted"=fitted,"plots"=plots,
               "plot.names"=plot.names,"syntax"=syntax,"candidate.paths"=candidate.paths,
               "count"=count)
  return(list)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
miSEMsemi <- function (vars,
                       continue=1,
                       varnames,
                       syntax,
                       out.files,
                       candidate.paths,
                       subjects,
                       count.group.paths,
                       count) {
  param=NULL
  Freq=NULL
  mi.index <- matrix(1:((vars-1)*vars),nrow=((vars-1)*vars),ncol=1)
  count.group.paths <- 0 + count
  cutoff <- qchisq(.95,subjects)
  cutoffgroup <- qchisq(1-.05/subjects,1)
  while (continue == 1) {
    mi.list <- matrix(0, nrow=vars*(vars-1)*subjects, ncol=6)
    colnames(mi.list) <- c("subject", "index", "lhs", "op", "rhs", "mi")
    for (k in 1:subjects)
    {
      data.file <- read.table(out.files[k],sep=" ")
      time <- nrow(data.file)
      colnames(data.file) <- c(varnames)
      fit <- fit.model(varnames = varnames,  
                       syntax = syntax,
                       data.file = data.file)
      singular <- tryCatch(modindices(fit),error=function(e) e)
      check.singular <- any(grepl("singular",singular)==TRUE)
      converge <- lavInspect(fit, "converged")
      if (check.singular==FALSE & converge==TRUE) {
        print(k)
        mi <- as.matrix(singular[singular$op == "~",])[,1:4]
      }
      # if it doesn't converge or is computationally singular
      if (converge==FALSE | check.singular==TRUE) {
        mi <- matrix(NA, nrow=((vars-1)*vars), ncol=4)}
      #stacking matrices
      mi.subject <- matrix(k,nrow=((vars-1)*vars),ncol=1)
      mi.list[(((nrow(mi)*k)-nrow(mi))+1):(nrow(mi)*k),(1:6)] <- 
        as.matrix(cbind(mi.subject,mi.index,mi),rownames.force=FALSE) 
    }
    
    mi.all <- as.data.frame(mi.list[complete.cases(mi.list),])
    mi.all$param <- paste(mi.all$lhs, mi.all$op, mi.all$rhs, sep="")
    mi.all <- mi.all[-c(3:5)]
    mi.all <- subset(mi.all,param %in% candidate.paths) 
    mi.all[,3] <- as.numeric(as.character(mi.all[,3]))
    mi.all[,1] <- as.numeric(as.character(mi.all[,1]))
    mi.all <- transform(mi.all, sum = ave(mi, param, FUN=sum))
    mi.high <- subset(mi.all, mi>cutoffgroup)
    mi.count <- subset(as.data.frame(table(mi.high$param)),Freq>0)
    mi.high.count <- subset(mi.high, !duplicated(param))
    mi.merge <- merge(x=mi.high.count, y=mi.count, by.x="param", by.y="Var1")
    paramadd <- mi.merge[which.max(mi.merge$Freq),1]
    paramsum <- as.numeric(mi.merge[which.max(mi.merge$Freq),5])
    y <- subset(mi.all, param==paramadd, select=mi)$mi
    x <- seq(0, max(y), len=(subjects)) 
    drop.element <- ifelse(max(y)<cutoffgroup,TRUE,FALSE)
    n_elements <- append(0,hist(y,x,plot=FALSE)$counts)
    p <- rev(coef(lm.fit(outer(x,0:4,'^'),n_elements)))
    integratelow <- as.numeric((p[1]/(5))*cutoffgroup^(4+1)+(p[2]/4)*cutoffgroup^(3+1)+
                                 (p[3]/3)*cutoffgroup^(2+1)+(p[4]/2)*cutoffgroup^(1+1)+
                                 (p[5]/1)*cutoffgroup^(1))
    integratehigh <- as.numeric(((p[1]/(5))*max(y)^(4+1)+(p[2]/4)*max(y)^(3+1)+
                                   (p[3]/3)*max(y)^(2+1)+(p[4]/2)*max(y)^(1+1)+
                                   (p[5]/1*max(y)^(1)))-integratelow)
    
    if (integratelow >= integratehigh | drop.element==TRUE) {continue <-0 
    } else if (paramsum >= cutoff) {
      continue <- 1
      syntax <- paste(syntax,as.name(paramadd),sep="\n")
      count.group.paths <- count.group.paths +1          
    }
  }  #end of while continue>0
  list <- list("syntax"=syntax,
               "count.group.paths"=count.group.paths)
  return(list)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
evalbetassemi <- function (bad,
                           count.group.paths,
                           syntax,
                           out.files,
                           subjects,
                           rois,
                           varnames,
                           paths) {
  op=NULL
  sig=NULL
  param=NULL
  z=NULL
  #getting coefficients for final model
  while (bad == 1) {
    list.all <- matrix(NA, nrow=(rois+count.group.paths)*subjects, ncol=2)
    colnames(list.all) <- c("param", "z")
    for (k in 1:subjects)
    {
      data.file <- read.table(out.files[k],sep=" ")
      time <- nrow(data.file)
      colnames(data.file) <- c(varnames)
      fit <- fit.model(varnames = varnames, 
                       syntax = syntax,
                       data.file = data.file)
      singular <- tryCatch(modindices(fit),error=function(e) e)
      check.singular <- any(grepl("singular",singular)==TRUE)
      converge <- lavInspect(fit, "converged")
      if (check.singular==FALSE & converge==TRUE) {
        print(k)
        # new code to store coefficients and t values
        z.list <- subset(standardizedSolution(fit),op=="~")
        z.list$param <- paste(z.list$lhs, z.list$op, z.list$rhs, sep="")
        z.list <- z.list[c("param","z")]
        z.list <- z.list[complete.cases(z.list),]
      }
      # if it doesn't converge
      if (converge==FALSE | check.singular==TRUE) {
        z.list <- matrix(NA, nrow=(rois+count.group.paths), ncol=2) #if AR
        print("nonconvergence")
      }  
      list.all[(((nrow(z.list)*k)-nrow(z.list)+1):(nrow(z.list)*k)),] <- as.matrix(z.list)
    }
    list.all <- as.data.frame(list.all)
    list.all$z <- as.numeric(as.character(list.all$z))
    list.all$sig <- ifelse(abs(list.all$z)>1.96, 1, 0)
    list.all <- transform(list.all, sum=ave(sig, param, FUN=sum))
    list.all <- subset(list.all,!(param %in% paths))
    if ((nrow(list.all)==0)==TRUE) bad <- 0
    paramdrop <- as.character(list.all[which.min(list.all$sum),1])
    #create histogram of t/z values
    y2 <- abs(subset(list.all, param==paramdrop, select=z)$z)  
    x2 <-  seq(0, max(y2), len=(subjects)) 
    drop.element <- ifelse(max(y2)<1.96,TRUE,FALSE)
    n_elements2 <- append(0,hist(y2,x2,plot=FALSE)$counts)
    p2 <- rev(coef(lm.fit(outer(x2,0:4,'^'),n_elements2)))
    integratelow2 <- as.numeric((p2[1]/(5))*1.96^(4+1)+(p2[2]/4)*1.96^(3+1)+
                                  (p2[3]/3)*1.96^(2+1)+(p2[4]/2)*1.96^(1+1)+
                                  (p2[5]/1)*1.96^(1))
    integratehigh2 <- as.numeric(((p2[1]/(5))*max(y2)^(4+1)+(p2[2]/4)*max(y2)^(3+1)+
                                    (p2[3]/3)*max(y2)^(2+1)+(p2[4]/2)*max(y2)^(1+1)+
                                    (p2[5]/1*max(y2)^(1)))-integratelow2)
    
    if (integratelow2 >= integratehigh2 | drop.element==TRUE) {
      syntax <- unlist(strsplit(syntax, "[\n]"))
      syntax <- syntax[!syntax %in% paramdrop]
      syntax <- paste(syntax,sep="",collapse="\n")
      bad <- 1
      count.group.paths <- count.group.paths - 1
      if ((count.group.paths==0)==TRUE) bad <- 0
    } else {bad <- 0
    }
    
  }
  list <- list("syntax"=syntax,"count.group.paths"=count.group.paths)
  return(list)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
addind <- function (done,
                    evaluate,
                    varnames,
                    syntax,
                    cutoffind,
                    data.file,
                    candidate.paths) {
  param=NULL
  while (done==0) {
    count.ind.paths <- 0
    vec.MI <- character()
    time <- nrow(data.file)
    colnames(data.file) <- c(varnames)
    fit <- fit.model(varnames = varnames,
                     syntax = syntax,
                     data.file = data.file)
    converge <- lavInspect(fit, "converged")
    singular <- tryCatch(modindices(fit),error=function(e) e)
    check.singular <- any(grepl("singular",singular)==TRUE)
    if (check.singular==FALSE & converge==TRUE) {
      indMI <- as.matrix(singular[singular$op == "~",])[,1:4]
      indMI <- as.data.frame(indMI[complete.cases(indMI),])
      indMI$param <- paste(indMI$lhs, indMI$op, indMI$rhs, sep="")
      indMI <- subset(indMI,param %in% candidate.paths)
      indMI$mi <- as.numeric(as.character(indMI$mi))
      indparamadd <- indMI[which.max(indMI$mi),5]
      indparamaddval <- indMI[which.max(indMI$mi),4]
      indices <- fitMeasures(fit,c("chisq","df","pvalue","rmsea","srmr",
                                   "nnfi","cfi"))  
      rmsea <- indices[4] ; srmr <- indices[5]
      cfi <- indices[6]; nnfi <- indices[7]
      rmseaE <- ifelse(rmsea<.05, 1, 0); srmrE <- ifelse(srmr<.05, 1, 0)
      cfiE <- ifelse(cfi>.95, 1, 0);nnfiE <- ifelse(nnfi>.95, 1, 0)
      excellent <- sum(rmseaE, srmrE, cfiE, nnfiE)
      if (excellent >= 2) {
        done <- 1 
        fixfit <- 0
      } else if (excellent <2) {
        done <- 0
        if (indparamaddval >= cutoffind) {
          syntax <- paste(syntax,as.name(indparamadd),sep="\n")
          count.ind.paths <- count.ind.paths + 1
          vec.MI <- append(vec.MI,indparamadd)
        } else {
          done <- 1 
          fixfit <- 1
        }
      }
    } else {done=1;fixfit=0}
    fit <- fit.model(varnames = varnames,
                     syntax = syntax,
                     data.file = data.file)
    converge <- lavInspect(fit, "converged")
    singular <- tryCatch(modindices(fit),error=function(e) e)
    check.singular <- any(grepl("singular",singular)==TRUE)
    if (converge==TRUE & check.singular==FALSE){
      indices <- fitMeasures(fit,c("chisq","df","pvalue","rmsea","srmr",
                                   "nnfi","cfi"))  
      rmsea <- indices[4] ; srmr <- indices[5]
      cfi <- indices[6]; nnfi <- indices[7]
      rmseaE <- ifelse(rmsea<.05, 1, 0);  srmrE <- ifelse(srmr<.05, 1, 0)
      cfiE <- ifelse(cfi>.95, 1, 0);  nnfiE <- ifelse(nnfi>.95, 1, 0)
      excellent <- sum(rmseaE, srmrE, cfiE, nnfiE)
      if (excellent >= 2) {done <- 1; fixfit <- 0}
    } else {done <- 1;  fixfit <- 0;evaluate <- 0}  
    if (converge==FALSE) done <- 1
    if (check.singular==TRUE) done <-1
  }
  if (count.ind.paths==0) evaluate <- 0
  list <- list("evaluate"=evaluate,"fixfit"=fixfit,"syntax"=syntax,"vec.MI"=vec.MI)
  return(list)
}  
########################################################################################
########################################################################################
########################################################################################
########################################################################################
evalind <- function (evaluate,
                     varnames,
                     syntax,
                     data.file,
                     fixfit,
                     vec.MI) {
  op=NULL
  param=NULL
  while (evaluate==1) {
    fit <- fit.model(varnames = varnames, 
                     syntax = syntax,
                     data.file = data.file)  
    converge <- lavInspect(fit, "converged")
    singular <- tryCatch(modindices(fit),error=function(e) e)
    check.singular <- any(grepl("singular",singular)==TRUE)
    if (check.singular==FALSE & converge==TRUE) {
      indlist <- subset(standardizedSolution(fit),op=="~")
      indlist$param <- paste(indlist$lhs, indlist$op, indlist$rhs, sep="")
      indlist <- indlist[c("param","z")]
      indlist <- indlist[complete.cases(indlist),]
      indlist <- as.data.frame(indlist)
      indlist$z <- as.numeric(as.character(indlist$z))
      indlist <- subset(indlist,param %in% vec.MI)
      parampruneposs <- as.character(indlist[which.min(indlist$z),1])
      parampruneval <- as.numeric(indlist[which.min(indlist$z),2])
      prune <- ifelse(parampruneval>1.96, 0, 1)
      if (nrow(indlist)==0) prune <- 0
      if (prune == 1) {
        syntax <- unlist(strsplit(syntax, "[\n]"))
        syntax <- syntax[!syntax %in% parampruneposs]
        syntax <- paste(syntax,sep="",collapse="\n")
        paramprune <- parampruneposs
        evaluate <- 1
      } else {evaluate <- 0;fixfit <- 0}
      fit <- fit.model(varnames=varnames,  
                       syntax=syntax,data.file=data.file) 
      converge <- lavInspect(fit, "converged")
      singular <- tryCatch(modindices(fit),error=function(e) e)
      check.singular <- any(grepl("singular",singular)==TRUE)
      if (converge==TRUE & check.singular==FALSE){
        indices <- fitMeasures(fit,c("chisq","df","pvalue","rmsea","srmr",
                                     "nnfi","cfi"))  
        rmsea <- indices[4] ; srmr <- indices[5]
        cfi <- indices[6]; nnfi <- indices[7]
        rmseaE <- ifelse(rmsea<.05, 1, 0);srmrE <- ifelse(srmr<.05, 1, 0)
        cfiE <- ifelse(cfi>.95, 1, 0);nnfiE <- ifelse(nnfi>.95, 1, 0)
        excellent <- sum(rmseaE, srmrE, cfiE, nnfiE)
        if (excellent >= 2) fixfit<-0 else fixfit<-1
        # if it doesn't converge
      }
      if (converge==FALSE) {evaluate <- 0; fixfit <- 0}
      if (check.singular==TRUE) {evaluate <- 0; fixfit <- 0}
    }
  }
  list <- list("fixfit"=fixfit,"syntax"=syntax)
  return(list)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
fixfitind <- function (fixfit,
                       varnames,
                       syntax,
                       data.file,
                       candidate.paths) {
  param=NULL
  while (fixfit==1) {
    fit <- fit.model(varnames = varnames,  
                     syntax = syntax,
                     data.file = data.file)
    converge <- lavInspect(fit, "converged")
    singular <- tryCatch(modindices(fit),error=function(e) e)
    check.singular <- any(grepl("singular",singular)==TRUE)
    if (check.singular==FALSE & converge==TRUE) {
      indMI <- as.matrix(singular[singular$op == "~",])[,1:4]
      indMI <- as.data.frame(indMI[complete.cases(indMI),])
      indMI$param <- paste(indMI$lhs, indMI$op, indMI$rhs, sep="")
      indMI <- subset(indMI,param %in% candidate.paths)
      indMI$mi <- as.numeric(as.character(indMI$mi))
      indparamadd <- indMI[which.max(indMI$mi),5]
      syntax <- paste(syntax,as.name(indparamadd),sep="\n")     
      fit <- fit.model(varnames = varnames,  
                       syntax = syntax,
                       data.file = data.file)
      converge <- lavInspect(fit, "converged")
      singular <- tryCatch(modindices(fit),error=function(e) e)
      check.singular <- any(grepl("singular",singular)==TRUE)
      if (check.singular==FALSE & converge==TRUE) {
        indices <- fitMeasures(fit,c("chisq","df","pvalue","rmsea","srmr",
                                     "nnfi","cfi"))  
        rmsea <- indices[4] ; srmr <- indices[5]
        cfi <- indices[6]; nnfi <- indices[7]
        df <- indices[2]
        rmseaE <- ifelse(rmsea<.05, 1, 0);  srmrE <- ifelse(srmr<.05, 1, 0)
        cfiE <- ifelse(cfi>.95, 1, 0); nnfiE <- ifelse(nnfi>.95, 1, 0)
        if (df == 0) fixfit <- 0
        excellent <- sum(rmseaE, srmrE, cfiE, nnfiE)
        if (excellent >= 2) fixfit<-0 else fixfit<-1    
      }  
      # if it doesn't converge
      if (converge==FALSE) fixfit <- 0
      if (check.singular==TRUE) fixfit <- 0 
    }
  }
  return(syntax)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
final.fit <- function(varnames,
                      syntax,
                      data.file,
                      plot,
                      trackparts,
                      rois,
                      plot.names,
                      k,
                      agg,
                      x,
                      y,
                      betas,
                      SEs,
                      plots,
                      fitted){
  op=NULL
  fit <- fit.model(varnames = varnames,
                   syntax = syntax,
                   data.file = data.file)
  
  converge <- lavInspect(fit, "converged")
  if (converge==TRUE) last.converge <- TRUE
  singular <- tryCatch(modindices(fit),error=function(e) e)
  check.singular <- any(grepl("singular",singular)==TRUE)
  
  ## update code here to remove last element
  if (converge==FALSE) {
    syntax <- unlist(strsplit(syntax, "[\n]"))
    syntax <- syntax[-length(syntax)]
    syntax <- paste(syntax,sep="",collapse="\n")
    fit <- fit.model(varnames = varnames,  
                     syntax = syntax,
                     data.file = data.file)
    converge <- lavInspect(fit, "converged")
    singular <- tryCatch(modindices(fit),error=function(e) e)
    check.singular <- any(grepl("singular",singular)==TRUE)
    last.converge <- FALSE
  }
  
  ind.fit <- matrix(NA,nrow=1,ncol=9)
  if (converge==TRUE & check.singular==FALSE) {
    indices <- fitMeasures(fit,c("chisq","df","pvalue","rmsea","srmr",
                                 "nnfi","cfi"))  
    rmsea <- indices[4] ; srmr <- indices[5]
    cfi <- indices[6]; nnfi <- indices[7]
    chisq <- indices[1]; df <- indices[2];  pval <- indices[3]
    rmseaE <- ifelse(rmsea<.05, 1, 0)
    srmrE <- ifelse(srmr<.05, 1, 0)
    cfiE <- ifelse(cfi>.95, 1, 0)
    nnfiE <- ifelse(nnfi>.95, 1, 0)
    # insert indfit
    ind.fit[1,2] <- round(chisq,digits=4)
    ind.fit[1,3] <- df; ind.fit[1,4] <- round(pval,digits=4)
    ind.fit[1,5] <- round(rmsea,digits=4); ind.fit[1,6] <- round(srmr,digits=4)
    ind.fit[1,7] <- round(nnfi,digits=4); ind.fit[1,8] <- round(cfi,digits=4)
    if (last.converge==FALSE) {ind.fit[1,9] <- "last known convergence"
    } else {ind.fit[1,9] <- "converged normally"}
    
    indlist <- subset(standardizedSolution(fit),op=="~")
    indlist$param <- paste(indlist$lhs, indlist$op, indlist$rhs, sep="")
    indlist <- indlist[complete.cases(indlist),]
    indlist <- as.data.frame(indlist)
    indsubject <- matrix(k,nrow=(nrow(indlist)), ncol=1)
    colnames(indsubject) <- c("subject")
    indlist <- cbind(indsubject, indlist)
    #creating individual-level beta matrices
    individual.paths <- matrix(0,nrow=(rois*2), ncol=(rois*2))
    individual.SEs <- matrix(0,nrow=(rois*2), ncol=(rois*2))
    indlist$row <- substring(indlist$lhs,4)
    indlist$col <- substring(indlist$rhs,4)
    indrows <- indlist$row
    indcols <- indlist$col
    indcols <- as.numeric(recoderFunc(indcols,y,x))
    indrows <- as.numeric(recoderFunc(indrows,y,x))
    indbetas <- as.numeric(as.character(indlist$est.std))
    indSEs <- as.numeric(as.character(indlist$se))
    indelements <- nrow(indlist)
    for (s in 1:indelements){
      individual.paths[indrows[s], indcols[s]] <- indbetas[s]
    }
    indoutput <- paste (betas,"/",trackparts[k,2],".csv",sep="")
    if (agg==TRUE) indoutput <- paste (betas,"/all.csv",sep="")
    individual.paths[is.na(individual.paths)] <- 0 
    write.table(individual.paths, file=indoutput, sep=",", 
                row.names=FALSE, col.names=FALSE)
    for (q in 1:indelements){
      individual.SEs[indrows[q], indcols[q]] <- indSEs[q]
    }
    seoutput <- paste (SEs,"/",trackparts[k,2],".csv",sep="")
    if (agg==TRUE) seoutput <- paste (SEs,"/all.csv",sep="")
    individual.SEs[is.na(individual.SEs)] <- 0 
    write.table(individual.SEs, file=seoutput, sep=",", 
                row.names=FALSE, col.names=FALSE) 
    
    individual.paths.dich <- matrix(0,nrow=(rois*2), ncol=(rois*2))
    indbetas.dich <- ifelse(indbetas!=0,1,0)
    for (r in 1:indelements){
      individual.paths.dich[indrows[r], indcols[r]] <- indbetas.dich[r]
    }
    inddichoutput <- paste (fitted,"/",trackparts[k,2],".txt",sep="")
    if (agg==TRUE) inddichoutput <- paste (fitted,"/all.csv",sep="")
    individual.paths.dich[is.na(individual.paths.dich)] <- 0 
    write.table(individual.paths.dich, file=inddichoutput, sep=",", 
                row.names=FALSE, col.names=FALSE)
    
    if (plot==TRUE){
      individual.paths.t <- t(individual.paths)
      Lagged <- individual.paths.t[1:(rois),(rois+1):(rois*2)]
      Contemporaneous <- individual.paths.t[(rois+1):(rois*2),(rois+1):(rois*2)]
      eLagged <- W2E(Lagged)
      eContemporaneous <- W2E(Contemporaneous)
      isLagged <- c(rep(TRUE,nrow(eLagged)), rep(FALSE,nrow(eContemporaneous)))
      plotind <- paste(plots,"/",trackparts[k,2],".pdf",sep="")
      if (agg==TRUE) plotind <- paste(plots,"/all.pdf",sep="")
      pdf(plotind)
      qgraph(rbind(eLagged,eContemporaneous),
             layout="circle", lty = ifelse(isLagged,2, 1),
             edge.labels=F, curve = TRUE, fade=FALSE,
             posCol="red",negCol="blue", labels=plot.names,
             label.cex=3, edge.label.cex=1.5,edge.label.position=.3)
      dev.off()
    } 
  }
  if (converge==FALSE) {
    ind.fit[1,9] <- "nonconvergence"
    indlist <- data.frame()
  }
  if (check.singular==TRUE) {
    ind.fit[1,9] <- "computationally singular"
    indlist <- data.frame()
  }
  list <- list("ind.fit"=ind.fit,"ind.elements"=indlist)
  return(list)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
indsem.internal <- function(subjects,
                            varnames,
                            trackparts,
                            evalbetas.out,
                            setup.out,
                            plot){
  
  all.elements <- data.frame()
  all.fit <- matrix(NA, nrow=subjects, ncol=9)
  colnames(all.fit) <- c("subject", "chisq", "df", "pval", 
                         "rmsea", "srmr", "nnfi", "cfi", "status")
  for (k in 1:subjects) {
    
    print(k)
    
    data.file <- read.table(setup.out$out.files[k],sep=" ")
    syntax <- evalbetas.out$syntax  
    
    
    colnames(data.file) <- c(varnames)
    
    addind.out <- addind(done = 0,
                         evaluate = 1,
                         varnames = setup.out$varnames,
                         syntax = syntax,
                         cutoffind = setup.out$cutoffind,
                         data.file = data.file,
                         candidate.paths = setup.out$candidate.paths)
    
    evalind.out <- evalind(evaluate = addind.out$evaluate,
                           varnames = setup.out$varnames,
                           syntax = addind.out$syntax,
                           data.file = data.file,
                           fixfit = addind.out$fixfit,
                           vec.MI = addind.out$vec.MI)
    
    fixfitind.out <- fixfitind(fixfit = evalind.out$fixfit,
                               varnames = setup.out$varnames,
                               syntax = evalind.out$syntax,
                               data.file = data.file,
                               candidate.paths = setup.out$candidate.paths)
    
    final.fit.out <- final.fit(varnames = setup.out$varnames,
                               syntax = fixfitind.out,
                               data.file = data.file,
                               plot = plot,
                               trackparts = setup.out$trackparts,
                               rois = setup.out$rois,
                               plot.names = setup.out$plot.names,
                               k = k,
                               agg = FALSE,
                               x = setup.out$x,
                               y = setup.out$y,
                               betas = setup.out$betas,
                               SEs = setup.out$SEs,
                               plots = setup.out$plots,
                               fitted = setup.out$fitted)
    
    all.elements <-rbind(all.elements,final.fit.out$ind.elements)
    all.fit[k,] <- as.matrix(final.fit.out$ind.fit)
    all.fit[,1] <- trackparts[,2]
    
  }
  list <- list("all.elements"=all.elements,"all.fit"=all.fit)
  return(list)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
wrapup <- function(all.elements,
                   all.fit,
                   rois,
                   out,
                   plot,
                   plot.names,
                   subjects,
                   plots,
                   agg,
                   x,
                   y,
                   header,
                   varnames,
                   lvarnames){
  present=NULL
  sig=NULL
  est.std=NULL
  param=NULL
  all.elements.final <- as.data.frame(all.elements)
  all.elements.final$est.std <- as.numeric(as.character(all.elements.final$est.std))
  if (agg==FALSE) {
    all.elements.final <- transform(all.elements, 
                                    mean.beta = (ave(est.std, param, FUN=sum))/subjects)
  }
  all.elements.final$sig <- ifelse(abs(all.elements.final$z)>1.96, 1, 0)
  all.elements.final$present <- ifelse(all.elements.final$sig<2,1,0)
  all.elements.final <- transform(all.elements.final, sum.sig = ave(sig, param, FUN=sum))
  all.elements.final <- transform(all.elements.final, count = ave(present, param, FUN=sum))
  all.elements.final <- subset(all.elements.final, !duplicated(param))
  
  if (agg==FALSE){
    final.paths <- matrix(0,nrow=(rois*2), ncol=(rois*2))
    final.betas <- matrix(0,nrow=(rois*2),ncol=(rois*2))
    final.paths.dich <- matrix(0,nrow=(rois*2), ncol=(rois*2))
    elements <- nrow(all.elements.final)
    rows <- substring(all.elements.final$lhs, 4)
    cols <- substring(all.elements.final$rhs, 4)
    cols <- as.numeric(recoderFunc(cols,y,x))
    rows <- as.numeric(recoderFunc(rows,y,x))
    counts <- as.numeric(all.elements.final$count)
    beta.weights <- as.numeric(all.elements.final$mean.beta)
    for (m in 1:elements){
      final.paths[rows[m], cols[m]] <- counts[m]
      final.betas[rows[m], cols[m]] <- beta.weights[m]
    }
    final.paths[is.na(final.paths)] <- 0 
    finalpaths <- paste (out,"/","finalpaths.csv",sep="")
    write.table(final.paths, file=finalpaths, sep=",", row.names=FALSE, col.names=FALSE)
    final.betas[is.na(final.betas)] <- 0 
    finalbetas <- paste (out,"/","finalbetas.csv",sep="")
    write.table(final.betas, file=finalbetas, sep=",", row.names=FALSE, col.names=FALSE)
    all.elements.summary <- all.elements.final[c("param","mean.beta","sum.sig","count")]
    
    if (plot==TRUE){
      final.paths.t <- t(final.paths)
      final.paths.t <- final.paths.t/subjects
      Lagged <- final.paths.t[1:(rois),(rois+1):(rois*2)]
      Contemporaneous <- final.paths.t[(rois+1):(rois*2),(rois+1):(rois*2)]
      eLagged <- W2E(Lagged)
      eContemporaneous <- W2E(Contemporaneous)
      isLagged <- c(rep(TRUE,nrow(eLagged)), rep(FALSE,nrow(eContemporaneous)))
      plotfinal <- paste(plots,"/finalpaths.pdf",sep="")
      pdf(plotfinal)
      qgraph(rbind(eLagged,eContemporaneous),
             layout="circle", lty = ifelse(isLagged,2, 1),
             edge.labels=FALSE, curve = TRUE,
             labels=plot.names, fade=FALSE,
             label.cex=3, edge.label.cex=1.5,edge.label.position=.3)
      dev.off()
    }
  }
  
  allfit <- paste(out,"/","allfit.csv",sep="")
  write.table(all.fit,file=allfit,sep=",",row.names=FALSE)
  all.elements <- all.elements[c("subject","param","est.std","se","pvalue","z")]
  if (agg==TRUE) {all.elements[,1] <- "all"; all.elements.summary <- data.frame()}

  if (header==TRUE){
    for (v in 1:length(varnames)){
      all.elements$param <- gsub(lvarnames[v],varnames[v],all.elements$param)
      if (agg==FALSE) all.elements.summary$param <- gsub(lvarnames[v],varnames[v],all.elements.summary$param)
    }
  }
  row.names(all.elements.summary) <- NULL
  row.names(all.elements) <- NULL
  all.fit <- as.data.frame(all.fit)
  list <- list("all.elements.summary"=all.elements.summary,"all.elements"=all.elements,
               "all.fit"=all.fit)
  return(list)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
gimme <- function(data,
                  sep,
                  header,
                  out,
                  plot = TRUE,
                  ar = FALSE){
  
  setup.out <- setup(data = data, 
                     sep = sep,
                     header = header, 
                     out = out,
                     plot = plot,
                     ar = ar)
  
  miSEM.out <- miSEM(vars = setup.out$vars,
                     continue = 1,
                     varnames = setup.out$varnames,
                     syntax = setup.out$syntax,
                     out.files = setup.out$out.files,
                     candidate.paths = setup.out$candidate.paths,
                     subjects = setup.out$subjects,
                     count.group.paths = 0)
  
  evalbetas.out <- evalbetas(bad = ifelse(identical(setup.out$syntax,miSEM.out$syntax),0,1),
                             count.group.paths = miSEM.out$count.group.paths,
                             syntax = miSEM.out$syntax,
                             out.files = setup.out$out.files, 
                             subjects = setup.out$subjects,
                             rois = setup.out$rois,
                             varnames = setup.out$varnames)
  
  indsem.internal.out <- indsem.internal(subjects = setup.out$subjects,
                                         varnames = setup.out$varnames,
                                         trackparts = setup.out$trackparts,
                                         setup.out = setup.out,
                                         evalbetas.out = evalbetas.out,
                                         plot = plot)
  
  wrapup.out <- wrapup(all.elements = indsem.internal.out$all.elements,
                       all.fit = indsem.internal.out$all.fit,
                       rois = setup.out$rois,
                       out = out,
                       plot = plot,
                       plot.names = setup.out$plot.names,
                       subjects = setup.out$subjects,
                       plots = setup.out$plots,
                       agg = FALSE,
                       x = setup.out$x,
                       y = setup.out$y,
                       header = header,
                       varnames = setup.out$varnames,
                       lvarnames = setup.out$lvarnames)
  
  return(wrapup.out)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
semigimme <- function(paths,
                      data,
                      sep,
                      header,
                      out,
                      plot=TRUE,
                      ar=FALSE){
  
  setup.out <- setupsemi(paths = paths,
                         data = data, 
                         sep = sep, 
                         header = header, 
                         out = out,
                         plot = plot,
                         ar = ar)
  
  miSEM.out <- miSEMsemi(vars = setup.out$vars,
                         continue = 1,
                         varnames = setup.out$varnames,
                         syntax = setup.out$syntax,
                         out.files = setup.out$out.files,
                         candidate.paths = setup.out$candidate.paths,
                         subjects = setup.out$subjects,
                         count.group.paths = 0,
                         count = setup.out$count)
  
  evalbetas.out <- evalbetassemi(bad = ifelse(identical(setup.out$syntax,miSEM.out$syntax),0,1),
                                 count.group.paths = miSEM.out$count.group.paths,
                                 syntax = miSEM.out$syntax,
                                 out.files = setup.out$out.files, 
                                 subjects = setup.out$subjects,
                                 rois = setup.out$rois,
                                 varnames = setup.out$varnames,
                                 paths = paths)
  
  indsem.internal.out <- indsem.internal(subjects = setup.out$subjects,
                                         varnames = setup.out$varnames,
                                         trackparts = setup.out$trackparts,
                                         setup.out = setup.out,
                                         evalbetas.out = evalbetas.out,
                                         plot = plot)
  
  wrapup.out <- wrapup(all.elements = indsem.internal.out$all.elements,
                       all.fit = indsem.internal.out$all.fit,
                       rois = setup.out$rois,
                       out = out,
                       plot = plot,
                       plot.names = setup.out$plot.names,
                       subjects = setup.out$subjects,
                       plots = setup.out$plots,
                       agg = FALSE,
                       x = setup.out$x,
                       y = setup.out$y,
                       header = header,
                       varnames = setup.out$varnames,
                       lvarnames = setup.out$lvarnames)
  
  return(wrapup.out)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
indsem.standalone <- function(subjects,
                              varnames,
                              trackparts,
                              syntax,
                              out.files,
                              setup.out,
                              plot){
  
  all.elements <- data.frame()
  all.fit <- matrix(NA, nrow=subjects, ncol=9)
  colnames(all.fit) <- c("subject", "chisq", "df", "pval", 
                         "rmsea", "srmr", "nnfi", "cfi", "status")
  for (k in 1:subjects) {
    
    print(k)
    data.file <- read.table(out.files[k],sep=" ")    
    colnames(data.file) <- c(varnames)
    
    addind.out <- addind(done = 0,
                         evaluate = 1,
                         varnames = setup.out$varnames,
                         syntax = syntax,
                         cutoffind = setup.out$cutoffind,
                         data.file = data.file,
                         candidate.paths = setup.out$candidate.paths)
    
    evalind.out <- evalind(evaluate = addind.out$evaluate,
                           varnames = setup.out$varnames,
                           syntax = addind.out$syntax,
                           data.file = data.file,
                           fixfit = addind.out$fixfit,
                           vec.MI = addind.out$vec.MI)
    
    fixfitind.out <- fixfitind(fixfit = evalind.out$fixfit,
                               varnames = setup.out$varnames,
                               syntax = evalind.out$syntax,
                               data.file = data.file,
                               candidate.paths = setup.out$candidate.paths)
    
    final.fit.out <- final.fit(varnames = setup.out$varnames,
                               syntax = fixfitind.out,
                               data.file = data.file,
                               plot = plot,
                               trackparts = setup.out$trackparts,
                               rois = setup.out$rois,
                               plot.names = setup.out$plot.names,
                               k = k,
                               agg = FALSE,
                               x = setup.out$x,
                               y = setup.out$y,
                               betas = setup.out$betas,
                               SEs = setup.out$SEs,
                               plots = setup.out$plots,
                               fitted = setup.out$fitted)
    
    
    all.elements <-rbind(all.elements,final.fit.out$ind.elements)
    all.fit[k,] <- as.matrix(final.fit.out$ind.fit)
    all.fit[,1] <- trackparts[,2]
  }
  list <- list("all.elements"=all.elements,"all.fit"=all.fit)
  return(list)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
indSEM <- function(data,
                   sep,
                   header,
                   out,
                   ar=FALSE,
                   plot=TRUE){
  
  setup.out <- setup(data = data, 
                     sep = sep, 
                     header = header, 
                     out = out,
                     plot = plot,
                     ar = ar)
  
  indsem.standalone.out <- indsem.standalone(setup.out = setup.out,
                                             subjects = setup.out$subjects,
                                             varnames = setup.out$varnames,
                                             trackparts = setup.out$trackparts,
                                             syntax = setup.out$syntax,
                                             out.files = setup.out$out.files,
                                             plot = plot)
  
  wrapup.out <- wrapup(all.elements = indsem.standalone.out$all.elements,
                       all.fit = indsem.standalone.out$all.fit,
                       rois = setup.out$rois,
                       out = out,
                       plot = plot,
                       plot.names = setup.out$plot.names,
                       subjects = setup.out$subjects,
                       plots = setup.out$plots,
                       agg = FALSE,
                       x = setup.out$x,
                       y = setup.out$y,
                       header = header,
                       varnames = setup.out$varnames,
                       lvarnames = setup.out$lvarnames)
  return(wrapup.out)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
agg.internal <- function(setup.out,
                         subjects,
                         varnames,
                         trackparts,
                         syntax,
                         out.files,
                         plot){
  
  data.all <- data.frame()
  for (k in 1:subjects){
    data.file <- read.table(out.files[k],sep=" ") 
    data.all <- rbind(data.all,data.file)
  }
  colnames(data.all) <- c(varnames)
  
  addind.out <- addind(done = 0,
                       evaluate = 1,
                       varnames = setup.out$varnames,
                       syntax = syntax,
                       cutoffind = setup.out$cutoffind,
                       data.file = data.all,
                       candidate.paths = setup.out$candidate.paths)
  
  evalind.out <- evalind(evaluate = addind.out$evaluate,
                         varnames = setup.out$varnames,
                         syntax = addind.out$syntax,
                         data.file = data.all,
                         fixfit = addind.out$fixfit,
                         vec.MI = addind.out$vec.MI)
  
  fixfitind.out <- fixfitind(fixfit = evalind.out$fixfit,
                             varnames = setup.out$varnames,
                             syntax = evalind.out$syntax,
                             data.file = data.all,
                             candidate.paths = setup.out$candidate.paths)
  
  final.fit.out <- final.fit(varnames = setup.out$varnames,
                             syntax = fixfitind.out,
                             data.file = data.all,
                             plot = plot,
                             trackparts = setup.out$trackparts,
                             rois = setup.out$rois,
                             plot.names = setup.out$plot.names,
                             k = k,
                             agg = TRUE,
                             x = setup.out$x,
                             y = setup.out$y,
                             betas = setup.out$betas,
                             SEs = setup.out$SEs,
                             plots = setup.out$plots,
                             fitted = setup.out$fitted)
  
  
  all.elements <- final.fit.out$ind.elements
  all.fit <- as.matrix(final.fit.out$ind.fit)
  all.fit[,1] <- "all"
  colnames(all.fit) <- c("subject", "chisq", "df", "pval", 
                         "rmsea", "srmr", "nnfi", "cfi", "status") 
  row.names(all.elements) <- NULL
  list <- list("all.elements"=all.elements,"all.fit"=all.fit)
  return(list)
}
########################################################################################
########################################################################################
########################################################################################
########################################################################################
aggregate <- function(data,
                      sep,
                      header,
                      out,
                      ar = FALSE,
                      plot = TRUE){
  
  setup.out <- setup(data = data, 
                     sep = sep, 
                     header = header, 
                     out = out,
                     plot = plot,
                     ar = ar)
  
  agg.internal.out <- agg.internal(setup.out = setup.out,
                                   subjects = setup.out$subjects,
                                   varnames = setup.out$varnames,
                                   trackparts = setup.out$trackparts,
                                   syntax = setup.out$syntax,
                                   out.files = setup.out$out.files,
                                   plot = plot)
  
  wrapup.out <- wrapup(all.elements = agg.internal.out$all.elements,
                       all.fit = agg.internal.out$all.fit,
                       rois = setup.out$rois,
                       out = out,
                       plot = plot,
                       plot.names = setup.out$plot.names,
                       subjects = setup.out$subjects,
                       plots = setup.out$plots,
                       agg = TRUE,
                       x = setup.out$x,
                       y = setup.out$y,
                       header = header,
                       varnames = setup.out$varnames,
                       lvarnames = setup.out$lvarnames)
  
  list <- list("all.elements"=wrapup.out$all.elements,"all.fit"=wrapup.out$all.fit)
  return(list)
  
}
########################################################################################
########################################################################################
############ END OF DEFINED FUNCTIONS ##################################################
