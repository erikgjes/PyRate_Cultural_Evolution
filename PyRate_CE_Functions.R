#LTT Plots
plot_LTT<-function(x,y){
  a<-rev(sort(x$ts))
  a<-c(a,0)
  n<-length(a)
  z<-max(a)-(a)
  vn<-1:n
  f<-log(vn+1) ~ z
  plot(f,xlab="Time \n (0 is Most Recent)", ylab="Log of Lineages",type="l",xaxt="n",col="dodgerblue",lwd=2)
  year.labels<-seq(round(min(z),digits=0),round(max(z),digits=0),by=y)
  axis(1,at=seq(0,round(max(z),digits=0),by=y),labels=rev(year.labels))
  segments(min(z),min(log(vn+1)),max(z),max(log(vn+1)),lty=2)
}

###Histogram Functions
lifespan_hist<-function(x,y){
  lifespans <-x %>% group_by(species) %>% summarize(lifespan=mean(ts)-mean(te))
  #lifespans <- lifespans %>% group_by(species,lifespan) %>% data.frame(species=species,lifespan=lifespan)
  hist(lifespans$lifespan,breaks=seq(floor(min(lifespans$lifespan)),ceiling(max(lifespans$lifespan))),
       right=TRUE,xlab="Lifespans",main="",xaxt="n")
  lifespan.labels<-seq(round(min(lifespans$lifespan)),round(max(lifespans$lifespan)),by=y)
  axis(1,at=seq(round(min(lifespans$lifespan)),round(max(lifespans$lifespan)),by=y),labels=lifespan.labels)
}

###Average Lifespan Through Time
average_lifespan<-function(x,y){
  rounded<-x %>% mutate(round_ts=round(ts,digits=0), round_te=round(te,digits=0))
  average_lifespan<-rounded %>% group_by(round_ts) %>% summarize(diff=mean(ts)-mean(te))
  z<-average_lifespan$round_ts
  plot(average_lifespan,xlab="Time \n (0 is Most Recent)", ylab="Mean Lifespan",col="dodgerblue",lwd=2,
       type="l",xaxt="n",xlim=c(max(z),min(z)))
  year.labels<-seq(floor(min(z)),ceiling(max(z)),by=y)
  axis(1,at=seq(floor(min(z)),ceiling(max(z)),by=y),labels=year.labels)
}

###Diversity Through Time
diversity_time<-function(x,y){
  x$ts<-round(x$ts,digits=0)
  x$te<-round(x$te,digits=0)
  ts.counts<-x %>% group_by(ts) %>% summarize(total_ts=n_distinct(species)) %>% mutate(cumsum = cumsum(total_ts)) %>% mutate(ID=ts)
  plot(ts.counts$ts,ts.counts$cumsum,type="l",xlab="Time \n (0 is Most Recent)",col="dodgerblue",lwd=2,
       ylab="Cumulative Diversity",xaxt="n")
  year.labels<-seq(floor(min(ts.counts$ts)),ceiling(max(ts.counts$ts)),by=y)
  axis(1,at=seq(floor(min(ts.counts$ts)),ceiling(max(ts.counts$ts)),by=y),labels=rev(year.labels))
}

summary_plots<-function(x,y){
  par(mfrow=c(2,2))
  plot_LTT(x,y)
  diversity_time(x,y)
  average_lifespan(x,y)
  lifespan_hist(x,y)
  par(mfrow=c(1,1))
}

occur_to_lifespan<-function(x){
  x.1<-x
  colnames(x.1)<-c("year","make","model")
  x.1[,2]<-str_c(tolower(str_trim(x[,2])))
  x.1[,3]<-str_c(tolower(str_trim(x[,3])))
  x.1[,2]<-str_replace_all(string=x[,2], pattern=" ", repl="_")
  x.1[,3]<-str_replace_all(string=x[,3], pattern=" ", repl="_")
  lifespan <-  x.1 %>%
    arrange(year) %>%
    group_by(make,model) %>%
    summarise(first_year=first(year), last_year=last(year))
  colnames(lifespan)<-c(names(x[,2:3]),"first_year","last_year")
  raw_lifespan<<-lifespan
  }

pyrate_fixed<-function(x){
  x.1<-x
  colnames(x.1)<-c("make","model","first_year","last_year")
  fixed_pyrate_data<<-data.frame(clade=as.numeric(factor(x.1$make)),
                               species=as.numeric(factor(x.1$model)),
                               ts=(max(x.1$last_year)-x.1$first_year),
                               te=(max(x.1$last_year)-x.1$last_year))
  fixed_make_key<-data.frame(make_name=x.1$make,make_no=as.numeric(factor(x.1$make))) %>% group_by(make_name) %>% summarize(No=first(make_no), count=n()) %>% arrange(No)
  colnames(fixed_make_key)<-c(names(x[,1]),"No","Count")
  fixed_model_key<-data.frame(model_name=x.1$model,model_no=as.numeric(factor(x.1$model))) %>% arrange(model_no)
  colnames(fixed_model_key)<-c(names(x[,2]),"No")
  write.table(fixed_make_key,"~/Desktop/PyRate_CE_Tutorial/pyrate_fixed_make_key.csv",quote=FALSE,sep=",",row.names = FALSE)
  write.table(fixed_model_key,"~/Desktop/PyRate_CE_Tutorial/pyrate_fixed_model_key.csv",quote=FALSE,sep=",",row.names = FALSE)
}

pyrate_preservation<-function(x){
  x$status<-"extinct"
  x$status[which(x$last_year>=max(x$last_year))]<-"extant"
  preservation_pyrate_data<<-data.frame(Species=x[,2],
                                        Status=x$status,
                                        min_age=max(x$last_year)-x$last_year,
                                        max_age=max(x$last_year)-x$first_year)
  write.table(preservation_pyrate_data,"~/Desktop/PyRate_CE_Tutorial/preservation_pyrate_data.txt",sep="\t",quote=FALSE,row.names=FALSE)
  }




