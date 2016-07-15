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




