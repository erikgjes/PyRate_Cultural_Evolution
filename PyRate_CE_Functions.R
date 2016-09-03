lineage_pyrate<-function(x){
  x.1<-x
  colnames(x.1)<-c("clade","species","min_age","max_age")
  species_pyrate_data<<-data.frame(clade=as.numeric(factor(x.1$clade)),
                               species=as.numeric(factor(x.1$species)),
                               ts=(max(x.1$max_age)-x.1$min_age),
                               te=(max(x.1$max_age)-x.1$max_age))
                               
                               
  fixed_clade_key<-data.frame(clade_name=x.1$clade,clade_no=as.numeric(factor(x.1$clade))) %>% 
  group_by(clade_name) %>% 
  summarize(No=first(clade_no), count=n()) %>% 
  arrange(No)
  colnames(fixed_clade_key)<-c("Names","ID Number","Count")
  write.table(fixed_clade_key,"~/Desktop/PyRate_CE_Tutorial/lineage_clade_key.csv",quote=FALSE,sep=",",row.names=FALSE)
  fixed_species_key<-data.frame(species_name=x.1$species,species_no=as.numeric(factor(x.1$species))) %>% arrange(species_no)
  colnames(fixed_species_key)<-c("Names","ID Number")
  write.table(fixed_species_key,"~/Desktop/PyRate_CE_Tutorial/lineage_species_key.csv",quote=FALSE,sep=",",row.names=FALSE)
}

occurrence_pyrate<-function(x){
	x$status<-"extinct"
	colnames(x)<-c("species","min_age","max_age","trait","status")
	x$status[which(x$min_age==0)]<-"extant"
	#x$number<-as.numeric(factor(x$species))
	#x <- x %>% group_by(species) %>% mutate(new2 = if(n( ) > 1) {paste(species, row_number( ),sep="_")} else {paste0(species)})
  	occurrence_pyrate_data<<-data.frame(Species=x$species,
                                        Status=x$status,
                                        min_age=x$min_age,
                                        max_age=x$max_age,
                                        trait=x$trait)
    write.table(occurrence_pyrate_data,"~/Desktop/PyRate_CE_Tutorial/occurrence_pyrate_data.txt",sep="\t",quote=FALSE,row.names=FALSE)
}


#LTT Plots
plot_LTT<-function(x,y){
  a<-rev(sort(x$ts))
  a<-c(a,0)
  n<-length(a)
  z<-max(a)-(a)
  vn<-1:n
  f<-log(vn+1) ~ z
  plot(f,xlab="Time \n (0 is Most Recent)", ylab="Log of speciess",type="l",xaxt="n",col="dodgerblue",lwd=2,
       main="Log species \n Through Time")
  year.labels<-seq(round(min(z),digits=0),round(max(z),digits=0),by=y)
  axis(1,at=seq(0,round(max(z),digits=0),by=y),labels=rev(year.labels))
  segments(min(z),min(log(vn+1)),max(z),max(log(vn+1)),lty=2)
}

plot_LTT_occurrence<-function(x,y){
  a<-rev(sort(x$max_age))
  a<-c(a,0)
  n<-length(a)
  z<-max(a)-(a)
  vn<-1:n
  f<-log(vn+1) ~ z
  plot(f,xlab<-"Time \n (0 is Most Recent)", ylab="Log of speciess",type="l",xaxt="n",col="dodgerblue",lwd=2,
       main="Log species \n Through Time")
  year.labels<-seq(round(min(z),digits=0),round(max(z),digits=0),by=y)
  axis(1,at<-seq(0,round(max(z),digits=0),by=y),labels=rev(year.labels))
  segments(min(z),min(log(vn+1)),max(z),max(log(vn+1)),lty=2)
}

###Histogram Functions
lifespan_hist<-function(x,y){
  lifespans <-x %>% group_by(species) %>% summarize(lifespan=mean(ts)-mean(te))
  #lifespans <- lifespans %>% group_by(species,lifespan) %>% data.frame(species=species,lifespan=lifespan)
  hist(lifespans$lifespan,breaks=seq(floor(min(lifespans$lifespan)),ceiling(max(lifespans$lifespan))),
       right=TRUE,xlab="Lifespans",main="Histogram of \n Lifespans",xaxt="n")
  lifespan.labels<-seq(round(min(lifespans$lifespan)),round(max(lifespans$lifespan)),by=y)
  axis(1,at=seq(round(min(lifespans$lifespan)),round(max(lifespans$lifespan)),by=y),labels=lifespan.labels)
}

lifespan_hist_occurrence<-function(x,y){
  lifespans <-x %>% group_by(Species) %>% summarize(lifespan=mean(max_age)-mean(min_age))
  #lifespans <- lifespans %>% group_by(species,lifespan) %>% data.frame(species=species,lifespan=lifespan)
  hist(lifespans$lifespan,breaks=seq(floor(min(lifespans$lifespan)),ceiling(max(lifespans$lifespan))),
       right=TRUE,xlab="Lifespans",main="Histogram of Species \n Lifespans",xaxt="n")
  lifespan.labels<-seq(round(min(lifespans$lifespan)),round(max(lifespans$lifespan)),by=y)
  axis(1,at=seq(round(min(lifespans$lifespan)),round(max(lifespans$lifespan)),by=y),labels=lifespan.labels)
}

###Average Lifespan Through Time
average_lifespan<-function(x,y){
  rounded<-x %>% mutate(round_ts=round(ts,digits=0), round_te=round(te,digits=0))
  average_lifespan<-rounded %>% group_by(round_ts) %>% summarize(diff=mean(ts)-mean(te))
  z<-average_lifespan$round_ts
  plot(average_lifespan,xlab="Time \n (0 is Most Recent)", ylab="Mean Lifespan",col="dodgerblue",lwd=2,
       type="l",xaxt="n",main="Average Lifespan \n through Time",xlim=c(max(z),min(z)))
  year.labels<-seq(floor(min(z)),ceiling(max(z)),by=y)
  axis(1,at=seq(floor(min(z)),ceiling(max(z)),by=y),labels=year.labels)
}

average_lifespan_occurrence<-function(x,y){
  rounded<-x %>% mutate(round_max_age=round(max_age,digits=0), round_min_age=round(min_age,digits=0))
  average_lifespan<-rounded %>% group_by(round_max_age) %>% summarize(diff=mean(max_age)-mean(min_age))
  z<-average_lifespan$round_max_age
  plot(average_lifespan,xlab="Time \n (0 is Most Recent)", ylab="Mean Lifespan",col="dodgerblue",lwd=2,
       type="l",xaxt="n",main="Average Lifespan \n through Time",xlim=c(max(z),min(z)))
  year.labels<-seq(floor(min(z)),ceiling(max(z)),by=y)
  axis(1,at=seq(floor(min(z)),ceiling(max(z)),by=y),labels=year.labels)
}


###Diversity Through Time
diversity_time<-function(x,y){
  x$ts<-round(x$ts,digits=0)
  x$te<-round(x$te,digits=0)
  ts.counts<-x %>% group_by(ts) %>% summarize(total_ts=n_distinct(species)) %>% mutate(cumsum = cumsum(total_ts)) %>% mutate(ID=ts)
  plot(ts.counts$ts,ts.counts$cumsum,type="l",xlab="Time \n (0 is Most Recent)",col="dodgerblue",lwd=2,
       ylab="Cumulative Diversity",xaxt="n",main="Cumulative Diversity \n Through Time")
  year.labels<-seq(floor(min(ts.counts$ts)),ceiling(max(ts.counts$ts)),by=y)
  axis(1,at=seq(floor(min(ts.counts$ts)),ceiling(max(ts.counts$ts)),by=y),labels=rev(year.labels))
}

diversity_time_occurrence<-function(x,y){
  x$max_age<-round(x$max_age,digits=0)
  x$min_age<-round(x$min_age,digits=0)
  max_age.counts<-x %>% group_by(max_age) %>% summarize(total_max_age=n_distinct(Species)) %>% 
  mutate(cumsum = cumsum(total_max_age)) %>% mutate(ID=max_age)
  plot(max_age.counts$max_age,max_age.counts$cumsum,type="l",xlab="Time \n (0 is Most Recent)",col="dodgerblue",lwd=2,
       ylab="Cumulative Diversity",xaxt="n",main="Cumulative Diversity \n Through Time")
  year.labels<-seq(floor(min(max_age.counts$max_age)),ceiling(max(max_age.counts$max_age)),by=y)
  axis(1,at=seq(floor(min(max_age.counts$max_age)),ceiling(max(max_age.counts$max_age)),by=y),labels=rev(year.labels))
}

summary_plots<-function(x,y){
  par(mfrow=c(2,2))
  plot_LTT(x,y)
  diversity_time(x,y)
  average_lifespan(x,y)
  lifespan_hist(x,y)
  par(mfrow=c(1,1))
}

summary_plots_occurrence<-function(x,y){
  par(mfrow=c(2,2))
  plot_LTT_occurrence(x,y)
  diversity_time_occurrence(x,y)
  average_lifespan_occurrence(x,y)
  lifespan_hist_occurrence(x,y)
  par(mfrow=c(1,1))
}
  
 plot_RTT <- function (age,hpd_M,hpd_m,mean_m,color,RM){
  N=100
  rem = 1
  for (i in 1:(N-1)){
    trans=1/N#-0.0045
    polygon(c(age, rev(age)), c(hpd_M-((hpd_M-mean_m)*rem), rev(hpd_m+((mean_m-hpd_m)*rem))), col = alpha(color,trans), border = NA)
    rem = rem-(1/N)		
    #print(c(rem,trans))
  }
  lines(rev(age), rev(mean_m), col = color, lwd=3)  
} 
   
 shift.points.func<-function(x,y){
  shift.data<-hist(x,breaks=(max(round(x))-min(round(x))),plot=FALSE)
  shift.data<-data.frame(breaks=shift.data$mids,counts=shift.data$counts)
  summary<-data.frame(max_age=(y-shift.data$breaks),
                      min_age=(y-shift.data$breaks),
                      counts=shift.data$counts,
                      proportion=(shift.data$counts/sum(shift.data$counts)))
  summary<-arrange(summary,desc(counts))
  summary_top_3<-head(summary,3)
  print(summary_top_3)
}

