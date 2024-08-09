#npstat_plot_windows.R
#options: outgroup: yes/no annot_file: yes/no
#ncols <- 18
# R --vanilla < ./npstat_plot_windows.R [infile] [outfile]

if(!require("ggplot2")) install.packages("ggplot2")
library("ggplot2")

args = commandArgs(trailingOnly=TRUE)

if(length(args)<2) {
  #  show(sprintf("Please specify arguments: \n R --vanilla < ./npstat_plot_windows.R [infile] [outfile]"))
  #  exit()
  args <- array(dim=c(2))
  args[1] <- "Pool_seq2.mel_mpileup_2L.txt.stats.ex1.txt"#"inputfile.stats.txt"
  args[2] <- "Pool_seq2.mel_mpileup_2L.txt.stats.ex1.pdf"#"outputfile.pdf" 
  #  args[3] <- "20"
}
#interv <- as.numeric(args[3])
infile <- read.table(file=args[1],header=T)
#outgroup
outg <- 0
if(sum(infile$length_outgroup))  outg <- 1
#annot_file
gff <- 0
if(sum(!is.na(infile$nonsyn_pol))>0) gff <- 1

ncols <- c(2,4:8)
if(outg==0 && gff==1) ncols <- c(ncols,c(14,15))
if(outg) {
  ncols <- c(2:13)
  if(gff) ncols <- c(ncols,c(14:18))
}

pdf(file=args[2],width=10,height=5)
par(mfrow=c(1,1))

#plot sliding windows
#  len.win <- length(infile[,1])
#  x<-seq(from=1,to=len.win,by=interv)
for(nc in ncols) {
  if(nc<18) {
    plot(infile[,nc],type="p",pch=20,main=sprintf("Sliding Windows: \n%s",colnames(infile)[nc]),xlab="Scaffold",ylab=sprintf("%s",colnames(infile)[nc]),ylim=c(min(0,infile[,nc],na.rm=T),max(infile[,nc],na.rm=T)))
  }
  if(nc==18) {
    plot(infile[,nc],type="p",pch=20,main=sprintf("Sliding Windows: \n%s",colnames(infile)[nc]),xlab="Scaffold",ylab=sprintf("%s",colnames(infile)[nc]))
  }
  abline(h=0,col="grey")
  #  mi <- NULL
  #  for(y1 in x) {
  #    y2 <- seq(y1,y1+interv-1,1)
  #    mi <- c(mi,mean(infile[y2,nc],na.rm=T))
  #  }
  #  lines(x=x,y=mi,col="blue",lwd=1)
}
#plot densities
plot_list = list(); i<-1
for(nc in ncols) {
  if(sum(!is.na(infile[,nc]))>10) {
    x <- data.frame(Values=infile[!(is.na(infile[,nc]) | is.infinite(infile[,nc])),nc],Stat=colnames(infile)[nc])
    p <- ggplot(data.frame(x), aes(x=Stat,y=Values, color=Stat)) + geom_violin(trim=FALSE) + labs(title=colnames(infile)[nc])
    p <- p + geom_jitter(shape=16, position=position_jitter(0.2)) + geom_boxplot(width=0.1,col="black") + theme_classic()
    plot_list[[i]] = p; i <- i + 1

    # plot(density(infile[,nc],na.rm=T),main=sprintf("Density Distribution: %s \n mean = %.3e",colnames(infile)[nc],mean(infile[,nc],na.rm=T)),xlab="Values")
  }
}
if(!outg) {
  #sliding windows
  plot(infile$Pi,type="p",pch=20,main=sprintf("Sliding Windows: \nPi vs theta Watt"),xlab="Scaffold",ylab=sprintf("Variability"),ylim=c(0,max(infile$Pi,infile$Watterson,na.rm=T)))
  lines(infile$Watterson,type="p",col="red",pch=20)
  legend("topleft",legend=c("Pi","Watt"),col=c("black","red"),pch=c(20,20))
  abline(h=0,col="grey")
  
  #plot densities
  if(sum(!is.na(infile$Pi))>10) {
    x1 <- data.frame(Values=infile$Pi,Stat="Pi")
    x2 <- data.frame(Values=infile$Watterson,Stat="Watterson")
    x <- rbind(x1,x2)
    p <- ggplot(data.frame(x), aes(x=Stat,y=Values, color=Stat)) + geom_violin(trim=FALSE) + labs(title="Pi vs theta Watt") + scale_x_discrete(limits=c("Pi", "Watterson"))
    p <- p + geom_jitter(shape=16, position=position_jitter(0.2))  + geom_boxplot(width=0.1,col="black") + theme_classic()
    plot_list[[i]] = p; i <- i + 1
    
    # plot(density(infile$Watterson,na.rm=T),col="red",main=sprintf("Density Distribution:\n Pi vs theta Watt"),xlim=c(0,max(infile$Pi,infile$Watterson,na.rm=T)),xlab="Values")
    # lines(density(infile$Pi,na.rm=T))
    # legend("topright",legend=c("Pi","Watt"),col=c("black","red"),lty=c(1,1))
    # abline(v=mean(infile$Pi,na.rm=T),col="black")
    # abline(v=mean(infile$Watterson,na.rm=T),col="red")
  }
}
if(outg) {
  #sliding windows
  plot(infile$Pi,type="p",pch=20,main=sprintf("Sliding Windows: \nPi, thetaWatt and thetaH"),xlab="Scaffold",ylab=sprintf("Variability"),ylim=c(0,max(infile$Pi,infile$Watterson,infile$Pi-infile$unnorm_FayWu_H,na.rm=T)))
  lines(infile$Watterson,type="p",col="red",pch=20)
  lines(infile$Pi-infile$unnorm_FayWu_H,type="p",col="blue",pch=20)
  legend("topleft",legend=c("Pi","Watt","thetaH"),col=c("black","red","blue"),pch=c(20,20,20))
  abline(h=0,col="grey")
  
  plot(infile$Pi,type="p",pch=20,main=sprintf("Sliding Windows:\n Pi vs div"),xlab="Scaffold",ylab=sprintf("Pi vs div"),ylim=c(0,max(infile$Pi,infile$div,na.rm=T)))
  lines(infile$div,type="p",col="red",pch=20)
  legend("topleft",legend=c("Pi","div"),col=c("black","red"),pch=c(20,20))
  abline(h=0,col="grey")
  
  plot(infile$Tajima_D,type="p",pch=20,main=sprintf("Sliding Windows:\n Tajima's D vs FayWu's H"),xlab="Scaffold",ylab=sprintf("Neutrality Tests"),ylim=c(min(infile$Tajima_D,infile$FayWu_H,na.rm=T),max(infile$Tajima_D,infile$FayWu_H,na.rm=T)))
  lines(infile$FayWu_H,type="p",col="red",pch=20)
  legend("topleft",legend=c("Tajima's D","FayWu's H"),col=c("black","red"),pch=c(20,20))
  abline(h=0,col="grey")
  
  #plot densities
  if(sum(!is.na(infile$Pi))>10) {
    x1 <- data.frame(Values=infile$Pi,Stat="Pi")
    x2 <- data.frame(Values=infile$Watterson,Stat="Watterson")
    x3 <- data.frame(Values=infile$Pi-infile$unnorm_FayWu_H,Stat="theta H")
    x <- rbind(x1,x2,x3)
    p <- ggplot(data.frame(x), aes(x=Stat,y=Values, color=Stat)) + geom_violin(trim=FALSE) + labs(title="Pi, theta Watt and theta H") + scale_x_discrete(limits=c("Pi", "Watterson", "theta H"))
    p <- p + geom_jitter(shape=16, position=position_jitter(0.2)) + geom_boxplot(width=0.1,col="black") + theme_classic()
    plot_list[[i]] = p; i <- i + 1
    
    # plot(density(infile$Pi,na.rm=T),main=sprintf("Density Distribution:\n Pi, theta Watt and theta FayWu"),xlim=c(min(infile$Pi,infile$Watterson,infile$unnorm_FayWu_H+infile$Pi,na.rm=T),max(infile$Pi,infile$Watterson,infile$unnorm_FayWu_H+infile$Pi,na.rm=T)),xlab="Values") 
    # lines(density(infile$Watterson,na.rm=T),col="red")
    # lines(density(infile$unnorm_FayWu_H+infile$Pi,na.rm=T),col="blue")
    # legend("topright",legend=c("Pi","Watt","FayWu"),col=c("black","red","blue"),lty=c(1,1,1))
    # abline(v=mean(infile$Pi,na.rm=T),col="black")
    # abline(v=mean(infile$Watterson,na.rm=T),col="red")
    # abline(v=mean(infile$unnorm_FayWu_H+infile$Pi,na.rm=T),col="blue")
    
    x1 <- data.frame(Values=infile$Pi,Stat="Pi")
    x2 <- data.frame(Values=infile$div,Stat="div")
    x <- rbind(x1,x2)
    p <- ggplot(data.frame(x), aes(x=Stat,y=Values, color=Stat)) + geom_violin(trim=FALSE) + labs(title="Pi vs Divergence") + scale_x_discrete(limits=c("Pi", "div"))
    p <- p + geom_jitter(shape=16, position=position_jitter(0.2)) + geom_boxplot(width=0.1,col="black") + theme_classic()
    plot_list[[i]] = p; i <- i + 1
    
    # plot(density(infile$Pi,na.rm=T),main=sprintf("Density Distribution:\n Pi vs Divergence"),xlim=c(0,max(infile$Pi,infile$div,na.rm=T)),xlab="Values")
    # lines(density(infile$div,na.rm=T),col="red")
    # legend("topright",legend=c("Pi","div"),col=c("black","red"),lty=c(1,1))
    # abline(v=mean(infile$Pi,na.rm=T),col="black")
    # abline(v=mean(infile$div,na.rm=T),col="red")
    
    x1 <- data.frame(Values=infile$Tajima_D,Stat="TajimaD")
    x2 <- data.frame(Values=infile$FayWu_H,Stat="FaywuH")
    x <- rbind(x1,x2)
    p <- ggplot(data.frame(x), aes(x=Stat,y=Values, color=Stat)) + geom_violin(trim=FALSE) + labs(title="Tajima's D vs Fay&Wu's H") + scale_x_discrete(limits=c("TajimaD", "FaywuH"))
    p <- p + geom_jitter(shape=16, position=position_jitter(0.2)) + geom_boxplot(width=0.1,col="black") + theme_classic()
    plot_list[[i]] = p; i <- i + 1
    
    # plot(density(infile$Tajima_D,na.rm=T),main=sprintf("Density Distribution:\n Tajima's D vs Fay&Wu's H"),xlim=c(min(infile$Tajima_D,infile$FayWu_H,na.rm=T),max(infile$Tajima_D,infile$FayWu_H,na.rm=T)),xlab="Values")
    # lines(density(infile$FayWu_H,na.rm=T),col="blue")
    # legend("topright",legend=c("Tajima's D","Fay & Wu's H"),col=c("black","blue"),lty=c(1,1))
    # abline(v=mean(infile$Tajima_D,na.rm=T),col="black")
    # abline(v=mean(infile$FayWu_H,na.rm=T),col="blue")
  }
}
if(outg && gff) {
  #14:nonsyn_pol	15:syn_pol	16:nonsyn_div	17:syn_div
  #sliding windows
  plot(infile$syn_pol,type="p",pch=20,main=sprintf("Sliding Windows:\n PolSyn vs PolNSyn"),xlab="Scaffold",ylab=sprintf("Pi(Syn) vs Pi(NSyn)"),ylim=c(0,max(infile$syn_pol,infile$nonsyn_pol,na.rm=T)))
  lines(infile$nonsyn_pol,type="p",col="red",pch=20)
  legend("topleft",legend=c("Pi_Syn","Pi_Nsyn"),col=c("black","red"),pch=c(20,20))
  abline(h=0,col="grey")
  
  plot(infile$syn_div,type="p",pch=20,main=sprintf("Sliding Windows:\n DivSyn vs DivNSyn"),xlab="Scaffold",ylab=sprintf("Div(Syn) vs Div(NSyn)"),ylim=c(0,max(infile$syn_div,infile$nonsyn_div,na.rm=T)))
  lines(infile$nonsyn_div,type="p",col="red",pch=20)
  legend("topleft",legend=c("Div_Syn","Div_Nsyn"),col=c("black","red"),pch=c(20,20))
  abline(h=0,col="grey")
  
  plot(infile$syn_pol,type="p",pch=20,main=sprintf("Sliding Windows:\n PolSyn vs DivSyn"),xlab="Scaffold",ylab=sprintf("Pol(Syn) vs Div(Syn)"),ylim=c(0,max(infile$syn_pol,infile$syn_div,na.rm=T)))
  lines(infile$syn_div,type="p",col="red",pch=20)
  legend("topleft",legend=c("Pol_Syn","Div_Syn"),col=c("black","red"),pch=c(20,20))
  abline(h=0,col="grey")
  
  plot(infile$nonsyn_pol,type="p",pch=20,main=sprintf("Sliding Windows:\n PolNsyn vs DivNsyn"),xlab="Scaffold",ylab=sprintf("Pol(Nsyn) vs Div(Nsyn)"),ylim=c(0,max(infile$nonsyn_pol,infile$nonsyn_div,na.rm=T)))
  lines(infile$nonsyn_div,type="p",col="red",pch=20)
  legend("topleft",legend=c("Pol_Nsyn","Div_Nsyn"),col=c("black","red"),pch=c(20,20))
  abline(h=0,col="grey")
  
  #plot densities
  if(sum(!is.na(infile$nonsyn_pol))>10) {
    x1 <- data.frame(Values=infile$nonsyn_pol,Stat="NSyn_Pol")
    x2 <- data.frame(Values=infile$syn_pol,Stat="Syn_Pol")
    x <- rbind(x1,x2)
    p <- ggplot(data.frame(x), aes(x=Stat,y=Values, color=Stat)) + geom_violin(trim=FALSE) + labs(title="Syn_Pol vs Nsyn_Pol")
    p <- p + geom_jitter(shape=16, position=position_jitter(0.2)) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") + theme_classic()
    plot_list[[i]] = p; i <- i + 1
    
    # plot(density(infile$nonsyn_pol,na.rm=T),main=sprintf("Density Distribution:\n PolSyn vs PolNSyn"),xlim=c(0,max(infile$syn_pol,infile$nonsyn_pol,na.rm=T)),xlab="Values")
    # lines(density(infile$syn_pol,na.rm=T),col="red")
    # legend("topright",legend=c("PolSyn","PolNSyn"),col=c("red","black"),lty=c(1,1))
    # abline(v=mean(infile$syn_pol,na.rm=T),col="black")
    # abline(v=mean(infile$nonsyn_pol,na.rm=T),col="red")
    
    x1 <- data.frame(Values=infile$nonsyn_div,Stat="NSyn_Div")
    x2 <- data.frame(Values=infile$syn_div,Stat="Syn_Div")
    x <- rbind(x1,x2)
    p <- ggplot(data.frame(x), aes(x=Stat,y=Values, color=Stat)) + geom_violin(trim=FALSE) + labs(title="Syn_Div vs Nsyn_Div")
    p <- p + geom_jitter(shape=16, position=position_jitter(0.2)) + geom_boxplot(width=0.1,col="black") + theme_classic()
    plot_list[[i]] = p; i <- i + 1
    
    # plot(density(infile$nonsyn_div,na.rm=T),main=sprintf("Density Distribution:\n DivSyn vs DivNSyn"),xlim=c(0,max(infile$syn_div,infile$nonsyn_div,na.rm=T)),xlab="Values")
    # lines(density(infile$syn_div,na.rm=T),col="red")
    # legend("topright",legend=c("DivSyn","DivNSyn"),col=c("red","black"),lty=c(1,1))
    # abline(v=mean(infile$syn_div,na.rm=T),col="black")
    # abline(v=mean(infile$nonsyn_div,na.rm=T),col="red")
    
    x1 <- data.frame(Values=infile$syn_pol,Stat="Syn_Pol")
    x2 <- data.frame(Values=infile$syn_div,Stat="Syn_Div")
    x <- rbind(x1,x2)
    p <- ggplot(data.frame(x), aes(x=Stat,y=Values, color=Stat)) + geom_violin(trim=FALSE) + labs(title="Syn_Pol vs Syn_Div")
    p <- p + geom_jitter(shape=16, position=position_jitter(0.2)) + geom_boxplot(width=0.1,col="black") + theme_classic()
    plot_list[[i]] = p; i <- i + 1
    
    # plot(density(infile$syn_pol,na.rm=T),main=sprintf("Density Distribution:\n DivSyn vs PolSyn"),xlim=c(0,max(infile$syn_div,infile$syn_pol,na.rm=T)),xlab="Values")
    # lines(density(infile$syn_div,na.rm=T),col="red")
    # legend("topright",legend=c("PolSyn","DivSyn"),col=c("black","red"),lty=c(1,1))
    # abline(v=mean(infile$syn_pol,na.rm=T),col="black")
    # abline(v=mean(infile$syn_div,na.rm=T),col="red")
    
    x1 <- data.frame(Values=infile$nonsyn_pol,Stat="Nsyn_Pol")
    x2 <- data.frame(Values=infile$nonsyn_div,Stat="Nsyn_Div")
    x <- rbind(x1,x2)
    p <- ggplot(data.frame(x), aes(x=Stat,y=Values, color=Stat)) + geom_violin(trim=FALSE) + labs(title="Nsyn_Pol vs Nsyn_Div")
    p <- p + geom_jitter(shape=16, position=position_jitter(0.2)) + geom_boxplot(width=0.1,col="black") + theme_classic()
    plot_list[[i]] = p; i <- i + 1
    
    # plot(density(infile$nonsyn_pol,na.rm=T),main=sprintf("Density Distribution:\n DivNSyn vs PolNSyn"),xlim=c(0,max(infile$nonsyn_div,infile$nonsyn_pol,na.rm=T)),xlab="Values")
    # lines(density(infile$nonsyn_div,na.rm=T),col="red")
    # legend("topright",legend=c("PolNSyn","DivNSyn"),col=c("black","red"),lty=c(1,1))
    # abline(v=mean(infile$nonsyn_pol,na.rm=T),col="black")
    # abline(v=mean(infile$nonsyn_div,na.rm=T),col="red")
  }
}
dev.off()
pdf(file=sprintf("%s.density.pdf",args[2]))
for(ii in c(1:(i-1))) {
  print(plot_list[[ii]])
}
dev.off()
