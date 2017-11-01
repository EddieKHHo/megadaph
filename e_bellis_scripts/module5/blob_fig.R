library(ggplot2)
library(viridis)

daphFA <- read.table("/Users/weissem/Desktop/Daphnia/FA_SC.blob.txt", header=F, fill=T)
colnames(daphFA) <- c("Contig","GC","Length","Cov","HitID","HitLength","Eval","PID","HSPCov")
daph2 <- daphFA[order(daphFA$HitLength),]
daph2 <- subset(daph2, Length>5000)

pdf(file="/Users/weissem/Desktop/Daphnia/MS_draft/blob_hist.pdf", width=3.4, height=3, pointsize=10, colormodel='cmyk')
hist(daph2$Cov,xlim=c(0,100),breaks=700, prob=TRUE, xlab="Coverage", main="FA Contigs >5kb", col="grey70")
d <- density(daph2$Cov)
min1 <- optimize(approxfun(d$x,d$y),interval=c(20,60))$minimum  #37
min2 <- optimize(approxfun(d$x,d$y),interval=c(60,80))$minimum  #78
optimize(approxfun(d$x,d$y),interval=c(0,100), maximum=TRUE)$maximum ##55
lines(d, col="red")
abline(v=min1, col="black", lty=2)
abline(v=min2, col="black", lty=2)
dev.off()

pFA <- ggplot(daph2, aes(x=GC, y=Cov, size=Length, col=HitLength))+geom_point(alpha=0.1)+ geom_hline(yintercept=min1, lty=2) + geom_hline(yintercept=min2, lty=2)+ theme_classic() + xlab("% GC") + ylab("Coverage")+ labs(col='HSP Length', size="Contig Length") + ylim(c(0,150))+ scale_color_viridis(option="plasma")+ annotate("text", label="FA", x=30, y=145, size = 6)+ scale_size(guide = 'none') + xlim(c(25,60))+ theme(legend.position="none")+xlab("")+ylab("")

####FB_SC
daphFB <- read.table("/Users/weissem/Desktop/Daphnia/FB_SC.blob.txt", header=F, fill=T)
colnames(daphFB) <- c("Contig","GC","Length","Cov","HitID","HitLength","Eval","PID","HSPCov")
daph2 <- daphFB[order(daphFB$HitLength),]
daph2 <- subset(daph2, Length>5000)

hist(daph2$Cov,xlim=c(0,100),breaks=700, prob=TRUE, xlab="Coverage", main="", col="grey70")
d <- density(daph2$Cov)
min1 <- optimize(approxfun(d$x,d$y),interval=c(20,60))$minimum  #50
min2 <- optimize(approxfun(d$x,d$y),interval=c(60,100))$minimum  #94
optimize(approxfun(d$x,d$y),interval=c(0,100), maximum=TRUE)$maximum ##70
lines(d, col="red")
abline(v=min1, col="black", lty=2)
abline(v=min2, col="black", lty=2)


pFB <- ggplot(daph2, aes(x=GC, y=Cov, size=Length, col=HitLength))+geom_point(alpha=0.1)+ geom_hline(yintercept=min1, lty=2) + geom_hline(yintercept=min2, lty=2)+ theme_classic() + xlab("% GC") + ylab("Coverage")+ labs(col='HSP Length', size="Contig Length") + ylim(c(0,150))+ scale_color_viridis(option="plasma")+ annotate("text", label="FB", x=30, y=145, size = 6)+ theme(legend.position="none")+xlim(c(25,60))+xlab("")+ylab("")

####FC_SC
daphFC <- read.table("/Users/weissem/Desktop/Daphnia/FC_SC.blob.txt", header=F, fill=T)
colnames(daphFC) <- c("Contig","GC","Length","Cov","HitID","HitLength","Eval","PID","HSPCov")
daph2 <- daphFC[order(daphFC$HitLength),]
daph2 <- subset(daph2, Length>5000)

hist(daph2$Cov,xlim=c(0,100),breaks=200, prob=TRUE, xlab="Coverage", main="", col="grey70")
d <- density(daph2$Cov)
min1 <- optimize(approxfun(d$x,d$y),interval=c(20,50))$minimum  #32
min2 <- optimize(approxfun(d$x,d$y),interval=c(50,100))$minimum  #70
optimize(approxfun(d$x,d$y),interval=c(20,100), maximum=TRUE)$maximum ##70
lines(d, col="red")
abline(v=min1, col="black", lty=2)
abline(v=min2, col="black", lty=2)


pFC <- ggplot(daph2, aes(x=GC, y=Cov, size=Length, col=HitLength))+geom_point(alpha=0.1)+ geom_hline(yintercept=min1, lty=2) + geom_hline(yintercept=min2, lty=2)+ theme_classic() + xlab("% GC") + ylab("Coverage")+ labs(col='HSP Length', size="Contig Length") + ylim(c(0,150))+ scale_color_viridis(option="plasma")+ annotate("text", label="FC", x=30, y=145, size = 6)+xlim(c(25,60))+theme(legend.position="none")+xlab("")+ylab("")

####IA_SC
daphIA <- read.table("/Users/weissem/Desktop/Daphnia/IA_SC.blob.txt", header=F, fill=T)
colnames(daphIA) <- c("Contig","GC","Length","Cov","HitID","HitLength","Eval","PID","HSPCov")
daph2 <- daphIA[order(daphIA$HitLength),]
daph2 <- subset(daph2, Length>5000)

hist(daph2$Cov,xlim=c(0,100),breaks=200, prob=TRUE, xlab="Coverage", main="", col="grey70")
d <- density(daph2$Cov)
min1 <- optimize(approxfun(d$x,d$y),interval=c(20,50))$minimum
min2 <- optimize(approxfun(d$x,d$y),interval=c(50,100))$minimum
optimize(approxfun(d$x,d$y),interval=c(20,100), maximum=TRUE)$maximum
lines(d, col="red")
abline(v=min1, col="black", lty=2)
abline(v=min2, col="black", lty=2)

pIA <- ggplot(daph2, aes(x=GC, y=Cov, size=Length, col=HitLength))+geom_point(alpha=0.1)+ geom_hline(yintercept=min1, lty=2) + geom_hline(yintercept=min2, lty=2)+ theme_classic() + xlab("% GC") + ylab("Coverage")+ labs(col='HSP Length', size="Contig Length") + ylim(c(0,150))+ scale_color_viridis(option="plasma")+ annotate("text", label="IA", x=30, y=145, size = 6)+xlim(c(25,60))+theme(legend.position="none")+xlab("")+ylab("")

####IB_SC
daphIB <- read.table("/Users/weissem/Desktop/Daphnia/IB_SC.blob.txt", header=F, fill=T)
colnames(daphIB) <- c("Contig","GC","Length","Cov","HitID","HitLength","Eval","PID","HSPCov")
daph2 <- daphIB[order(daphIB$HitLength),]
daph2 <- subset(daph2, Length>10000)

hist(daph2$Cov,xlim=c(0,100),breaks=100, prob=TRUE, xlab="Coverage", main="", col="grey70")
d <- density(daph2$Cov)
min1 <- optimize(approxfun(d$x,d$y),interval=c(30,50))$minimum
min2 <- optimize(approxfun(d$x,d$y),interval=c(50,80))$minimum
optimize(approxfun(d$x,d$y),interval=c(20,100), maximum=TRUE)$maximum
lines(d, col="red")
abline(v=min1, col="black", lty=2)
abline(v=min2, col="black", lty=2)

pIB <- ggplot(daph2, aes(x=GC, y=Cov, size=Length, col=HitLength))+geom_point(alpha=0.1)+ geom_hline(yintercept=min1, lty=2) + geom_hline(yintercept=min2, lty=2)+ theme_classic() + xlab("% GC") + ylab("Coverage")+ labs(col='HSP Length', size="Contig Length") + ylim(c(0,150))+ scale_color_viridis(option="plasma")+ annotate("text", label="IB", x=30, y=145, size = 6)+xlim(c(25,60))+theme(legend.position="none")+xlab("%GC")+ylab("")

####IC_SC
daphIC <- read.table("/Users/weissem/Desktop/Daphnia/IC_SC.blob.txt", header=F, fill=T)
colnames(daphIC) <- c("Contig","GC","Length","Cov","HitID","HitLength","Eval","PID","HSPCov")
daph2 <- daphIC[order(daphIC$HitLength),]
daph2 <- subset(daph2, Length>5000)

hist(daph2$Cov,xlim=c(0,100),breaks=150, prob=TRUE, xlab="Coverage", main="", col="grey70")
d <- density(daph2$Cov)
min1 <- optimize(approxfun(d$x,d$y),interval=c(30,50))$minimum
min2 <- optimize(approxfun(d$x,d$y),interval=c(50,80))$minimum
optimize(approxfun(d$x,d$y),interval=c(20,100), maximum=TRUE)$maximum
lines(d, col="red")
abline(v=min1, col="black", lty=2)
abline(v=min2, col="black", lty=2)

pIC <- ggplot(daph2, aes(x=GC, y=Cov, size=Length, col=HitLength))+geom_point(alpha=0.1)+ geom_hline(yintercept=min1, lty=2) + geom_hline(yintercept=min2, lty=2)+ theme_classic() + xlab("% GC") + ylab("Coverage")+ labs(col='HSP Length', size="Contig Length") + ylim(c(0,150))+ scale_color_viridis(option="plasma")+ annotate("text", label="IC", x=30, y=145, size = 6)+xlim(c(25,60))+theme(legend.position="none")+xlab("")+ylab("")

####GA_SC
daphGA <- read.table("/Users/weissem/Desktop/Daphnia/GA_SC.blob.txt", header=F, fill=T)
colnames(daphGA) <- c("Contig","GC","Length","Cov","HitID","HitLength","Eval","PID","HSPCov")
daph2 <- daphGA[order(daphGA$HitLength),]
daph2 <- subset(daph2, Length>5000)

hist(daph2$Cov,xlim=c(0,100),breaks=300, prob=TRUE, xlab="Coverage", main="", col="grey70")
d <- density(daph2$Cov)
min1 <- optimize(approxfun(d$x,d$y),interval=c(30,50))$minimum
min2 <- optimize(approxfun(d$x,d$y),interval=c(50,70))$minimum
optimize(approxfun(d$x,d$y),interval=c(20,100), maximum=TRUE)$maximum
lines(d, col="red")
abline(v=min1, col="black", lty=2)
abline(v=min2, col="black", lty=2)

pGA <- ggplot(daph2, aes(x=GC, y=Cov, size=Length, col=HitLength))+geom_point(alpha=0.1)+ geom_hline(yintercept=min1, lty=2) + geom_hline(yintercept=min2, lty=2)+ theme_classic() + xlab("% GC") + ylab("Coverage")+ labs(col='HSP Length', size="Contig Length") + ylim(c(0,150))+ scale_color_viridis(option="plasma")+ annotate("text", label="GA", x=30, y=145, size = 6)+xlim(c(25,60))+theme(legend.position="none")+xlab("")+ylab("Coverage")

####GB_SC
daphGB <- read.table("/Users/weissem/Desktop/Daphnia/GB_SC.blob.txt", header=F, fill=T)
colnames(daphGB) <- c("Contig","GC","Length","Cov","HitID","HitLength","Eval","PID","HSPCov")
daph2 <- daphGB[order(daphGB$HitLength),]
daph2 <- subset(daph2, Length>5000)

hist(daph2$Cov,xlim=c(0,100),breaks=300, prob=TRUE, xlab="Coverage", main="", col="grey70")
d <- density(daph2$Cov)
min1 <- optimize(approxfun(d$x,d$y),interval=c(30,50))$minimum
min2 <- optimize(approxfun(d$x,d$y),interval=c(50,100))$minimum
optimize(approxfun(d$x,d$y),interval=c(20,100), maximum=TRUE)$maximum
lines(d, col="red")
abline(v=min1, col="black", lty=2)
abline(v=min2, col="black", lty=2)

pGB <- ggplot(daph2, aes(x=GC, y=Cov, size=Length, col=HitLength))+geom_point(alpha=0.1)+ geom_hline(yintercept=min1, lty=2) + geom_hline(yintercept=min2, lty=2)+ theme_classic() + xlab("% GC") + ylab("Coverage")+ labs(col='HSP Length', size="Contig Length") + ylim(c(0,150))+ scale_color_viridis(option="plasma")+ annotate("text", label="GB", x=30, y=145, size = 6)+xlim(c(25,60))+theme(legend.position="none")+xlab("")+ylab("")

####GC_SC
daphGC <- read.table("/Users/weissem/Desktop/Daphnia/GC_SC.blob.txt", header=F, fill=T)
colnames(daphGC) <- c("Contig","GC","Length","Cov","HitID","HitLength","Eval","PID","HSPCov")
daph2 <- daphGC[order(daphGC$HitLength),]
daph2 <- subset(daph2, Length>5000)

hist(daph2$Cov,xlim=c(0,100),breaks=300, prob=TRUE, xlab="Coverage", main="", col="grey70")
d <- density(daph2$Cov)
min1 <- optimize(approxfun(d$x,d$y),interval=c(30,50))$minimum
min2 <- optimize(approxfun(d$x,d$y),interval=c(50,80))$minimum
optimize(approxfun(d$x,d$y),interval=c(20,100), maximum=TRUE)$maximum
lines(d, col="red")
abline(v=min1, col="black", lty=2)
abline(v=min2, col="black", lty=2)

pGC <- ggplot(daph2, aes(x=GC, y=Cov, size=Length, col=HitLength))+geom_point(alpha=0.1)+ geom_hline(yintercept=min1, lty=2) + geom_hline(yintercept=min2, lty=2)+ theme_classic() + xlab("% GC") + ylab("Coverage")+ labs(col='HSP Length', size="Contig Length") + ylim(c(0,150))+ scale_color_viridis(option="plasma")+ annotate("text", label="GC", x=30, y=145, size = 6)+xlim(c(25,60))+theme(legend.position="none")+xlab("")+ylab("")


#+ scale_size(guide = 'none')+xlim(c(25,60))



####composite figure
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

pdf(file="/Users/weissem/Desktop/Daphnia/MS_draft/blobs.pdf", width=6, height=6.5, pointsize=10, colormodel='cmyk')
multiplot(pFA,pGA,pIA, pFB, pGB, pIB, pFC, pGC,pIC, cols=3)
dev.off()

pdf(file="/Users/weissem/Desktop/Daphnia/MS_draft/blob_legends.pdf", width=3, height=3, pointsize=10, colormodel='cmyk')
ggplot(daph2, aes(x=GC, y=Cov, size=Length, col=HitLength))+geom_point(alpha=0.1)+ geom_hline(yintercept=min1, lty=2) + geom_hline(yintercept=min2, lty=2)+ theme_classic() + xlab("% GC") + ylab("Coverage")+ labs(col='HSP Length', size="Contig Length") + ylim(c(0,150))+ scale_color_viridis(option="plasma")+ annotate("text", label="GC", x=30, y=145, size = 6)+xlim(c(25,60))+xlab("")+ylab("")
dev.off()
