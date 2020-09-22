

p<-ggplot(plot_dat,aes(x=State,y=Value,fill=Tissue_state),show.legend=FALSE)+
     geom_boxplot(outlier.shape=18)+
     stat_compare_means(aes(group=Tissue),label="p.format",hide.ns=TRUE,label.y=5.5)+
     scale_x_discrete(labels=lbs)+theme(legend.position="none")+xlab("")+ylab("Fraction of Overlap")+scale_fill_manual(values=plotcol)




#Rotate axis labels
#+theme(axis.text.x=element_text(angle=45))


