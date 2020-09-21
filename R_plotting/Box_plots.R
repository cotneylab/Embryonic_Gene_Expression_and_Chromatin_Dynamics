p<-ggplot(plot_dat,aes(x=State,y=Value,fill=Tissue_state),show.legend=FALSE)+
     geom_boxplot(outlier.shape=18)+
     stat_compare_means(aes(group=Tissue),label="p.format",hide.ns=TRUE,label.y=5.5)+
     scale_x_discrete(labels=lbs)+theme(legend.position="none")+xlab("")+ylab("Fraction of Overlap")+scale_fill_manual(values=plotcol)


lbs<-c("1_TssA","2_PromU","3_PromD1","4_PromD2","5_Tx5'","6_Tx",
               "7_Tx3'","8_TxWk","9_TxReg","10_TxEnh5'","11_TxEnh3'","12_TxEnhW",
               "13_EnhA1","14_EnhA2","15_EnhAF","16_EnhW1","17_EnhW2","18_EnhAc",
               "19_Dnase","20_ZNF_Rpts","21_Het","22_PromP","23_PromBiv","24_ReprPC",
               "25_Quies")
