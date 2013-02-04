formatPlot <-
    function(dftailA,mismatchcolor='',title='Alignment Plot'){
        ## TO SATISFY R CMD CHECK
        posvec2 <- afvec2 <- ypostail <- NULL #l <-  NULL
        assign("l",c('E','C','A','G','T','-','M','A,C,G,T',
                     paste('Ref',c('C','A','G','T'))),envir=.GlobalEnv)
                                        #l=c('E','C','A','G','T','-','M')

        dftailA$posvec2 <- dftailA$posvec2/1e6
        if(mismatchcolor != ''){
            dftailA$afvec2<-as.character(dftailA$afvec2)
            dftailA$afvec2[dftailA$afvec2 %in% c("A","C","G","T") & dftailA$ypostail > 0] <- "A,C,G,T"
            dftailA$afvec2[dftailA$afvec2 %in% c("A")] <- "Ref A"
            dftailA$afvec2[dftailA$afvec2 %in% c("G")] <- "Ref C"
            dftailA$afvec2[dftailA$afvec2 %in% c("C")] <- "Ref G"
            dftailA$afvec2[dftailA$afvec2 %in% c("T")] <- "Ref T"
        }
        dtailA  <-  ggplot(dftailA,aes(x=posvec2,y=ypostail,fill=factor(afvec2,levels=get("l")),height=.85)) + geom_tile()+facet_grid(.~ref,scales="free");#stat_bin2d(binwidth=c(1,1))
        xlims <- c(min(dftailA$posvec2), max(dftailA$posvec2))
        base_colors=c('red','yellow','blue','green')
        trimmedA <-  dtailA +scale_x_continuous("Position (Mb)" )+scale_fill_manual( values=c('black',base_colors,'#99FFFF','grey',mismatchcolor,base_colors)[!is.na(match(factor(l),levels(factor(dftailA$afvec2))))])+ scale_y_continuous("Coverage")
        trimA <- trimmedA+labs(fill="",title=title)+
            scale_linetype(guide = "none")+theme(panel.background = element_blank())

        return(trimA)
    }
