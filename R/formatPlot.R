formatPlot <-
    function(dftailA,title='Alignment Plot'){
        ## TO SATISFY R CMD CHECK
        posvec2 <- afvec2 <- ypostail <- NULL #l <-  NULL


        assign("l",c('E','C','A','G','T','-','M'),envir=.GlobalEnv)
                                        #l=c('E','C','A','G','T','-','M')

        dftailA$posvec2 <- dftailA$posvec2/1e6
        dtailA  <-  ggplot(dftailA,aes(x=posvec2,y=ypostail,fill=factor(afvec2,levels=get("l")),height=.85)) + geom_tile()+facet_grid(.~ref,scales="free");#stat_bin2d(binwidth=c(1,1))

        xlims <- c(min(dftailA$posvec2), max(dftailA$posvec2))
        trimmedA <-  dtailA +scale_x_continuous("Position (Mb)" )+scale_fill_manual( values=c('black','red','yellow','blue','green', '#99FFFF','grey')[!is.na(match(factor(l),levels(dftailA$afvec2)))])+ scale_y_continuous("Coverage")


        trimA <- trimmedA+labs(fill="",title=title)+
            scale_linetype(legend = FALSE)+opts(panel.background = theme_blank())

        return(trimA)
    }
