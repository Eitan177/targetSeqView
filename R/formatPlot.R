formatPlot <-
    function(dfplot,mismatchcolor='red',title='',covMaxDisplay=Inf,flipLeftandRight=FALSE,cropGaps=TRUE,plotJunctions=NULL,nolegend=FALSE){
        ## TO SATISFY R CMD CHECK
        posvec2 <- afvec2 <- ypostail <- NULL #l <-  NULL
        assign("l",c('readEnd','C','A','G','T','Indel','Match','Mismatch',
                     paste('Ref',c('C','A','G','T'))),envir=.GlobalEnv)
                                        #l=c('E','C','A','G','T','-','M')

        if(max(dfplot$ypostail) > covMaxDisplay){
            dfm <- list()
            dfl=split(dfplot,dfplot$ref)
            for(ii in seq_along(dfl)){
                dfm[[ii]] <- dfl[[ii]]
                allreadpos <- dfm[[ii]]$posvec2[dfm[[ii]]$ypostail != -1 & dfm[[ii]]$ypostail*ifelse(covMaxDisplay<0,-1,1) <= covMaxDisplay]
                dfm[[ii]] <- dfm[[ii]][-which(dfm[[ii]]$ypostail*ifelse(covMaxDisplay<0,-1,1) > covMaxDisplay | (dfm[[ii]]$ypostail == -1 &
                                      (dfm[[ii]]$posvec2 > max(allreadpos) | dfm[[ii]]$posvec2 < min(allreadpos)))),]
            }
            dfplot <- do.call(rbind,dfm)
        }

        if(cropGaps == TRUE){
            tt=table(dfplot$posvec2,dfplot$afvec2)
            locsOnlyGaps=apply(tt,1,function(x){x[which(names(x)=='-')]==sum(x)-1})
            if(any(locsOnlyGaps)){
                gapPos <- as.numeric(names(locsOnlyGaps))
                dfm <- list()
                dfl=split(dfplot,dfplot$ref)
                for(ii in seq_along(dfl)){
                dfm[[ii]] <- dfl[[ii]]
                matchpos <- dfm[[ii]]$posvec2[dfm[[ii]]$afvec2=='M']
                torem <- which(dfm[[ii]]$posvec2 %in% gapPos & (dfm[[ii]]$posvec2 > max(matchpos)+20 | dfm[[ii]]$posvec2 < min(matchpos)-20))
                if(length(torem)>0){
                    dfm[[ii]] <- dfm[[ii]][-torem,]
                }
            }
                dfplot <- do.call(rbind,dfm)
            }
        }

        #dfplot$posvec2 <- dfplot$posvec2/1e6
        if(mismatchcolor != ''){
            dfplot$afvec2<-as.character(dfplot$afvec2)
            dfplot$afvec2[dfplot$afvec2 %in% c("A","C","G","T") & dfplot$ypostail > 0] <- "Mismatch"
            dfplot$afvec2[dfplot$afvec2 == "M"] <- "Match"
            dfplot$afvec2[dfplot$afvec2 == "E"] <- "readEnd"
            dfplot$afvec2[dfplot$afvec2 == "-"] <- "Indel"
            dfplot$afvec2[dfplot$afvec2 %in% c("A")] <- "Ref A"
            dfplot$afvec2[dfplot$afvec2 %in% c("G")] <- "Ref C"
            dfplot$afvec2[dfplot$afvec2 %in% c("C")] <- "Ref G"
            dfplot$afvec2[dfplot$afvec2 %in% c("T")] <- "Ref T"
        }
        if(flipLeftandRight){
            #dfplot$ref=factor(dfplot$ref,levels=rev(unique(dfplot$ref)))
            level1<-levels(factor(dfplot$ref))[2]
            sA<-split(dfplot,dfplot$ypostail)
            mins<-sapply(sA,function(x){min(x$posvec2[as.character(x$ref)==level1])})
            sB<-sA[order(mins)][-1]
            reord<-mapply(function(x,y){rep(y,nrow(x))},sB,as.list(1:length(sB)))
            dfplot<-do.call(rbind,sB)
            dfplot$ypostail<-unlist(reord)
            dfplot<-rbind(dfplot,sA[[1]])
            dfplot$ref=factor(dfplot$ref,levels=rev(unique(dfplot$ref)))
        }

        levelsUse<-levels(factor(dfplot$ref))
        innerFirst<-max(dfplot$posvec2[dfplot$ref == levelsUse[1]])
        outerFirst<-min(dfplot$posvec2[dfplot$ref == levelsUse[1]])
        outerSecond<-max(dfplot$posvec2[dfplot$ref == levelsUse[2]])
        innerSecond<-min(dfplot$posvec2[dfplot$ref == levelsUse[2]])
        if(length(levelsUse)==1){
        dfplot$ref2[as.character(dfplot$ref) == levelsUse[1]]<-paste('Contiguous Fragment Alignment Side',
                                 which.min(abs(as.numeric(plotJunctions)-as.numeric(levelsUse))))
    }else{
        dfplot$ref2[as.character(dfplot$ref) == levelsUse[1]]<-paste('Side 1 Chr',dfplot$chr[as.character(dfplot$ref) == levelsUse[1]],
                                 ':',outerFirst,'-',innerFirst,sep='')
        dfplot$ref2[as.character(dfplot$ref) == levelsUse[2]]<-paste('Side 2 Chr',dfplot$chr[as.character(dfplot$ref) == levelsUse[2]],
                                 ':',innerSecond,'-',outerSecond,sep='')
    }
        dftail<-dftailB<-dftailC<-dftailD<-dfplot
        dfplot$ypostail<-dftail$ypostail*4
        dftailB$ypostail<- (dftail$ypostail*4)-1
        dftailC$ypostail<- (dftail$ypostail*4)-2
        dftailD$ypostail<- (dftail$ypostail*4)-3
        dfplot$hei<-rep(0,nrow(dftailB))
        dftailB$hei<-rep(1,nrow(dftailB))
        dftailC$hei<-rep(1,nrow(dftailB))
        dftailD$hei<-rep(1,nrow(dftailB))

        splitdftailC<-split(dftailC[dftailC$ypostail>0,],paste(dftailC$ypostail[dftailC$ypostail>0],dftailC$ref[dftailC$ypostail>0],sep='.'))
        if(length(levelsUse)==1){
            makedashed<-lapply(splitdftailC,function(x){
            diffs<-diff(x[,2]);
            if(max(diffs)>1){
            xs<-seq(x[which.max(diffs),2],x[which.max(diffs)+1,2],4);xs<-c(xs,xs+1);

            toret<-data.frame(afvec2=rep('readEnd',length(xs)),posvec2=xs,ref=rep(x$ref[1],length(xs)),chr=rep(x$chr[1],length(xs)),ypostail=rep(x$ypostail[1],
length(xs)),refname=rep(x$refname[1],length(xs)),hei=rep(.5,length(xs)),ref2=rep(x$ref2[1],length(xs)))}})
        }else{

        makedashed<-lapply(splitdftailC,function(x){
            if(as.character(x$ref[1])==levelsUse[1]){xs<-seq(max(x$posvec2),innerFirst,4);xs<-c(xs,xs+1)
                                   }else{xs<-seq(innerSecond,min(x$posvec2),4);xs<-c(xs,xs+1)
                                     };
            toret<-data.frame(afvec2=rep('readEnd',length(xs)),posvec2=xs,ref=rep(x$ref[1],length(xs)),chr=rep(x$chr[1],length(xs)),ypostail=rep(x$ypostail[1],
length(xs)),refname=rep(x$refname[1],length(xs)),hei=rep(.5,length(xs)),ref2=rep(x$ref2[1],length(xs)))})
    }
        dashedLineRows<-do.call(rbind,makedashed)
        dftailE<-rbind(dfplot,dftailB,dftailC,dftailD,dashedLineRows)
        dftailE$hei<-dftailE$hei/4
        dftailE$ypostail<-dftailE$ypostail/4
        dtailA  <-  ggplot(dftailE,aes(x=posvec2,y=ypostail,fill=factor(afvec2,levels=get("l")),height=hei)) + geom_tile()+facet_grid(.~ref2,scales="free")#,space="free");
        xlims <- c(min(dfplot$posvec2), max(dfplot$posvec2))
        base_colors=c('red','yellow','blue','green')
        trimmedA <-  dtailA +scale_x_continuous("Position" )+
            scale_fill_manual( values=c('black',base_colors,'#99FFFF','grey',mismatchcolor,base_colors)[!is.na(match(factor(l),levels(factor(dfplot$afvec2))))])+
         scale_y_continuous("Coverage")
        trimA <- trimmedA+labs(fill="",title=title)+
            scale_linetype(guide = "none")+theme(panel.background = element_blank(),text = element_text(size=16),axis.text.x = element_text(angle=30, vjust=1),
                           legend.position=ifelse(nolegend,"none","top"))

        return(list(trimA,levelsUse))
    }
