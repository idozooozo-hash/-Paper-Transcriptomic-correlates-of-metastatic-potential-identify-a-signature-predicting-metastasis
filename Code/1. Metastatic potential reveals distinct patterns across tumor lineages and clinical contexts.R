########################################
## Directory
########################################
dir_data="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/data/"
dir_fig="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/figure/"



########################################
## Library
########################################
library(ComplexHeatmap)
library(colorRamp2)
library(stringr)
library(car)
library(dplyr)
library(ggplot2)



########################################
## Load Data
########################################
load(file=paste0(dir_data,"0. filtered sample information(DepMap).Rdata")) # variable name: sinfo
fsinfo=sinfo
load(file=paste0(dir_data,"0. all sample information(DepMap).Rdata")) # variable name: sinfo



########################################
## Metastasis Potential distribution by target organ
########################################
ttypes=colnames(fsinfo)[grep("meta to ",colnames(fsinfo))]

for(ttype in ttypes) {
    breaks=seq(floor(min(fsinfo[,ttype])) , ceiling(max(fsinfo[,ttype])), by=0.5)
    freqs=table(cut(fsinfo[,ttype], breaks=breaks, right=F))
    if(min(freqs)==0) {
        freqs=freqs[-1]
        breaks=breaks[-1]
    }
    if(max(freqs)==0) {
        freqs=freqs[-length(freqs)]
        breaks=breaks[-length(breaks)]
    }
    counts=as.numeric(freqs)

    width=1.7
    ymin=0
    ymax=250
    yseq=seq(ymin,ymax,50)
    ylim=range(ymin,ymax)
    tcols=setNames(c("#3d3aff","#ff3a3a") , c(-4,-2))

    tiff(filename=paste0(dir_fig,"Supplementary Figure 1(A)_",gsub("meta to ","",ttype),".tif"), width=23, height=20, units="cm", res=300)
    par(mfrow=c(1,1), mar=c(10,5,5,0.3), las=1, bty="l") 
    bp=barplot(freqs, main=gsub("meta to ","",ttype), xaxt="n", yaxt='n', xlab="", col="#ffe100", border="#000000", width=width, space=0, lwd=4, ylim=ylim, cex.main=2.5)
    xaxis=breaks
    sxseq=c(bp[breaks==-4]-width/2 , bp[breaks==-2]-width/2)
    text("Metastasis Potential", x=mean(bp), y=ymin-(ymax-ymin)/5, cex=2, xpd=T)
    text(xaxis, x=c(bp-width/2,max(bp)+width/2), y=ymin-(ymax-ymin)/23, xpd=T, srt=0, cex=ifelse(breaks %in% c(-4,-2) , 1.8 , 1.5), col=ifelse(breaks %in% c(-4,-2) , tcols[match(breaks,names(tcols))] , "black"))
    axis(side=2, at=yseq, label=NA, xpd=T, adj=1, line=-1.1, lwd=1.2, lwd.ticks=1.2, tcl=-0.3)
    mtext(yseq, at=yseq, side=2, line=-0.6, cex=1.7)
    text("Frequency", x=min(bp)-4.5, y=mean(yseq), cex=2, srt=90, xpd=T)
    segments(x0=sxseq , x1=sxseq , y0=ymin, y1=ymax, lty=2, lwd=3.5, col=c('#000dff','#ff0000'), xpd=T)
    text(counts, x=bp, y=freqs+(ymax-ymin)/25, cex=1.8, xpd=T)
    dev.off()
}




########################################
## Metastasis Potential Summary per target organ
########################################

## filter cancer type
ectypes=names(table(sinfo$`Cancer Type`)[table(sinfo$`Cancer Type`)>=10])
sinfo=sinfo[sinfo$`Cancer Type` %in% ectypes,]
sinfo$`Cancer Type`=factor(sinfo$`Cancer Type`, levels=names(sort(table(sinfo$`Cancer Type`)[table(sinfo$`Cancer Type`)>=10],decreasing = T)))


## correlation test
ttypes=colnames(sinfo)[grepl("meta to",colnames(sinfo))]
ctypes=names(sort(table(sinfo$`Cancer Type`),decreasing = T))

hmtx=t(as.data.frame(lapply(split(sinfo[,ttypes] , sinfo$`Cancer Type`), function(df) apply(df,2,mean))))
rownames(hmtx)=gsub("\\."," ",rownames(hmtx))
colnames(hmtx)=str_to_title(gsub("meta to ","",colnames(hmtx)))

tmtx=round(hmtx,3)

counts=table(sinfo$`Cancer Type`)
counts=as.numeric(counts[match(rownames(hmtx),names(counts))])


## viz
color.hp=colorRamp2(c(min(hmtx[hmtx!=100]),-3,max(hmtx[hmtx!=100])), c('#403dff','#d4d4d4','#ff3d3d'))
ran=rowAnnotation(`Count of\nSample` = anno_barplot(counts, width = unit(5, "cm"), gp=gpar(fill="#fff239", col="#ffee00", lwd=2.5), border=F, axis=F, add_numbers=T, numbers_gp=gpar(fontsize=30, col='#000000'), numbers_rot=0, numbers_offset=unit(0.3,"cm")),
                            annotation_name_offset=c(`Count of\nSample`="0.5cm"), annotation_name_gp=gpar(fontsize=30, col="#000000"), annotation_name_side='top', annotation_name_rot=0)
ran.gap = rowAnnotation(gap_dummy = anno_barplot(rep(0, length(counts)), height = unit(2.5, "cm"), width = unit(1, "cm"), gp = gpar(fill = "#ffffff", col = "#ffffff"), border = FALSE, axis=F),show_annotation_name = FALSE)
hp=Heatmap(hmtx, show_heatmap_legend=F, col=color.hp, border=F, cluster_columns=T, cluster_rows=T, show_row_names=T, show_column_names=T, cluster_column_slices=F, show_row_dend=T, column_names_rot=90, column_gap=unit(0.7,"cm"), row_labels=rownames(hmtx), column_labels=gsub("meta to ","",colnames(hmtx)),
        width=ncol(hmtx)*unit(3.5,"cm"), height=nrow(hmtx)*unit(2,"cm"), rect_gp = gpar(col = "#ffffff", lwd = 2.2, lty=1), 
        row_names_gp = gpar(fontsize=42), column_title_gp = gpar(fontsize=40, fontface="bold"), column_names_gp = gpar(fontsize=47),
        column_dend_height=unit(3,"cm"), row_dend_width=unit(3,"cm"), row_dend_gp=gpar(lwd=3, col="#808080"), column_dend_gp=gpar(lwd=3, col="#808080"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(tmtx[i, j], x, y, gp = gpar(fontsize = 28, col = "#ffffff")
        )
        },
        right_annotation = c(ran,ran.gap))
tiff(filename=paste0(dir_fig,"Main Figure 1(A).tif"), width=50, height=50, units = 'cm',res=300)
draw(hp, ht_gap=unit(0.7,"cm"), padding=unit(c(0, 0, 0, 7), "cm")) 
dev.off()


# legend
lg.hp=Legend(at=c(min(hmtx[hmtx!=100]),-3,max(hmtx[hmtx!=100])), labels=c("Minimum",-3,"Maximum"), col_fun = color.hp, title = "Average\n metastasis potential", border='#808080', title_position="topcenter", size=unit(4,"cm"))
tiff(filename=paste0(dir_fig,"Main Figure 1(A)_legend.tiff"), width=5, height=8, units = 'cm',res=300)
draw(packLegend(list=list(lg.hp))) 
dev.off() 



########################################
## Difference in metastasis potential by target organ
########################################
psinfo=as.data.frame(sinfo)

ctypes=names(table(psinfo$`Cancer Type`))[table(psinfo$`Cancer Type`)>=10]
psinfo=psinfo[psinfo$`Cancer Type` %in% ctypes,]

cols=colnames(psinfo)[grepl("^meta to ",colnames(psinfo))]
psinfo=psinfo[,cols,drop=F]

pdf=do.call(rbind,lapply(colnames(psinfo), function(ttype) data.frame(target=str_to_title(gsub("meta to ","",ttype)),mp=psinfo[,ttype])))
pdf$target=ifelse(pdf$target=="All5","all5",pdf$target)
pdf=pdf[!is.na(pdf$mp),]
pdf$mtype=ifelse(pdf$mp<=(-4),"Non Metastatic",ifelse(pdf$mp>=(-2),"Metastatic","Weakly Metastatic"))
mtypes=unique(pdf$mtype)

# sorting
mres=aggregate(mp ~ target, data = pdf, median)
targets=mres$target[order(mres$mp)]
pdf$target=factor(pdf$target, levels=targets)
pdf=pdf[order(pdf$target),]
targets=levels(pdf$target)
xseq=seq(1,length(targets),1)


## test
lm=lm(mp~target, data=pdf)
ares=Anova(lm, type=2)
pv=ares[1,"Pr(>F)"]

tres=TukeyHSD(aov(mp~target, data=pdf))


## viz
pcols=setNames(c("#000dff","#a4a4a4","#ff0000"),c("Non Metastatic","Weakly Metastatic","Metastatic"))
pfils=setNames(c("#403dff73","#d4d4d473","#ff3d3d73"),c("Non Metastatic","Weakly Metastatic","Metastatic"))

ymin=min(pdf$mp)
ymax=max(pdf$mp)
yrange=ymax-ymin
ylim=c(ymin,ymax+yrange*0.05)

tiff(filename=paste0(dir_fig,"Main Figure 1(B).tif"), width=7.5, height=9.5, units="cm", res=300)
par(mfrow=c(1,1), mar=c(5,4.5,5.5,1), mgp=c(2.8,0.5,0), tck=-0.03, las=1, bty="l") 
bp=boxplot(pdf$mp~pdf$target, cex.lab=1.1, cex.main=1.3, cex.axis=0.9, lwd=1.2, boxwex=0.6, frame=T, xaxt='n', yaxt='n', border="#5a5a5aff", col=NA, outline=F, main="", xlab="", ylab="Metastasis potential", ylim=ylim)
title(main="Target organ", cex.main=0.85, line=4.3)
axis(2, at=c(-4,-2,0,2), labels=c(-4,-2,0,2), las=1, cex.axis=0.9)
# for(target in targets) {
#     for(mtype in mtypes) {
#         stripchart(pdf$mp[pdf$target==target & pdf$mtype==mtype]~pdf$target[pdf$target==target & pdf$mtype==mtype], at=which(targets==target), cex=1, lwd=1.4, pch=21, col=pcols[mtype], bg=pfils[mtype], method="jitter", jitter=0.24, vertical=TRUE, add = TRUE)
#     }
# }
pvlab=if(pv<2.2e-16 | pv==0) expression(italic(p) < 2.2 %*% 10^-16) else if(pv<0.00001) bquote(italic(p) == .(round(pv / 10^floor(log10(pv)), 4)) %*% 10^.(floor(log10(pv)))) else bquote(italic(p) == .(round(pv,4)))
ytop=par("usr")[4]
segments(y0=ytop+yrange*0.3,y1=ytop+yrange*0.3,x0=min(xseq),x1=max(xseq),col="#000000ff",lty=1,lwd=1.2,xpd=NA)
segments(y0=rep(ytop+yrange*0.27,2),y1=rep(ytop+yrange*0.3,2),x0=c(min(xseq),max(xseq)),x1=c(min(xseq),max(xseq)),col="#000000ff",lty=1,lwd=1.2,xpd=NA)
text(labels=pvlab, x=mean(xseq), adj=0.5, cex=0.7, y=ytop+yrange*0.37, xpd=NA, srt=0)
segments(y0=ytop+yrange*0.17,y1=ytop+yrange*0.17,x0=min(xseq),x1=max(xseq),col="#000000ff",lty=1,lwd=1.2,xpd=NA)
text(labels="**", x=mean(xseq), adj=0.5, cex=0.85, y=ytop+yrange*0.2, xpd=NA, srt=0)
segments(y0=ytop+yrange*0.05,y1=ytop+yrange*0.05,x0=min(xseq),x1=max(xseq)-1,col="#000000ff",lty=1,lwd=1.2,xpd=NA)
text(labels="**", x=mean(c(min(xseq),max(xseq)-1)), adj=0.5, cex=0.85, y=ytop+yrange*0.08, xpd=NA, srt=0)
text(labels=targets, x=xseq, adj=1, cex=0.7, y=ymin-yrange*0.08, xpd=NA, srt=60)
dev.off()



########################################
## Number of metastatic target organs
########################################

ctypes=names(table(sinfo$`Cancer Type`)[table(sinfo$`Cancer Type`)>=10])
psinfo=sinfo[sinfo$`Cancer Type` %in% ctypes , ] 

ctypes=names(sort(table(psinfo$`Cancer Type`),decreasing = T))
ttypes=c(0:5)

hmtx=matrix(NA, nrow=length(ctypes),ncol=length(ttypes), dimnames=list(ctypes,ttypes))

for(ctype in ctypes) {
        for(ttype in ttypes) {
                esinfo=psinfo[psinfo$`Cancer Type`==ctype , ]
                counts=apply(esinfo[,grep("meta type to ",colnames(psinfo))], 1, function(t) sum(t=="Metastatic"))  
                values=esinfo$`meta to all5`[counts==as.character(ttype)]                     
                hmtx[ctype,as.character(ttype)]=ifelse(sum(!is.na(values))==0,100,mean(values,na.rm=T))
        }
}

tmtx=round(hmtx,3)
tmtx[tmtx==100]="None"
tmtx_b=tmtx
tmtx_w=tmtx
tmtx_b[tmtx_b!="None"]=""
tmtx_w[tmtx_w=="None"]=""

color.hp=colorRamp2(c(min(hmtx[hmtx!=100]),-3,max(hmtx[hmtx!=100]),100), c('#403dff','#d4d4d4','#ff3d3d','#ffffffff'))
counts=table(psinfo$`Cancer Type`)
counts=as.numeric(counts[match(rownames(hmtx),names(counts))])
ran=rowAnnotation(`Count of\nSample` = anno_barplot(counts, width = unit(5, "cm"), gp=gpar(fill="#fff239", col="#ffee00", lwd=2.5), border=F, axis=F, add_numbers=T, numbers_gp=gpar(fontsize=30, col='#000000'), numbers_rot=0, numbers_offset=unit(0.3,"cm")),
                            annotation_name_offset=c(`Couunt of\nSample`="0.5cm"), annotation_name_gp=gpar(fontsize=26, col="#000000"), annotation_name_side='top', annotation_name_rot=0)
ran.gap = rowAnnotation(gap_dummy = anno_barplot(rep(0, length(counts)), height = unit(2.5, "cm"), width = unit(1, "cm"),, gp = gpar(fill = "#ffffff", col = "#ffffff"), border = FALSE, axis=F),show_annotation_name = FALSE)
hp=Heatmap(hmtx, show_heatmap_legend=F, col=color.hp, border=F, cluster_columns=F, cluster_rows=F, show_row_names=T, show_column_names=T, cluster_column_slices=F, show_column_dend=F, row_title=NULL, column_title="Metastasis to",
        width=ncol(hmtx)*unit(3.5,"cm"), height=nrow(hmtx)*unit(2,"cm"), rect_gp = gpar(col = "#808080", lwd = 2.2, lty=1), column_names_rot=0,
        row_names_gp = gpar(fontsize=42), column_names_gp = gpar(fontsize=47), column_title_gp = gpar(fontsize=40),
        row_dend_gp=gpar(lwd=3, col="#808080"), row_dend_width=unit(4,"cm"), column_dend_gp=gpar(lwd=3, col="#808080"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(tmtx_b[i, j], x, y, gp = gpar(fontsize = 33, col = "#000000ff"))
                grid.text(tmtx_w[i, j], x, y, gp = gpar(fontsize = 33, col = "#ffffffff"))
                },
        right_annotation = c(ran,ran.gap)
        )

tiff(filename=paste0(dir_fig,"Main Figure 1(D).tif"), width=50, height=45, units='cm',res=300)
draw(hp, ht_gap=unit(0.7,"cm"), padding=unit(c(0, 0, 0, 6), "cm")) 
dev.off()


## legend
lg.hp=Legend(at=c(min(hmtx[hmtx!=100]),-3,max(hmtx[hmtx!=100])), labels=c(round(min(hmtx[hmtx!=100]),1),"-3",round(max(hmtx[hmtx!=100]),1)), col_fun = color.hp, title = "Average\n metastasis potential", border='#808080', title_position="topcenter", size=unit(4,"cm"))
tiff(filename=paste0(dir_fig,"Main Figure 1(D)_legend.tif"), width=5, height=8, units = 'cm',res=300)
draw(packLegend(list=list(lg.hp))) 
dev.off() 



########################################
## Confounding factor Identification 
########################################
psinfo=sinfo

ctypes=names(table(psinfo$`Cancer Type`))[table(psinfo$`Cancer Type`)>=9]
psinfo=psinfo[psinfo$`Cancer Type` %in% ctypes,]

psinfo$`Age Group`=floor(psinfo$Age*0.1)*10
psinfo$`Age Group`=ifelse(is.na(psinfo$`Age Group`),psinfo$`Age Group`,paste0(psinfo$`Age Group`,"s"))
"Cancer Type"s=c("Cell Line Origin","Age Group","Sex","Race")


## viz
mdf=summarise(group_by(psinfo,`Cancer Type`), median=median(`meta to all5`))
psinfo$`Cancer Type`=factor(psinfo$`Cancer Type`, levels=mdf$`Cancer Type`[order(mdf$median, decreasing=T)])
psinfo$`meta type to all5`=factor(psinfo$`meta type to all5`, levels=c("Non Metastatic","Metastatic","Weakly Metastatic (Low Confidence)"))
psinfo$`Cell Line Origin`=factor(psinfo$`Cell Line Origin`, levels=c("Primary","Metastatic"))
mtypes=unique(psinfo$`meta type to all5`)
pcols=setNames(c("#000dff","#a4a4a4","#ff0000"),c("Non Metastatic","Weakly Metastatic (Low Confidence)","Metastatic"))
pfils=setNames(c("#403dff73","#d4d4d473","#ff3d3d73"),c("Non Metastatic","Weakly Metastatic (Low Confidence)","Metastatic"))

for("Cancer Type" in "Cancer Type"s) {
    for(ctype in ctypes) {
        gdf=psinfo[psinfo$`Cancer Type`==ctype,]
        gdf=gdf[,c("meta to all5","meta type to all5","Cancer Type")]
        gdf=gdf[!is.na(gdf[,"Cancer Type"]) & gdf[,"Cancer Type"]!="Unknown" & gdf[,"Cancer Type"]!="",]
        gdf[,"Cancer Type"]=factor(gdf[,"Cancer Type"], levels=sort(unique(gdf[,"Cancer Type"])))    
        gdf=gdf[order(gdf[,"Cancer Type"]),]
        ftypes=unique(gdf[,"Cancer Type"])

        if(length(unique(gdf[,"Cancer Type"]))<2) next

        ymin=min(gdf$`meta to all5`)
        ymax=max(gdf$`meta to all5`)*1.6
        ylim=c(ymin,ymax)
        xseq=seq(1,length(ftypes),1)

        tiff(filename=paste0(dir_fig,"1.14. EDA(pan)/1.14.2.1. confounding_factor_identification_per_cancer type/","Cancer Type","_",ctype,".tif"), width=5+length(ftypes)*2.2, height=27, units="cm", res=300)
        par(mfrow=c(1,1), mar=c(25,10,10,2), mgp=c(5,1,0), las=1, bty="l") 
        bp=boxplot(gdf$`meta to all5`~gdf[,"Cancer Type"],  cex.lab=2.5, cex.main=2.9,  cex.axis=2.2, lwd=3, boxwex=0.47, frame=T, xaxt='n', border="#5a5a5aff", col=NA, outline=F, main="Cancer Type", xlab="", ylab="Metastasis Potential", ylim=ylim)
        for(ftype in ftypes) {
            for(mtype in mtypes) {
                stripchart(gdf$`meta to all5`[gdf[,"Cancer Type"]==ftype & gdf$`meta type to all5`==mtype]~gdf[,"Cancer Type"][gdf[,"Cancer Type"]==ftype & gdf$`meta type to all5`==mtype], cex=1.5, lwd=1.4, pch=21, col=pcols[mtype], bg=pfils[mtype], method="jitter", jitter=0.19, vertical=TRUE, add = TRUE)
            }
        }
        lm=lm(gdf[,"meta to all5"]~gdf[,"Cancer Type"])
        ares=Anova(lm,type=2)    
        pv=ares[1,"Pr(>F)"]
        if (pv<0.0001) pv=bquote(.(round( pv / 10^floor(log10(pv)), 4)) %*% 10^.(floor(log10(pv)))) else pv=round(pv,4)
        segments(y0=ymax+(ymax-ymin)*0.08, y1=ymax+(ymax-ymin)*0.08 , x0=min(xseq) , x1=max(xseq),col="#000000ff",lty=1,lwd=3, xpd=T)
        segments(y0=rep(c(ymax+(ymax-ymin)*0.08),2),y1=rep(ymax+(ymax-ymin)*0.06,2),x0=c(min(xseq),max(xseq)),x1=c(min(xseq),max(xseq)),col="#000000ff",lty=1,lwd=3, xpd=T)
        text(labels=pv, x=mean(xseq), adj=0.5, cex=2.5, y=ymax+(ymax-ymin)*0.17, xpd=T, srt=0)
        text(labels=gsub("_"," ",ftypes), x=xseq, adj=1, cex=2, y=ymin-0.5, xpd=T, srt=60)
        dev.off()
    }
}



########################################
## Metastasis potential by Tumor lineage
########################################
psinfo=sinfo

ctypes=names(table(psinfo$`Cancer Type`))[table(psinfo$`Cancer Type`)>=10]
psinfo=psinfo[psinfo$`Cancer Type` %in% ctypes,]

psinfo$`Age Group`=floor(psinfo$Age*0.1)*10
psinfo$`Age Group`=ifelse(is.na(psinfo$`Age Group`),psinfo$`Age Group`,paste0(psinfo$`Age Group`,"s"))

ldf=psinfo[,c("meta to all5","meta type to all5","Cancer Type")] # data frame for linear modeling
ldf=ldf[!is.na(ldf[,"Cancer Type"]) & ldf[,"Cancer Type"]!="Unknown" & ldf[,"Cancer Type"]!="",]
ldf[,"Cancer Type"]=factor(ldf[,"Cancer Type"])
lm=lm(ldf[,"meta to all5"]~ldf[,"Cancer Type"]) # linear modeling
ares=Anova(lm, type=2)
pv=ares[1,"Pr(>F)"] # 6.899297e-08


## viz
mdf=summarise(group_by(psinfo,`Cancer Type`), median=median(`meta to all5`)) # median value
psinfo$`Cancer Type`=factor(psinfo$`Cancer Type`, levels=mdf$`Cancer Type`[order(mdf$median, decreasing=T)])
psinfo$`meta type to all5`=factor(psinfo$`meta type to all5`, levels=c("Non Metastatic","Metastatic","Weakly Metastatic (Low Confidence)"))
psinfo$`Cell Line Origin`=factor(psinfo$`Cell Line Origin`, levels=c("Primary","Metastatic"))
mtypes=unique(psinfo$`meta type to all5`)
pcols=setNames(c("#000dff","#a4a4a4","#ff0000"),c("Non Metastatic","Weakly Metastatic (Low Confidence)","Metastatic"))
pfils=setNames(c("#403dff73","#d4d4d473","#ff3d3d73"),c("Non Metastatic","Weakly Metastatic (Low Confidence)","Metastatic"))

gdf=psinfo[,c("meta to all5","meta type to all5","Cancer Type")]
gdf=gdf[!is.na(gdf[,"Cancer Type"]) & gdf[,"Cancer Type"]!="Unknown" & gdf[,"Cancer Type"]!="",]
gdf[,"Cancer Type"]=factor(gdf[,"Cancer Type"], levels=sort(unique(gdf[,"Cancer Type"])))    
gdf=gdf[order(gdf[,"Cancer Type"]),]
ftypes=unique(gdf[,"Cancer Type"])

ymin=min(gdf$`meta to all5`)
ymax=max(gdf$`meta to all5`)
ylim=c(ymin,ymax*1.2)
xseq=seq(1,length(ftypes),1)

tiff(filename=paste0(dir_fig,"Main Figure 1(E).tif"), width=5+length(ftypes)*2.2, height=27, units="cm", res=300)
par(mfrow=c(1,1), mar=c(25,10,10,2), mgp=c(5,1,0), las=1, bty="l") 
bp=boxplot(gdf$`meta to all5`~gdf[,"Cancer Type"],  cex.lab=2.5, cex.main=2.9,  cex.axis=2.2, lwd=3, boxwex=0.7, frame=T, xaxt='n', border="#5a5a5aff", col=NA, outline=F, main="Tumor Lineage", xlab="", ylab="Metastasis potential", ylim=ylim)
for(ftype in ftypes) {
    for(mtype in mtypes) {
        stripchart(gdf$`meta to all5`[gdf[,"Cancer Type"]==ftype & gdf$`meta type to all5`==mtype]~gdf[,"Cancer Type"][gdf[,"Cancer Type"]==ftype & gdf$`meta type to all5`==mtype], cex=1.5, lwd=1.4, pch=21, col=pcols[mtype], bg=pfils[mtype], method="jitter", jitter=0.17, vertical=TRUE, add = TRUE)
    }
}
if (pv < 0.00001) {pv_label = bquote(.(round(pv / 10^floor(log10(pv)), 4)) %*% 10^.(floor(log10(pv))))} else { pv_label = round(pv, 4)}
segments(y0=ymax*1.2,y1=ymax*1.2,x0=min(xseq),x1=max(xseq),col="#000000ff",lty=1,lwd=3)
segments(y0=rep(ymax*1.2,2),y1=rep(ymax*1.05,2),x0=c(min(xseq),max(xseq)),x1=c(min(xseq),max(xseq)),col="#000000ff",lty=1,lwd=3)
text(labels=pv_label, x=mean(xseq), adj=0.5, cex=2.8, y=ymax*1.5, xpd=T, srt=0)
text(labels=gsub("_"," ",ftypes), x=xseq, adj=1, cex=2, y=ymin-0.5, xpd=T, srt=60)
dev.off()



########################################
## Metastasis potential by Cell line origin
########################################
psinfo=sinfo

ctypes=names(table(psinfo$`Cell Line Origin`))[table(psinfo$`Cell Line Origin`)>=10]
psinfo=psinfo[psinfo$`Cell Line Origin` %in% ctypes,]

psinfo$`Age Group`=floor(psinfo$Age*0.1)*10
psinfo$`Age Group`=ifelse(is.na(psinfo$`Age Group`),psinfo$`Age Group`,paste0(psinfo$`Age Group`,"s"))

ldf=psinfo[,c("meta to all5","meta type to all5","Cell Line Origin")] # data frame for linear modeling
ldf=ldf[!is.na(ldf[,"Cell Line Origin"]) & ldf[,"Cell Line Origin"]!="Unknown" & ldf[,"Cell Line Origin"]!="",]
ldf[,"Cell Line Origin"]=factor(ldf[,"Cell Line Origin"])
lm=lm(ldf[,"meta to all5"]~ldf[,"Cell Line Origin"]) # linear modeling
ares=Anova(lm, type=2)
pv=ares[1,"Pr(>F)"] # 6.899297e-08


## viz
mdf=summarise(group_by(psinfo,`Cell Line Origin`), median=median(`meta to all5`)) # median value
psinfo$`Cell Line Origin`=factor(psinfo$`Cell Line Origin`, levels=mdf$`Cell Line Origin`[order(mdf$median, decreasing=T)])
psinfo$`meta type to all5`=factor(psinfo$`meta type to all5`, levels=c("Non Metastatic","Metastatic","Weakly Metastatic (Low Confidence)"))
psinfo$`Cell Line Origin`=factor(psinfo$`Cell Line Origin`, levels=c("Primary","Metastatic"))
mtypes=unique(psinfo$`meta type to all5`)
pcols=setNames(c("#000dff","#a4a4a4","#ff0000"),c("Non Metastatic","Weakly Metastatic (Low Confidence)","Metastatic"))
pfils=setNames(c("#403dff73","#d4d4d473","#ff3d3d73"),c("Non Metastatic","Weakly Metastatic (Low Confidence)","Metastatic"))

gdf=psinfo[,c("meta to all5","meta type to all5","Cell Line Origin")]
gdf=gdf[!is.na(gdf[,"Cell Line Origin"]) & gdf[,"Cell Line Origin"]!="Unknown" & gdf[,"Cell Line Origin"]!="",]
gdf[,"Cell Line Origin"]=factor(gdf[,"Cell Line Origin"], levels=sort(unique(gdf[,"Cell Line Origin"])))    
gdf=gdf[order(gdf[,"Cell Line Origin"]),]
ftypes=unique(gdf[,"Cell Line Origin"])

ymin=min(gdf$`meta to all5`)
ymax=max(gdf$`meta to all5`)
ylim=c(ymin,ymax*1.2)
xseq=seq(1,length(ftypes),1)

tiff(filename=paste0(dir_fig,"Main Figure 1(F).tif"), width=13, height=27, units="cm", res=300)
par(mfrow=c(1,1), mar=c(25,12,10,5), mgp=c(5,1,0), las=1, bty="l") 
bp=boxplot(gdf$`meta to all5`~gdf[,"Cell Line Origin"],  cex.lab=2.5, cex.main=2.9,  cex.axis=2.2, lwd=3, boxwex=0.7, frame=T, xaxt='n', border="#5a5a5aff", col=NA, outline=F, main="Cell line origin", xlab="", ylab="Metastasis potential", ylim=ylim)
for(ftype in ftypes) {
    for(mtype in mtypes) {
        stripchart(gdf$`meta to all5`[gdf[,"Cell Line Origin"]==ftype & gdf$`meta type to all5`==mtype]~gdf[,"Cell Line Origin"][gdf[,"Cell Line Origin"]==ftype & gdf$`meta type to all5`==mtype], cex=1.5, lwd=1.4, pch=21, col=pcols[mtype], bg=pfils[mtype], method="jitter", jitter=0.17, vertical=TRUE, add = TRUE)
    }
}
if (pv < 0.00001) {pv_label = bquote(.(round(pv / 10^floor(log10(pv)), 4)) %*% 10^.(floor(log10(pv))))} else { pv_label = format(round(pv,4), scientific = FALSE)}
segments(y0=ymax*1.2,y1=ymax*1.2,x0=min(xseq),x1=max(xseq),col="#000000ff",lty=1,lwd=3)
segments(y0=rep(ymax*1.2,2),y1=rep(ymax*1.05,2),x0=c(min(xseq),max(xseq)),x1=c(min(xseq),max(xseq)),col="#000000ff",lty=1,lwd=3)
text(labels=pv_label, x=mean(xseq), adj=0.5, cex=2.8, y=ymax*1.5, xpd=T, srt=0)
text(labels=gsub("_"," ",ftypes), x=xseq, adj=1, cex=2, y=ymin-0.5, xpd=T, srt=60)
dev.off()



########################################
## Fraction of Metastasis Potential Types by Cell Line Origin
########################################
pdf = data.frame(origin = as.character(sinfo[["Cell Line Origin"]]), potential = as.character(sinfo[["meta type to all5"]]))
pdf$origin[pdf$origin==""]="Unknown"
pdf = pdf[!is.na(pdf$origin) & !is.na(pdf$potential), ]

origin_order = c("Metastatic", "Unknown", "Primary")
potential_order = c("Metastatic", "Weakly Metastatic (Low Confidence)", "Non Metastatic")
pdf = pdf[pdf$origin %in% origin_order & pdf$potential %in% potential_order, ]

pdf$origin = factor(pdf$origin, levels = origin_order)
pdf$potential = factor(pdf$potential, levels = potential_order)

tab = as.data.frame(table(origin = pdf$origin, potential = pdf$potential), stringsAsFactors = FALSE)
tab$origin = as.character(tab$origin)
tab$potential = as.character(tab$potential)

origin_sum = aggregate(Freq ~ origin, data = tab, sum)
origin_sum$origin = factor(as.character(origin_sum$origin), levels = origin_order)
origin_sum = origin_sum[order(origin_sum$origin), ]
origin_sum = origin_sum[origin_sum$Freq > 0, ]
origin_sum$origin = as.character(origin_sum$origin)

total_n = sum(origin_sum$Freq)
gap = 0.035

origin_sum$pct = 100 * origin_sum$Freq / total_n
origin_sum$height = origin_sum$Freq / total_n * (1 - gap * (nrow(origin_sum) - 1))
origin_sum$ymin = rep(NA_real_, nrow(origin_sum))
origin_sum$ymax = rep(NA_real_, nrow(origin_sum))

y_now = 1
for(i in seq_len(nrow(origin_sum))) {
    origin_sum$ymax[i] = y_now
    origin_sum$ymin[i] = y_now - origin_sum$height[i]
    y_now = origin_sum$ymin[i] - gap
}
x_left0 = 0
x_left1 = 0.30
x_right0 = 1.05
x_right1 = 1.35

origin_cols = setNames(rep("grey90", length(origin_order)), origin_order)
origin_cols[names(origin_cols) %in% "Metastatic"] = "grey10"
origin_cols[names(origin_cols) %in% "Primary"] = "grey92"
origin_cols[names(origin_cols) %in% "Unknown"] = "white"

potential_cols = setNames(rep("grey80", length(potential_order)), potential_order)
potential_cols[1] = "#ff3b3b"
if(length(potential_order) >= 2) potential_cols[2] = "grey85"
if(length(potential_order) >= 3) potential_cols[3] = "#3448ff"

origin_rect = origin_sum
origin_rect$xmin = x_left0
origin_rect$xmax = x_left1
origin_rect$y = (origin_rect$ymin + origin_rect$ymax) / 2
origin_rect$label = paste0(origin_rect$origin, "\n", sprintf("%.1f%%", origin_rect$pct))
origin_rect$fill_col = origin_cols[origin_rect$origin]
origin_rect$label_x = (x_left0 + x_left1) / 2
origin_rect$hjust = 0.5
origin_rect$text_col = ifelse(origin_rect$fill_col == "grey10", "white", "black")

rdf = data.frame()
for(i in seq_len(nrow(origin_sum))) {
    tmp = data.frame(origin = origin_sum$origin[i], x = c(x_left1, x_right0, x_right0, x_left1), y = c(origin_sum$ymin[i], origin_sum$ymin[i], origin_sum$ymax[i], origin_sum$ymax[i]))
    rdf = rbind(rdf, tmp)
}
pot_rect = data.frame()
for(i in seq_len(nrow(origin_sum))) {
    os = origin_sum[i, ]
    tmp = tab[tab$origin == os$origin, ]
    tmp = tmp[match(potential_order, tmp$potential), ]
    tmp = tmp[!is.na(tmp$potential), ]
    y_top = os$ymax
    for(j in seq_len(nrow(tmp))) {
        h = os$height * tmp$Freq[j] / os$Freq
        y_bottom = y_top - h
        
        if(tmp$Freq[j] > 0) {
            add = data.frame(origin = os$origin, potential = tmp$potential[j], Freq = tmp$Freq[j], pct = 100 * tmp$Freq[j] / os$Freq, xmin = x_right0, xmax = x_right1, ymin = y_bottom, ymax = y_top, y = (y_bottom + y_top) / 2, fill_col=potential_cols[tmp$potential[j]])
            pot_rect = rbind(pot_rect, add)
        }
    y_top = y_bottom
    }
}

pot_rect$label = sprintf("%.1f%%", pot_rect$pct)
pot_rect$text_col = ifelse(pot_rect$fill_col %in% c("#ff3b3b", "#3448ff"), "white", "black")

gp = ggplot() +
    geom_polygon(data = rdf, aes(x = x, y = y, group = origin), fill = "grey78", alpha = 0.65, color = NA) +
    geom_rect(data = origin_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_col), color = "grey45", size = 0.3) +
    geom_rect(data = pot_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_col), color = "grey45", size = 0.3) +
    geom_text(data = origin_rect, aes(x = label_x, y = y, label = label, hjust = hjust, color = text_col), size = 3.2, lineheight = 0.95) +
    geom_text(data = pot_rect, aes(x = (xmin + xmax) / 2, y = y, label = label, color = text_col), size = 3.2) +
    annotate("text", x = (x_left0 + x_left1) / 2, y = 1.07, label = "Origin", size = 4) +
    annotate("text", x = (x_right0 + x_right1) / 2, y = 1.07, label = "Potential", size = 4) +
    scale_fill_identity() +
    scale_color_identity() +
    coord_cartesian(xlim = c(-0.05, 1.48), ylim = c(0, 1.11), clip = "off") +
    theme_void()

ggsave(paste0(dir_fig,"Main Figure 1(G).tiff"), plot = gp, width = 6, height = 6, units = "in", dpi = 600, bg = "white")



########################################
## Metastasis Potential Type per Cancer Type (split)
########################################
psinfo=sinfo
psinfo$`Cancer Type`=factor(psinfo$`Cancer Type`, levels=sort(unique(psinfo$`Cancer Type`), decreasing=T))
psinfo=psinfo[order(psinfo$`Cancer Type`, decreasing=T),]

cancer_types=names(table(psinfo$`Cancer Type`)[table(psinfo$`Cancer Type`)>=10])
gdf=table(psinfo[,c("meta type to all5","Cancer Type")])
gdf=gdf[,colnames(gdf) %in% cancer_types]
gdf=apply(gdf,2,function(c) c/sum(c)*100)
gdf=gdf[,order(gdf["Non Metastatic",], decreasing=T)]
sorts=names(sort(gdf["Metastatic",]))


for(origin in unique(psinfo$`Cell Line Origin`)) {
    esinfo=psinfo[psinfo$`Cell Line Origin`==origin & psinfo$`Cancer Type` %in% cancer_types,]
    esinfo[,"meta type to all5"]=factor(esinfo[,"meta type to all5"], levels=c("Non Metastatic","Weakly Metastatic (Low Confidence)","Metastatic"))
    gdf=table(esinfo[,c("meta type to all5","Cancer Type")])
    gdf=gdf[,colnames(gdf) %in% cancer_types]
    counts=apply(gdf,2,sum)
    gdf=apply(gdf,2,function(c) c/sum(c)*100)
    gdf=gdf[,sorts]
    ctypes=colnames(gdf)
    gdf[is.na(gdf)]=0
    ncounts=table(esinfo$`Cancer Type`[esinfo[,"meta type to all5"]=="Non Metastatic"])[ctypes] # non meta counts
    mcounts=table(esinfo$`Cancer Type`[esinfo[,"meta type to all5"]=="Metastatic"])[ctypes] # meta counts
    wcounts=table(esinfo$`Cancer Type`[esinfo[,"meta type to all5"]=="Weakly Metastatic (Low Confidence)"])[ctypes] # weakly meta counts

    xmax=100
    xseq=seq(0,xmax,20)
    xrange=c(min(xseq),max(xseq))
    yseq=1:length(ctypes)
    yrange=c(min(yseq),max(yseq))
    pcols=setNames(c("#ff3d3d","#c7c7c7","#403dff") , c("Metastatic","Weakly Metastatic (Low Confidence)","Non Metastatic"))

    tiff(filename=paste0(dir_fig,"Main Figure 1(H)_",origin,".tif"), width=13, height=10, units="cm", res=300)
    par(mfrow=c(1,1), mar=c(3,15,2.2,2.7), las=1, bty="l")
    bp=barplot(gdf, main=origin, horiz=T, xaxt="n", yaxt="n", xlab="", ylab="", col=pcols[match(rownames(gdf),names(pcols))], border='black', width=0.4, lwd=0.5, space=0.2, xlim=xrange)
    axis(1, at=xseq, label=NA, xpd=T, adj=1, line=-0.4, lwd=1.2, tcl=-0.2)
    mtext(xseq, at=xseq, side=1, line=-0.4, cex=0.85) 
    text(labels=paste0(ctypes," (",counts[ctypes],")"), x=rep(0,length(bp))-(xmax*0.2), y=bp, col='black', cex=1, adj=1, xpd=T)
    text(labels=ifelse(ncounts>0,ncounts,""), x=rep(0,length(bp))-(xmax*0.02), y=bp, col='#403dff', cex=1, adj=1, xpd=T)
    text(labels=ifelse(mcounts>0,mcounts,""), x=rep(xmax*1.02,length(bp)), y=bp, col='#ff3d3d', cex=1, adj=0, xpd=T)
    text(labels=ifelse(wcounts>0,wcounts,""), x=100-(gdf["Metastatic",])*0.98, y=bp, col='#dbdbdb', cex=1, adj=0, xpd=T)
    text(labels="Percentage (%)", adj=0.5, cex=1, x=mean(xseq), y=-1, xpd=T)
    dev.off()
}



########################################
## Classification of cancer cell lines by metastatic category
########################################
tb=sinfo[sinfo$`meta type to all5`!="Weakly Metastatic (Low Confidence)",c("meta type to all5","Cancer Type")]
tb=table(tb)[c("Non Metastatic","Metastatic"),names(sort(apply(table(tb),2,sum),decreasing = T))]
tb
#                  Cancer Type
# meta type to all5 Lung Skin Esophagus or Stomach Pancreas
#    Non Metastatic    3    0                    1        0
#    Metastatic       55   33                   25       23
#                  Cancer Type
# meta type to all5 Ovary or Fallopian Tube Bladder or Urinary Tract Uterus Bowel
#    Non Metastatic                       1                        1      1     2
#    Metastatic                          15                       14     14    12
#                  Cancer Type
# meta type to all5 Head and Neck CNS or Brain Bone Breast Soft Tissue Thyroid
#    Non Metastatic             0            0    0      0           0       0
#    Metastatic                14           12    7      7           6       6
#                  Cancer Type
# meta type to all5 Kidney Liver Peripheral Nervous System Biliary Tract
#    Non Metastatic      2     0                         1             0
#    Metastatic          3     4                         3             3
#                  Cancer Type
# meta type to all5 Ampulla of Vater Prostate
#    Non Metastatic                0        0
#    Metastatic                    1        1