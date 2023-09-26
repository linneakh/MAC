### FUNCTIONS FOR PCA

make_screeplot <- function(eigen, dimensions){
    scree <- ggplot(data=eigen, aes(x=as.factor(dimensions), y=round(variance.percent,digits=2))) +
        geom_bar(stat="identity", fill='steelblue', color="black") +
        theme_linedraw(base_size = size) +
        xlab('Dimensions') +
        ylab('Explained variance (%)') + ylim(0,max(eigen$variance.percent)*1.3) +
        geom_text(data=eigen, aes(label=paste0(round(variance.percent,digits=2),'%')), size=4, vjust=-0.7, hjust=-0.1, angle=30) +
        geom_point(data=eigen, aes(x=as.factor(dimensions), y=variance.percent), color='black') +
        geom_path(data=eigen, aes(x=as.numeric(dimensions), y=variance.percent), color='black') +
        theme(axis.title.x = element_text(size=size+2,face="bold"),
              axis.title.y = element_text(size=size+2,face="bold"),
              axis.text = element_text(size=size-2),
              plot.title = element_text(size=size+2,face="bold"))
    return(scree)
}

make_cumvar <- function(eigen, dimensions){
    cumvar <- ggplot(data=eigen, aes(x=as.factor(dimensions), y=cumulative.variance.percent/100, group=1)) +
        geom_line(size=1.2, color='black') +
        geom_point(size=3, color='black') +
        xlab('PC axes') + ylab('Amount of explained variance') +
        theme_linedraw(base_size = size) + ylim(0, 1.05) +
        theme(axis.title.x = element_text(size=size+2,face="bold"),
              axis.title.y = element_text(size=size+2,face="bold"),
              axis.text = element_text(size=size-2),
              plot.title = element_text(size=size+2,face="bold"))
    return(cumvar)
}

make_pca_plot <- function(pca_coordinates){
    pca_plot <-  ggplot(mapping = aes(x, y)) +
        geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Treatment, shape=Depth), 
                   size=size/3, show.legend = TRUE) +
        # stat_conf_ellipse(data=pca_coordinates, aes(x=PC1, y=PC2, color= XX), 
        #                   alpha=0, geom='polygon', show.legend = FALSE) +
        theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
        scale_color_manual(values= list_of_colors) +
        scale_shape_manual(values= list_of_shapes) +
        theme( legend.text = element_text(face="bold"),
               legend.title = element_blank(),
               legend.key.size = unit(0.6, "cm"),
               legend.key.width = unit(0.6,"cm"),
               axis.title.x = element_text(size=size+2,face="bold"),
               axis.title.y = element_text(size=size+2,face="bold"),
               plot.title = element_text(size=size+2,face="bold")) 
    return(pca_plot)
}

get_arrows <- function(pca, pca_coordinates){
    # prepare arrows
    # extract variable coordinates
    vars <- facto_summarize(pca, element = "var", result=c("coord","contrib","cos2"), axes=c(1,2))
    colnames(vars)[2:3] <-  c("xend", "yend")
    # rescale variable coordinates
    r <- min( (max(pca_coordinates[,"PC1"])-min(pca_coordinates[,"PC1"])/(max(vars[,"xend"])-min(vars[,"xend"]))),
              (max(pca_coordinates[,"PC2"])-min(pca_coordinates[,"PC2"])/(max(vars[,"yend"])-min(vars[,"yend"]))) )
    # multiply vars data by r before biplotting
    # in the original code (https://github.com/kassambara/factoextra/blob/master/R/fviz_pca.R), it calculates r for rescalling and
    # multiplies by 0.7 for plotting... IDK why...
    vars[,c("xend", "yend")] <- vars[,c("xend", "yend")]*r
    arrows <- as.data.frame(vars)
    return(arrows) # columns names, xend, yend for biplot :)
}

make_pca_biplot <- function(pca_coordinates, arrows){
    pca_biplot <- make_pca_plot(pca_coordinates) +
        new_scale_color() +
        geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend, color=contrib),
                     arrow=arrow(length = unit(0.1,"cm")), size=0.6, alpha = 0.5) +
        scale_color_gradient(low="#D9F0A3",high="#005A32") +
        geom_text_repel(data=arrows, aes(x=xend, y=yend, color=contrib),
                        label=arrows$name, size=size/4,show.legend = FALSE) +
        guides(color = guide_colourbar(barwidth = 1, barheight = 5))
    return(pca_biplot)
}
