##Linnea Honeker
##PCA on MNR data from root rhizosphere paper
##1/23/23

library(ggfortify)
library(tidyverse)
library(factoextra)
library(ggnewscale)
library(ggrepel)
library(viridis)
library(vegan)
library(nlme)
library(ecodist)


source('./RScripts/extra_functions.R')



list_of_shapes <- c(15,16,17,18,3,7,8,9,10,4,11,12,13,14)
list_of_colors <- c("blue", "red", "green", "orange")


# parameters for plots
h = 6
w = 7
res = 300
size = 10

###inport data
#NMR root
df.r0 = read.csv("./Data/NMR/MAC_root_NMR.csv", header = TRUE)
colnames(df.r0) <- paste(colnames(df.r0),"root",sep="_")
df.r0b <- df.r0 %>%
  #select(-X_root) %>%
  replace(is.na(.), 0)

rownames(df.r0b) <- df.r0b$X_root
df.r0b$X_root <- NULL


df.rt <- t(df.r0b)

df.r1 <- df.rt +1
df.r <- log(df.r1)
df.r <- as.data.frame(df.r)


#NMR soil
df.s0 = read.csv("./Data/NMR/MAC_soil_NMR.csv", header = TRUE)
colnames(df.s0) <- paste(colnames(df.s0),"soil",sep="_")
df.s0b <- df.s0 %>%
  replace(is.na(.), 0)

rownames(df.s0b) <- df.s0b$X_soil
df.s0b$X_soil <- NULL


df.st <- t(df.s0b)

df.s1 <- df.st +1
df.s <- log(df.s1)
df.s <- as.data.frame(df.s)

#df.s <- df.s %>%
#  select(-V44) #remove because only 2 positive values, and V44 random column

#combine soil and root
df.s.long <- df.s %>%
  rownames_to_column(var = "SampleID") %>%
  gather(key = "metabolite", value = "concentration", -SampleID )

df.r.long <- df.r %>%
  rownames_to_column(var = "SampleID") %>%
  gather(key = "metabolite", value = "concentration", -SampleID )

df.long <- rbind(df.s.long, df.r.long)



df <- df.long %>%
  spread(key = "metabolite", value = "concentration") %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(var = "SampleID")

df.delta.0 <- df.long %>%
  separate(., SampleID, into = c('MAC','Site',	'Treatment',	'Replicate', 'Date', 'Compartment'), sep = "_", remove = T) %>%
  mutate(Line = case_when(
    startsWith(Site, "P8") ~ "High",
    startsWith(Site, "P11") ~  "Low",
    startsWith(Site, "P4") ~ "Low",
    startsWith(Site, "P6") ~ "High",
    startsWith(Site, "P7") ~ "Low",
    startsWith(Site, "P10") ~ "High"
  )) %>%
  mutate(Date = case_when(
    startsWith(Date, "20210718") ~ "July",
    startsWith(Date, "20210909") ~ "Sept",
    startsWith(Date, "20210919") ~ "Sept"
  )) %>%
  spread(key= "metabolite", value = "concentration")

df.r.v0 <- read.csv("./Data/NMR/root_NMR_RVI.csv", header = TRUE) 
df.r.v <- df.r.v0 %>%
  select(-Name, -Formula) %>%
  na.omit() %>%
  arrange(Volatility) %>%
  filter (Category != "high")

df.s.v0 <- read.csv("./Data/NMR/soil_NMR_RVI.csv", header = TRUE) 
df.s.v <- df.s.v0 %>%
  select(-Name, -Formula) %>%
  na.omit()%>%
  filter (Category != "<NA>") %>%
  arrange(Volatility) %>%
  filter (Category != "high")

####NMR with volatiltiy
##outliers to remove:
r="P4_A_Low_S3_Sept_root"
s="P11_A_Low_S3_July_soil"
#soil
df.s.v.final <- df.s.long %>%
  mutate(name = metabolite) %>%
  separate(., SampleID, into = c('MAC','Site',	'Treatment',	'Replicate', 'Date', 'Compartment'), sep = "_", remove = T) %>%
  mutate(Line = case_when(
    startsWith(Site, "P8") ~ "High",
    startsWith(Site, "P11") ~  "Low",
    startsWith(Site, "P4") ~ "Low",
    startsWith(Site, "P6") ~ "High",
    startsWith(Site, "P7") ~ "Low",
    startsWith(Site, "P10") ~ "High"
  )) %>%
  mutate(Date = case_when(
    startsWith(Date, "20210718") ~ "July",
    startsWith(Date, "20210909") ~ "Sept",
    startsWith(Date, "20210919") ~ "Sept"
  )) %>%
  select(-metabolite) %>%
  merge(df.s.v) %>%
  unite(sample, c("Site", "Treatment", "Line", "Replicate", "Date", "Compartment")) %>%
  filter(sample != s) %>%
  separate(sample, c("Site", "Treatment", "Line", "Replicate", "Date", "Compartment")) %>%
  drop_na()

#df.s.v.names <- df.s.v.final %>%
#  group_by(name) %>%
#  summarize_all( mean) %>%
#  arrange(Volatility)

#order.s <- df.s.v.names$name #create list of compounds in order of volatility (soil)

df.s.v.avg <- df.s.v.final %>%
  group_by(Line, Date, name, Category) %>%
  summarize_all(mean) %>%
  select(-MAC, -Site, -Treatment, -Replicate, -Compartment, -Compounds) %>%
  arrange(Volatility)

#convert unsummarized long form back to short form that can be analyzed with PCA, then convert to presence absence

df.s.v.short <- df.s.v.final %>%
  unite(sample, c("Site", "Treatment", "Line", "Replicate", "Date", "Compartment"), sep = "_") %>%
  select(-MAC, -Volatility, -Compounds, -Category) %>%
  spread(key = name, value = concentration) %>%
  column_to_rownames(var = "sample")

######root
df.r.v.final <- df.r.long %>%
  mutate(name = metabolite) %>%
  separate(., SampleID, into = c('MAC','Site',	'Treatment',	'Replicate', 'Date', 'Compartment'), sep = "_", remove = T) %>%
  mutate(Line = case_when(
    startsWith(Site, "P8") ~ "High",
    startsWith(Site, "P11") ~  "Low",
    startsWith(Site, "P4") ~ "Low",
    startsWith(Site, "P6") ~ "High",
    startsWith(Site, "P7") ~ "Low",
    startsWith(Site, "P10") ~ "High"
  )) %>%
  mutate(Date = case_when(
    startsWith(Date, "20210718") ~ "July",
    startsWith(Date, "20210909") ~ "Sept",
    startsWith(Date, "20210919") ~ "Sept"
  )) %>%
  select(-metabolite) %>%
  merge(df.r.v) %>%
  unite(sample, c("Site", "Treatment", "Line", 
                  "Replicate", "Date", "Compartment")) %>%
  filter(sample != r) %>%
  separate(sample, c("Site", "Treatment", "Line", 
                     "Replicate", "Date", "Compartment")) %>%
  drop_na()


#df.r.v.names <- df.r.v.final %>%
#  group_by(name) %>%
#  summarize_all(mean_) %>%
#  arrange(Volatility) 

#order.r <- df.r.v.names$name #create list of compounds in order of volatility (soil)

df.r.v.avg <- df.r.v.final %>%
  group_by(Line, Date, name, Category) %>%
  summarize_all(mean) %>%
  select(-MAC, -Site, -Treatment, -Replicate, -Compartment, -Compounds) %>%
  arrange(Volatility)

#convert unsummarized long form back to short form that can be analyzed with PCA, then convert to presence absence

df.r.v.short <- df.r.v.final %>%
  unite(sample, c("Site", "Treatment", "Line", "Replicate", "Date", "Compartment"), sep = "_") %>%
  select(-MAC, -Volatility, -Compounds, -Category) %>%
  spread(key = name, value = concentration) %>%
  column_to_rownames(var = "sample")
  
###combine filtered root and soil tables (without high volatility)
#summarized table of combined root and soil tables
df.v.avg <- rbind(df.r.v.avg, df.s.v.avg)

#combined root and soil final tables
df.v.final <- rbind(df.r.v.final, df.s.v.final)

#convert unsummarized long form back to short form that can be analyzed with PCA, then convert to presence absence
df.v.short <- df.v.final %>%
  unite(sample, c("Site", "Treatment", "Line", "Replicate", "Date", "Compartment"), sep = "_") %>%
  select(-MAC, -Volatility, -Compounds, -Category) %>%
  spread(key = name, value = concentration) %>%
  column_to_rownames(var = "sample") %>%
  replace(is.na(.), 0)

df.v.short.l <- log(df.v.short +1) 



########################################## PCA on root  NMR#####################
# Calculate PCA with prcomp()
pca <- prcomp(df.r.v.short, scale = TRUE, center = TRUE)
# Extract eigenvalues and variances
eigen <- get_eigenvalue(pca) # this function is from factoextra
dimensions <- c(1:dim(eigen)[1]) # this is for the plot

scree <- make_screeplot(eigen, dimensions) + ggtitle('Scree plot, PCA Class')
filename <- paste0("./Figures/NMR/PCA_screeplot_root.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,scree)
cumvar <-  make_cumvar(eigen, dimensions) + ggtitle('Cumulative variance plot, PCA Class')
filename <- paste0("./Figures/NMR/PCA_cumulative_var_root.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,cumvar)

# extract coordinates for PC1 and PC2
pca_results <- get_pca_ind(pca) #pca[["x"]]
pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2)])
colnames(pca_coordinates) <- c('PC1','PC2')

# merge metadata with pca_coordinates:
pca_coordinates$SampleID <- rownames(pca_coordinates)
pca_coordinates <- pca_coordinates %>% 
  separate(., SampleID, into = c('Plant', 'Treatment', 'Line', 'Replicate', 'Date', 'Compartment'), sep = "_", remove = T)

pca_coordinates$Line <- factor(pca_coordinates$Line, levels = c("Low", "High"))

write.csv(pca_coordinates,file=paste0("./Output/NMR/PCA_individual_coordinates_root.csv"),row.names=TRUE)

# prepare label for graph
pc1 <- paste0('PC1 (',round(eigen$variance.percent[1],digits=1),'%)')
pc2 <- paste0('PC2 (',round(eigen$variance.percent[2],digits=1),'%)')

# arrows
arrows.r <- get_arrows(pca, pca_coordinates)
arrows.r.v <- arrows.r %>%
  left_join(.,df.r.v) 

write.csv(arrows.r,file=paste0("./Output/NMR/PCA_vector_coordinates_root.csv"), row.names=TRUE)

### pca -treatment color and line shape 
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Line, shape=Date), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  scale_color_manual(values = c("pink", "lightblue")) +
  scale_shape_manual(values= list_of_shapes) +
  
  theme( legend.text = element_text(size=size+7, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         panel.grid = element_blank(),
         axis.title.x = element_text(size=size+7,face="bold"),
         axis.title.y = element_text(size=size+7,face="bold"),
         plot.title = element_text(size=size+7,face="bold")) 
pca_plot
filename <- paste0("./Figures/NMR/NMR_PCA_Line_Date_root.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows.r.v, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, show.legend=FALSE) +
  #scale_color_manual(values = c("#EA87F5", "#F3BCF9", "#ACAAAD")) +
  geom_text_repel(data=arrows.r.v, aes(x=xend, y=yend), color = "black",
                  label=arrows.r.v$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/NMR/NMR_PCA_line_Date_root_biplot.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)


########################################## PCA on soil NMR#####################
# Calculate PCA with prcomp()
pca <- prcomp(df.s.v.short, scale = TRUE, center = TRUE)
# Extract eigenvalues and variances
eigen <- get_eigenvalue(pca) # this function is from factoextra
dimensions <- c(1:dim(eigen)[1]) # this is for the plot

scree <- make_screeplot(eigen, dimensions) + ggtitle('Scree plot, PCA Class')
filename <- paste0("./Figures/NMR/PCA_screeplot_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,scree)
cumvar <-  make_cumvar(eigen, dimensions) + ggtitle('Cumulative variance plot, PCA Class')
filename <- paste0("./Figures/NMR/PCA_cumulative_var_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,cumvar)

# extract coordinates for PC1 and PC2
pca_results <- get_pca_ind(pca) #pca[["x"]]
pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2)])
colnames(pca_coordinates) <- c('PC1','PC2')

# merge metadata with pca_coordinates:
pca_coordinates$SampleID <- rownames(pca_coordinates)
pca_coordinates <- pca_coordinates %>% 
  separate(., SampleID, into = c('Plant', 'Treatment', 'Line', 'Replicate', 'Date', 'Compartment'), sep = "_", remove = T)

pca_coordinates$Line <- factor(pca_coordinates$Line, levels = c("Low", "High"))

write.csv(pca_coordinates,file=paste0("./Output/NMR/PCA_individual_coordinates_soil.csv"),row.names=TRUE)

# prepare label for graph
pc1 <- paste0('PC1 (',round(eigen$variance.percent[1],digits=1),'%)')
pc2 <- paste0('PC2 (',round(eigen$variance.percent[2],digits=1),'%)')

# arrows
arrows.s <- get_arrows(pca, pca_coordinates)
arrows.s.v <- arrows.s %>%
  left_join(.,df.s.v) 

write.csv(arrows.s,file=paste0("./Output/NMR/PCA_vector_coordinates_soil.csv"), row.names=TRUE)

### pca -treatment color and line shape 
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Line, shape=Date), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  scale_color_manual(values = c("pink", "lightblue")) +
  scale_shape_manual(values= list_of_shapes) +
  
  theme( legend.text = element_text(size=size+7, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         panel.grid = element_blank(),
         axis.title.x = element_text(size=size+7,face="bold"),
         axis.title.y = element_text(size=size+7,face="bold"),
         plot.title = element_text(size=size+7,face="bold")) 
pca_plot
filename <- paste0("./Figures/NMR/NMR_PCA_Line_Date_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows.s, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, show.legend=FALSE) +
  #scale_color_manual(values = c("#EA87F5", "#F3BCF9", "#ACAAAD")) +
  geom_text_repel(data=arrows.s, aes(x=xend, y=yend), color = "black",
                  label=arrows.s$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/NMR/NMR_PCA_line_Date_biplot-soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)
### pca - date and treatment
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Treatment, shape=Date), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  scale_shape_manual(values= list_of_shapes) +
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         panel.grid = element_blank(),
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
pca_plot
filename <- paste0("./Figures/NMR/NMR_PCA_date_treatment_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows.s, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows.s, aes(x=xend, y=yend), color = "black",
                  label=arrows.s$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/NMR_PCA_date_treatment_biplot_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)


########################################## PCA on combined NMR#####################
# Calculate PCA with prcomp()
pca <- prcomp(df.v.short.l, scale = TRUE, center = TRUE)
# Extract eigenvalues and variances
eigen <- get_eigenvalue(pca) # this function is from factoextra
dimensions <- c(1:dim(eigen)[1]) # this is for the plot

scree <- make_screeplot(eigen, dimensions) + ggtitle('Scree plot, PCA Class')
filename <- paste0("./Figures/NMR/PCA_screeplot_comb.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,scree)
cumvar <-  make_cumvar(eigen, dimensions) + ggtitle('Cumulative variance plot, PCA Class')
filename <- paste0("./Figures/NMR/PCA_cumulative_var_comb.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,cumvar)

# extract coordinates for PC1 and PC2
pca_results <- get_pca_ind(pca) #pca[["x"]]
pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2)])
colnames(pca_coordinates) <- c('PC1','PC2')

# merge metadata with pca_coordinates:
pca_coordinates$SampleID <- rownames(pca_coordinates)
pca_coordinates <- pca_coordinates %>% 
  separate(., SampleID, into = c('Plant', 'Treatment', 'Line', 'Replicate', 'Date', 'Compartment'), sep = "_", remove = T)


write.csv(pca_coordinates,file=paste0("./Output/NMR/PCA_individual_coordinates_comb.csv"),row.names=TRUE)

# prepare label for graph
pc1 <- paste0('PC1 (',round(eigen$variance.percent[1],digits=1),'%)')
pc2 <- paste0('PC2 (',round(eigen$variance.percent[2],digits=1),'%)')

# arrows
arrows <- get_arrows(pca, pca_coordinates)
write.csv(arrows,file=paste0("./Output/NMR/PCA_vector_coordinates_comb.csv"), row.names=TRUE)

### pca -compartment 
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Compartment, shape=Date), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  scale_shape_manual(values= list_of_shapes) +
  
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         panel.grid = element_blank(),
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
pca_plot
filename <- paste0("./Figures/NMR/NMR_PCA_compartment_comb.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows, aes(x=xend, y=yend), color = "black",
                  label=arrows$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/NMR/NMR_PCA_compartment_biplot_scomb.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)


pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows, aes(x=xend, y=yend), color = "black",
                  label=arrows$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/NMR/NMR_PCA_date_compartment_biplot_scomb.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)




###heatmap####
df.s.v.avg$Line <- factor(df.s.v.avg$Line, levels = c("Low", "High"))
df.r.v.avg$Line <- factor(df.r.v.avg$Line, levels = c("Low", "High"))

heatmap_soil <- 
  ggplot(df.s.v.avg, aes(x=Line, y=name, fill=concentration)) + 
  geom_tile()+
  #scale_y_discrete(limits=order.s)+
  #geom_text(aes(label=round(mean)),size=2) +
  facet_wrap(~Date)+
  #scale_fill_gradientn(colours = colorspace::heat_hcl(20))+
  scale_fill_continuous(low = "#fee0d2", high = "purple")+
  labs(x="",y="",fill="log(concentration)")+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            colour = "black", fill = NA, inherit.aes = FALSE) +
  theme(strip.text.y=element_blank(),
        strip.placement="outside",strip.background=element_rect(fill="white"),axis.text.x = element_text(angle = 15, hjust=1),text=element_text(size=16),legend.position="left")
heatmap_soil
filename <- paste0("./Figures/heatmap_soil_log_conc.png")
ggsave(filename,units=c('in'),width=7,height=8,dpi=res,heatmap_soil)


heatmap_root <- 
  ggplot(df.r.v.avg, aes(x=Line, y=name, fill=concentration)) + 
  geom_tile()+
  #scale_y_discrete(limits=order.r)+
  #geom_text(aes(label=round(mean)),size=2) +
  facet_wrap(~Date)+
  #scale_fill_gradientn(colours = colorspace::heat_hcl(20))+
  scale_fill_continuous(low = "#fee0d2", high = "purple")+
  labs(x="",y="",fill="log(concentration)")+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            colour = "black", fill = NA, inherit.aes = FALSE) +
  theme(strip.text.y=element_blank(),
        strip.placement="outside",strip.background=element_rect(fill="white"),axis.text.x = element_text(angle = 15, hjust=1),text=element_text(size=16),legend.position="left")
heatmap_root
filename <- paste0("./Figures/heatmap_root_log_conc.png")
ggsave(filename,units=c('in'),width=7,height=8,dpi=res,heatmap_root)

###delta root - soil####
deltamap <- 
  ggplot(avgdeltametabsID, aes(x=mean, y=Metabolite, color=Line)) + 
  geom_point(size=2,position=position_dodge(width=0.8))+
  geom_errorbarh(aes(xmin = mean-se,xmax = mean+se),height=0, size=0.5,position=position_dodge(width=0.8))+
  scale_y_discrete(limits=rev)+
  geom_hline(yintercept = seq(0, length(unique(avgdeltametabsID$Metabolite))) + .5, color="gray30")+
  geom_vline(xintercept = 0)+
  #coord_trans(x="log2")+
  #scale_shape_manual(values=c(8, 17, 15, 16))+
  #scale_color_brewer(palette = "Set1")+
  #geom_text(aes(label=round(mean)),size=2) +
  facet_grid(~Date)+
  #xlim(
  labs(x="Mean delta (ln(Plant-soil uM))",y="Primary metabolite",fill="Genotype")+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            colour = "black", fill = NA, inherit.aes = FALSE) +
  theme(#strip.text.y=element_blank(),
    strip.placement="outside",strip.background=element_rect(fill="white"),text=element_text(size=16),legend.position="left")

##########NMDS###################
### Calculate nmds for all samples - bray curtis
set.seed(123)
nmds = metaMDS(df.v.short.l, distance = "bray")
plot(nmds)

#pull out envirnmental data
df.stat <- df.v.short.l %>%
  rownames_to_column(var = "sample") %>%
  separate(., sample, into = c("Plant", "Treatment", "Line", "Replicate", "Date", "Compartment"), sep = "_", remove = T) %>%
  select(Plant, Treatment, Line, Replicate, Date, Compartment)
  

#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores(nmds))

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores %>% 
  separate(., SampleID, into = c('Plant', 'Treatment', 'Line', 'Replicate', 'Date', 'Compartment'), sep = "_", remove = T)


write.csv(data.scores,file=paste0("./Output/NMR/nmds_individual_coordinates.csv"),row.names=TRUE)

### pca with no shape 
nmds_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=data.scores, aes(x=NMDS1, y=NMDS2, col=Compartment, shape=Line), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= "NMDS1", y="NMDS2") +
  scale_color_manual(values = list_of_colors) +
  scale_shape_manual(values= list_of_shapes) +
  
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
nmds_plot
filename <- paste0("./Figures/NMR/NMR_nmds_plot_combined.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

#add enviromenmental data as vectors
en = envfit(nmds, df.stat, permutations = 999, na.rm = TRUE)

en_coord_cat = as.data.frame(scores(en, "factors")) #* ordiArrowMul(en) 
en_coord_cat.f <- en_coord_cat[-c(1:9),]

nmds_fact <- nmds_plot + geom_point(data = en_coord_cat.f, aes(x = NMDS1, y = NMDS2), 
                                    shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cat.f, aes(x = NMDS1, y = NMDS2+0.0001 ), 
            label = row.names(en_coord_cat.f), colour = "navy", fontface = "bold")


nmds_fact
filename <- paste0("./Figures/NMR/metaG_nmds_comp.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_fact)

stressplot(nmds)
nmds$stress

### Calculate nmds for all samples - jaccard
set.seed(123)
nmds = metaMDS(df.v.short.l, distance = "jaccard")
plot(nmds)

#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores(nmds))

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores %>% 
  separate(., SampleID, into = c('Plant', 'Treatment', 'Line', 'Replicate', 'Date', 'Compartment'), sep = "_", remove = T)


write.csv(data.scores,file=paste0("./Output/NMR/nmds_individual_coordinates_jaccard.csv"),row.names=TRUE)

### pca with no shape 
nmds_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=data.scores, aes(x=NMDS1, y=NMDS2, col=Compartment, shape=Line), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= "NMDS1", y="NMDS2") +
  scale_color_manual(values = list_of_colors) +
  scale_shape_manual(values= list_of_shapes) +
  
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
nmds_plot
filename <- paste0("./Figures/NMR/NMR_nmds_plot_combined_jaccard.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

#add enviromenmental data as vectors
en = envfit(nmds, df.stat, permutations = 999, na.rm = TRUE)

en_coord_cat = as.data.frame(scores(en, "factors")) #* ordiArrowMul(en) 
en_coord_cat.f <- en_coord_cat[-c(1:9),]

nmds_fact <- nmds_plot + geom_point(data = en_coord_cat.f, aes(x = NMDS1, y = NMDS2), 
                                    shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cat.f, aes(x = NMDS1, y = NMDS2+0.0001 ), 
            label = row.names(en_coord_cat.f), colour = "navy", fontface = "bold")


nmds_fact
filename <- paste0("./Figures/NMR/metaG_nmds_comp_jaccard.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_fact)

stressplot(nmds)
nmds$stress

### Calculate nmds for soil samples - bray 

set.seed(123)
nmds = metaMDS(df.s.v.short, distance = "bray")
plot(nmds)

#pull out envirnmental data
df.stat.s <- df.s.v.short %>%
  rownames_to_column(var = "sample") %>%
  separate(., sample, into = c("Plant", "Treatment", "Line", "Replicate", "Date", "Compartment"), sep = "_", remove = T) %>%
  select(Plant, Treatment, Line, Replicate, Date, Compartment)


#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores(nmds))

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores %>% 
  separate(., SampleID, into = c('Plant', 'Treatment', 'Line', 'Replicate', 'Date', 'Compartment'), sep = "_", remove = T)


write.csv(data.scores,file=paste0("./Output/NMR/nmds_individual_coordinates_soil.csv"),row.names=TRUE)

### pca with no shape 
nmds_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=data.scores, aes(x=NMDS1, y=NMDS2, col=Line, shape=Date), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= "NMDS1", y="NMDS2") +
  scale_color_manual(values = list_of_colors) +
  scale_shape_manual(values= list_of_shapes) +
  
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
nmds_plot
filename <- paste0("./Figures/NMR/NMR_nmds_plot_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

#add enviromenmental data as vectors
en = envfit(nmds, df.stat.s, permutations = 999, na.rm = TRUE)

en_coord_cat = as.data.frame(scores(en, "factors")) #* ordiArrowMul(en) 
en_coord_cat.f <- en_coord_cat[-c(1:9, 12:14, 17),]

nmds_fact <- nmds_plot + geom_point(data = en_coord_cat.f, aes(x = NMDS1, y = NMDS2), 
                                    shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cat.f, aes(x = NMDS1, y = NMDS2+0.0001 ), 
            label = row.names(en_coord_cat.f), colour = "navy", fontface = "bold")


nmds_fact
filename <- paste0("./Figures/NMR/metaG_nmds_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_fact)

stressplot(nmds)
nmds$stress

### Calculate nmds for soil samples - jaccard
set.seed(123)
nmds = metaMDS(df.s.v.short, distance = "jaccard")
plot(nmds)


#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores(nmds))

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores %>% 
  separate(., SampleID, into = c('Plant', 'Treatment', 'Line', 'Replicate', 'Date', 'Compartment'), sep = "_", remove = T)


write.csv(data.scores,file=paste0("./Output/NMR/nmds_individual_coordinates_soil_jaccard.csv"),row.names=TRUE)

### pca with no shape 
nmds_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=data.scores, aes(x=NMDS1, y=NMDS2, col=Line, shape=Date), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= "NMDS1", y="NMDS2") +
  scale_color_manual(values = list_of_colors) +
  scale_shape_manual(values= list_of_shapes) +
  
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
nmds_plot
filename <- paste0("./Figures/NMR/NMR_nmds_plot_soil_jaccard.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

#add enviromenmental data as vectors
en = envfit(nmds, df.stat.s, permutations = 999, na.rm = TRUE)

en_coord_cat = as.data.frame(scores(en, "factors")) #* ordiArrowMul(en) 
en_coord_cat.f <- en_coord_cat[-c(1:9,12:14,17),]

nmds_fact <- nmds_plot + geom_point(data = en_coord_cat.f, aes(x = NMDS1, y = NMDS2), 
                                    shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cat.f, aes(x = NMDS1, y = NMDS2+0.0001 ), 
            label = row.names(en_coord_cat.f), colour = "navy", fontface = "bold")


nmds_fact
filename <- paste0("./Figures/NMR/metaG_nmds_soil_jaccard.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_fact)

stressplot(nmds)
nmds$stress

### Calculate nmds for root samples - bray 
set.seed(123)
nmds = metaMDS(df.r.v.short, distance = "bray")
plot(nmds)

#pull out envirnmental data
df.stat.r <- df.r.v.short %>%
  rownames_to_column(var = "sample") %>%
  separate(., sample, into = c("Plant", "Treatment", "Line", "Replicate", "Date", "Compartment"), sep = "_", remove = T) %>%
  select(Plant, Treatment, Line, Replicate, Date, Compartment)


#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores(nmds))

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores %>% 
  separate(., SampleID, into = c('Plant', 'Treatment', 'Line', 'Replicate', 'Date', 'Compartment'), sep = "_", remove = T)


write.csv(data.scores,file=paste0("./Output/NMR/nmds_individual_coordinates_root.csv"),row.names=TRUE)

### pca with no shape 
nmds_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=data.scores, aes(x=NMDS1, y=NMDS2, col=Line, shape=Date), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= "NMDS1", y="NMDS2") +
  scale_color_manual(values = list_of_colors) +
  scale_shape_manual(values= list_of_shapes) +
  
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
nmds_plot
filename <- paste0("./Figures/NMR/NMR_nmds_plot_root.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

#add enviromenmental data as vectors
en = envfit(nmds, df.stat.r, permutations = 999, na.rm = TRUE)

en_coord_cat = as.data.frame(scores(en, "factors")) #* ordiArrowMul(en) 
en_coord_cat.f <- en_coord_cat[-c(1:9, 12:14,17),]

nmds_fact <- nmds_plot + geom_point(data = en_coord_cat.f, aes(x = NMDS1, y = NMDS2), 
                                    shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cat.f, aes(x = NMDS1, y = NMDS2+0.0001 ), 
            label = row.names(en_coord_cat.f), colour = "navy", fontface = "bold")


nmds_fact
filename <- paste0("./Figures/NMR/metaG_nmds_root.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_fact)

stressplot(nmds)
nmds$stress

### Calculate nmds for root samples - jaccard
set.seed(123)
nmds = metaMDS(df.r.v.short, distance = "jaccard")
plot(nmds)

#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores(nmds))

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores %>% 
  separate(., SampleID, into = c('Plant', 'Treatment', 'Line', 'Replicate', 'Date', 'Compartment'), sep = "_", remove = T)


write.csv(data.scores,file=paste0("./Output/NMR/nmds_individual_coordinates_root_jaccard.csv"),row.names=TRUE)

### pca with no shape 
nmds_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=data.scores, aes(x=NMDS1, y=NMDS2, col=Line, shape=Date), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= "NMDS1", y="NMDS2") +
  scale_color_manual(values = list_of_colors) +
  scale_shape_manual(values= list_of_shapes) +
  
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
nmds_plot
filename <- paste0("./Figures/NMR/NMR_nmds_plot_root_jaccard.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

#add enviromenmental data as vectors
en = envfit(nmds, df.stat.r, permutations = 999, na.rm = TRUE)

en_coord_cat = as.data.frame(scores(en, "factors")) #* ordiArrowMul(en) 
en_coord_cat.f <- en_coord_cat[-c(1:9, 12:14, 17),]

nmds_fact <- nmds_plot + geom_point(data = en_coord_cat.f, aes(x = NMDS1, y = NMDS2), 
                                    shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cat.f, aes(x = NMDS1, y = NMDS2+0.0001 ), 
            label = row.names(en_coord_cat.f), colour = "navy", fontface = "bold")


nmds_fact
filename <- paste0("./Figures/NMR/metaG_nmds_root_jaccard.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_fact)

stressplot(nmds)
nmds$stress


####statistical tests#####
##NMR###
#add columns for compound count and total metabolites as well as categorical variables
#use filtered data without high volatile compounds
df.df <- as.data.frame(df.v.short.l) %>%
  replace(is.na(.), 0)

df.stat <- df.df %>%
  rownames_to_column(var = "SampleID") %>%
  separate(., SampleID, into = c('Plant', 'Treatment', 'Line', 'Replicate', 'Date', 'Compartment'), sep = "_") 

  
# separate into soil and root
df.stat.s <- df.stat %>%
  filter(Compartment == "soil") %>%
  unite(sample, c("Plant", "Treatment", "Line", "Replicate", "Date", "Compartment"), sep = "_") %>%
  column_to_rownames(var = "sample") %>%
  select_if(colSums(.) != 0)

df.s <- df.stat.s %>%
  rownames_to_column(var = "SampleID") %>%
  separate(., SampleID, into = c('Plant', 'Treatment', 'Line', 'Replicate', 'Date', 'Compartment'), sep = "_") 

df.stat.r <- df.stat %>%
  filter(Compartment == "root") %>%
  unite(sample, c("Plant", "Treatment", "Line", "Replicate", "Date", "Compartment"), sep = "_") %>%
  column_to_rownames(var = "sample") %>%
  select_if(colSums(.) != 0)

df.r <- df.stat.r %>%
  rownames_to_column(var = "SampleID") %>%
  separate(., SampleID, into = c('Plant', 'Treatment', 'Line', 'Replicate', 'Date', 'Compartment'), sep = "_") 

####permanova on bray-curtis distance matrix####
##root and soil##
dist<- distance(df.df, "bray-curtis")

adonis(dist~Compartment +  Date + Line + Treatment, df.stat)
#.              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Compartment  1    5.8937  5.8937 297.548 0.82144  0.001 ***
#  Date         1    0.0100  0.0100   0.507 0.00140  0.535    
#Line         1    0.0504  0.0504   2.543 0.00702  0.101    
#Treatment    2    0.0719  0.0359   1.814 0.01002  0.162    
#Residuals   58    1.1488  0.0198         0.16012           
#Total       63    7.1749                 1.00000         

##root and soil##
dist<- distance(df.df, "jaccard")

adonis(dist~Compartment +  Date + Line + Treatment , df.stat)
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Compartment  1    6.8581  6.8581 172.603 0.72936  0.001 ***
#  Date         1    0.0481  0.0481   1.210 0.00511  0.258    
#Line         1    0.0552  0.0552   1.390 0.00587  0.221    
#Treatment    2    0.1369  0.0684   1.723 0.01456  0.131    
#Residuals   58    2.3045  0.0397         0.24509           
#Total       63    9.4029                 1.00000    

#soil only
dist.soil <- distance(df.stat.s, method = "bray-curtis")

adonis(dist.soil~ Date + Line + Treatment +  Date*Line + Date*Treatment, df.s)
#                Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#Date            1   0.04182 0.041819  1.5839 0.04001  0.167  
#Line            1   0.07567 0.075670  2.8660 0.07240  0.027 *
#  Treatment       2   0.06972 0.034860  1.3203 0.06671  0.247  
#Date:Line       1   0.06707 0.067074  2.5405 0.06418  0.044 *
#  Date:Treatment  2   0.07802 0.039009  1.4775 0.07465  0.168  
#Residuals      27   0.71286 0.026402         0.68206         
#Total          34   1.04516                  1.00000    

dist.soil <- distance(df.stat.s, method = "jaccard")

adonis(dist.soil~ Date + Line + Treatment +  Date*Line + Date*Treatment, df.s)
#               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#Date            1   0.08232 0.082317  1.4111 0.03751  0.219  
#Line            1   0.06673 0.066729  1.1439 0.03041  0.346  
#Treatment       2   0.13664 0.068321  1.1712 0.06227  0.307  
#Date:Line       1   0.15097 0.150966  2.5880 0.06880  0.027 *
#  Date:Treatment  2   0.18273 0.091364  1.5662 0.08327  0.120  
#Residuals      27   1.57501 0.058334         0.71774         
#Total          34   2.19439                  1.00000          

#root only
dist.root <- distance(df.stat.r, method = "bray-curtis")

adonis(dist.root~ Date + Line + Treatment +  Date*Line + Date*Treatment, df.r)
#.             Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)  
#Date            1  0.010153 0.0101534 1.29477 0.04303  0.275  
#Line            1  0.017373 0.0173735 2.21547 0.07362  0.077 .
#Treatment       2  0.023128 0.0115640 1.47464 0.09801  0.188  
#Date:Line       1  0.006239 0.0062388 0.79558 0.02644  0.568  
#Date:Treatment  2  0.014400 0.0071999 0.91813 0.06102  0.551  
#Residuals      21  0.164680 0.0078419         0.69787         
#Total          28  0.235973                   1.00000   

dist.root <- distance(df.stat.r, method = "jaccard")

adonis(dist.root~ Date + Line + Treatment +  Date*Line + Date*Treatment, df.r)
#               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
#Date            1   0.01579 0.015787  1.4162 0.04506  0.251
#Line            1   0.02162 0.021625  1.9399 0.06172  0.125
#Treatment       2   0.03902 0.019511  1.7503 0.11138  0.104
#Date:Line       1   0.01452 0.014523  1.3028 0.04145  0.294
#Date:Treatment  2   0.02531 0.012657  1.1354 0.07225  0.367
#Residuals      21   0.23410 0.011147         0.66814       
#Total          28   0.35037                  1.00000  


