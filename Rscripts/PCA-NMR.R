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


# parameters for plots
h = 6
w = 7
res = 300
size = 10

###inport data
#NMR root
df.r0 = read.csv("./Data/MAC_root_NMR.csv", header = TRUE)
colnames(df.r0) <- paste(colnames(df.r0),"root",sep="_")
df.r0b <- df.r0 %>%
  select(-X.1_root) %>%
  replace(is.na(.), 0)

rownames(df.r0b) <- df.r0b$X_root
df.r0b$X_root <- NULL


df.rt <- t(df.r0b)

df.r1 <- df.rt +1
df.r <- log(df.r1)
df.r <- as.data.frame(df.r)


#NMR soil
df.s0 = read.csv("./Data/MAC_soil_NMR.csv", header = TRUE)
colnames(df.s0) <- paste(colnames(df.s0),"soil",sep="_")
df.s0b <- df.s0 %>%
  replace(is.na(.), 0)

rownames(df.s0b) <- df.s0b$X_soil
df.s0b$X_soil <- NULL


df.st <- t(df.s0b)

df.s1 <- df.st +1
df.s <- log(df.s1)
df.s <- as.data.frame(df.s)

df.s <- df.s %>%
  select(-V44) #remove because only 2 positive values, and V44 random column

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

df.r.v0 <- read.csv("./Data/root_NMR_RVI.csv", header = TRUE) 
df.r.v <- df.r.v0 %>%
  select(-Name, -Formula) %>%
  na.omit() %>%
  arrange(Volatility)

df.s.v0 <- read.csv("./Data/soil_NMR_RVI.csv", header = TRUE) 
df.s.v <- df.s.v0 %>%
  select(-Name, -Formula) %>%
  na.omit()%>%
  filter (Category != "<NA>") %>%
  arrange(Volatility)

####NMR with volatiltiy
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
  left_join(.,df.s.v) %>%
  na.omit()

df.s.v.names <- df.s.v.final %>%
  group_by(name) %>%
  summarize_all( mean) %>%
  arrange(Volatility)

order.s <- df.s.v.names$name #create list of compounds in order of volatility (soil)

df.s.v.avg <- df.s.v.final %>%
  group_by(Line, Date, name, Category) %>%
  summarize_all(concentration, mean) %>%
  select(-MAC, -Site, -Treatment, -Replicate, -Compartment, -Compounds) %>%
  arrange(Volatility)

#root
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
  left_join(.,df.r.v) %>%
  na.omit()

df.r.v.names <- df.r.v.final %>%
  group_by(name) %>%
  summarize_all(mean) %>%
  arrange(Volatility)

order.r <- df.r.v.names$name #create list of compounds in order of volatility (soil)

df.r.v.avg <- df.r.v.final %>%
  group_by(Line, Date, name, Category) %>%
  summarize_all(mean) %>%
  select(-MAC, -Site, -Treatment, -Replicate, -Compartment, -Compounds) %>%
  arrange(Volatility)

########################################## PCA on root  NMR#####################
# Calculate PCA with prcomp()
pca <- prcomp(df.r, scale = TRUE, center = TRUE)
# Extract eigenvalues and variances
eigen <- get_eigenvalue(pca) # this function is from factoextra
dimensions <- c(1:dim(eigen)[1]) # this is for the plot

scree <- make_screeplot(eigen, dimensions) + ggtitle('Scree plot, PCA Class')
filename <- paste0("./Figures/PCA_screeplot.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,scree)
cumvar <-  make_cumvar(eigen, dimensions) + ggtitle('Cumulative variance plot, PCA Class')
filename <- paste0("./Figures/PCA_cumulative_var.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,cumvar)

# extract coordinates for PC1 and PC2
pca_results <- get_pca_ind(pca) #pca[["x"]]
pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2)])
colnames(pca_coordinates) <- c('PC1','PC2')

# merge metadata with pca_coordinates:
pca_coordinates$SampleID <- rownames(pca_coordinates)
pca_coordinates <- pca_coordinates %>% 
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
  ))

pca_coordinates$Line <- factor(pca_coordinates$Line, levels = c("Low", "High"))

write.csv(pca_coordinates,file=paste0("./Output/PCA_individual_coordinates.csv"),row.names=TRUE)

# prepare label for graph
pc1 <- paste0('PC1 (',round(eigen$variance.percent[1],digits=1),'%)')
pc2 <- paste0('PC2 (',round(eigen$variance.percent[2],digits=1),'%)')

# arrows
arrows.r <- get_arrows(pca, pca_coordinates)
arrows.r.v <- arrows.r %>%
  left_join(.,df.r.v) 

write.csv(arrows.r,file=paste0("./Output/PCA_vector_coordinates.csv"), row.names=TRUE)

### pca -treatment color and line shape 
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Line, shape=Date), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +

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
filename <- paste0("./Figures/NMR_PCA_Line_label.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows.r.v, aes(x=0, y=0, xend=xend, yend=yend, color = Category),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, show.legend=FALSE) +
  scale_color_manual(values = c("#EA87F5", "#F3BCF9", "#ACAAAD")) +
  geom_text_repel(data=arrows.r.v, aes(x=xend, y=yend), color = "black",
                  label=arrows$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/NMR_PCA_line_label_biplot.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)

### pca -Date color and line shape 
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Label, shape=Line), 
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
filename <- paste0("./Figures/NMR_PCA_date_line.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows, aes(x=xend, y=yend), color = "black",
                  label=arrows$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/NMR_PCA_date_line_biplot.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)
########################################## PCA on soil NMR#####################
# Calculate PCA with prcomp()
pca <- prcomp(df.s, scale = TRUE, center = TRUE)
# Extract eigenvalues and variances
eigen <- get_eigenvalue(pca) # this function is from factoextra
dimensions <- c(1:dim(eigen)[1]) # this is for the plot

scree <- make_screeplot(eigen, dimensions) + ggtitle('Scree plot, PCA Class')
filename <- paste0("./Figures/PCA_screeplot_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,scree)
cumvar <-  make_cumvar(eigen, dimensions) + ggtitle('Cumulative variance plot, PCA Class')
filename <- paste0("./Figures/PCA_cumulative_var_sol.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,cumvar)

# extract coordinates for PC1 and PC2
pca_results <- get_pca_ind(pca) #pca[["x"]]
pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2)])
colnames(pca_coordinates) <- c('PC1','PC2')

# merge metadata with pca_coordinates:
pca_coordinates$SampleID <- rownames(pca_coordinates)
pca_coordinates <- pca_coordinates %>% 
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
  ))

pca_coordinates$Line <- factor(pca_coordinates$Line, levels = c("Low", "High"))


write.csv(pca_coordinates,file=paste0("./Output/PCA_individual_coordinates_soil.csv"),row.names=TRUE)

# prepare label for graph
pc1 <- paste0('PC1 (',round(eigen$variance.percent[1],digits=1),'%)')
pc2 <- paste0('PC2 (',round(eigen$variance.percent[2],digits=1),'%)')

# arrows
arrows.s <- get_arrows(pca, pca_coordinates)
arrows.s.v <- arrows.s %>%
  left_join(.,df.s.v)  %>%
  na.omit()

write.csv(arrows.s,file=paste0("./Output/PCA_vector_coordinates_soil.csv"), row.names=TRUE)

### pca -treatment color and line shape 
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Line, shape=Date), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
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
filename <- paste0("./Figures/NMR_PCA_label_line_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows, aes(x=xend, y=yend), color = "black",
                  label=arrows$name, size=size/3, show.legend = FALSE)


pca_biplot 


filename <- paste0("./Figures/NMR_PCA_treatment_line_biplot_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)

### pca - date and treatment
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Label, shape=Treatment), 
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
filename <- paste0("./Figures/NMR_PCA_date_treatment_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows, aes(x=xend, y=yend), color = "black",
                  label=arrows$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/NMR_PCA_date_treatment_biplot_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)

### pca - date and line
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Label, shape=Line), 
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
filename <- paste0("./Figures/NMR_PCA_date_tline_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows, aes(x=xend, y=yend), color = "black",
                  label=arrows$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/NMR_PCA_date_line_biplot_soil.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)

########################################## PCA on combined NMR#####################
# Calculate PCA with prcomp()
pca <- prcomp(df, scale = TRUE, center = TRUE)
# Extract eigenvalues and variances
eigen <- get_eigenvalue(pca) # this function is from factoextra
dimensions <- c(1:dim(eigen)[1]) # this is for the plot

scree <- make_screeplot(eigen, dimensions) + ggtitle('Scree plot, PCA Class')
filename <- paste0("./Figures/PCA_screeplot_comb.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,scree)
cumvar <-  make_cumvar(eigen, dimensions) + ggtitle('Cumulative variance plot, PCA Class')
filename <- paste0("./Figures/PCA_cumulative_var_comb.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,cumvar)

# extract coordinates for PC1 and PC2
pca_results <- get_pca_ind(pca) #pca[["x"]]
pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2)])
colnames(pca_coordinates) <- c('PC1','PC2')

# merge metadata with pca_coordinates:
pca_coordinates$SampleID <- rownames(pca_coordinates)
pca_coordinates <- pca_coordinates %>% 
  separate(., SampleID, into = c('MAC','Site',	'Treatment',	'Replicate', 'Date', 'Compartment'), sep = "_", remove = T) %>%
  mutate(Line = case_when(
    startsWith(Site, "P8") ~ "High",
    startsWith(Site, "P11") ~  "Low",
    startsWith(Site, "P4") ~ "Low",
    startsWith(Site, "P6") ~ "High",
    startsWith(Site, "P7") ~ "Low",
    startsWith(Site, "P10") ~ "High"
  )) %>%
  mutate(Label = case_when(
    startsWith(Date, "20210718") ~ "Label2",
    startsWith(Date, "20210909") ~ "Label1",
    startsWith(Date, "20210919") ~ "Label1"
  ))


write.csv(pca_coordinates,file=paste0("./Output/PCA_individual_coordinates_comb.csv"),row.names=TRUE)

# prepare label for graph
pc1 <- paste0('PC1 (',round(eigen$variance.percent[1],digits=1),'%)')
pc2 <- paste0('PC2 (',round(eigen$variance.percent[2],digits=1),'%)')

# arrows
arrows <- get_arrows(pca, pca_coordinates)
write.csv(arrows,file=paste0("./Output/PCA_vector_coordinates_comb.csv"), row.names=TRUE)

### pca -compartment 
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Compartment), 
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
filename <- paste0("./Figures/NMR_PCA_compartment_comb.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows, aes(x=xend, y=yend), color = "black",
                  label=arrows$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/NMR_PCA_compartment_biplot_scomb.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)


### pca -treatment color and line shape 
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Label, shape=Compartment), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  #scale_color_manual(values = colors_plant_species,
  #                   breaks = c("Clitoria", "Hibiscus", "Piper")) +
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
filename <- paste0("./Figures/NMR_PCA_date_compartment_comb.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows, aes(x=xend, y=yend), color = "black",
                  label=arrows$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/NMR_PCA_date_compartment_biplot_scomb.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)




###heatmap####
df.s.v.avg$Line <- factor(df.s.v.avg$Line, levels = c("Low", "High"))
df.r.v.avg$Line <- factor(df.r.v.avg$Line, levels = c("Low", "High"))

heatmap_soil <- 
  ggplot(df.s.v.avg, aes(x=Line, y=name, fill=concentration)) + 
  geom_tile()+
  scale_y_discrete(limits=order.s)+
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
  scale_y_discrete(limits=order.r)+
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

####statistical tests#####
##NMR###
#add columns for compound count and total metabolites as well as categorical variables
#df.df <- as.data.frame(df)
df.df<- as.data.frame(df)

df.stat <- df.df %>%
  rownames_to_column(var = "SampleID") %>%
  separate(., SampleID, into = c('MAC','Site',	'Treatment',	'Replicate', 'Date', 'Compartment'), sep = "_", remove = T) %>%
  mutate(Line = case_when(
    startsWith(Site, "P8") ~ "High",
    startsWith(Site, "P11") ~  "Low",
    startsWith(Site, "P4") ~ "Low",
    startsWith(Site, "P6") ~ "High",
    startsWith(Site, "P7") ~ "Low",
    startsWith(Site, "P10") ~ "High"
  )) %>%
  mutate(Label = case_when(
    startsWith(Date, "20210718") ~ "Label2",
    startsWith(Date, "20210909") ~ "Label1",
    startsWith(Date, "20210919") ~ "Label1"
  ))

  
# separate into soil and root
df.stat.s <- df.stat %>%
  filter(Compartment == "soil")

df.stat.r <- df.stat %>%
  filter(Compartment == "root")




####t-test####

##soil##
modelList <- list()
for(i in 7:67) {
  fmla <- formula(paste(names(df.stat.s) [i]," ~ Line"))
  modelList[[i]]<- t.test(fmla, data = df.stat.s, paired = FALSE)
}
modelList

sink("./Output/ttest-soil-line.txt")
print(modelList)
sink()

##hibiscus##
modelList <- list()
for(i in 5:51) {
  fmla <- formula(paste(names(df.s.hib) [i]," ~ Condition"))
  modelList[[i]]<- t.test(fmla, data = df.s.hib, paired = FALSE)
}
modelList

sink("./Output/ttest-hibiscus-results-29.txt")
print(modelList)
sink()

##piper##
modelList <- list()
for(i in 5:51) {
  fmla <- formula(paste(names(df.s.pip) [i]," ~ Condition"))
  modelList[[i]]<- t.test(fmla, data = df.s.pip, paired = FALSE)
}
modelList

sink("./Output/ttest-piper-result-29s.txt")
print(modelList)
sink()

#linear mixed effects
modelList <- list()
for(i in 5:49) {
  fmla <- formula(paste(names(df.s) [i]," ~ Condition"))
  lme<- lme(fmla,
            random = list(PlantSpecies = ~1),
            data = df.s,
            #weights =  varIdent(form = ~1|Condition),
            na.action=na.omit)
  modelList[[i]] <- summary(lme)
}
modelList

sink("./Output/NMR-LME-29.txt")
print(modelList)
sink()

####permanova on bray-curtis distance matrix####
##root and soil##
dist<- distance(df, "bray-curtis")

adonis(dist~Compartment + Label + Line + Treatment + Compartment*Label + Compartment*Line + Compartment*Treatment, df.stat)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Compartment            1    5.4718  5.4718 257.575 0.79244  0.001 ***
#  Label                  1    0.0331  0.0331   1.559 0.00480  0.212    
#Line                   1    0.0756  0.0756   3.560 0.01095  0.033 *  
#  Treatment              2    0.0797  0.0399   1.877 0.01155  0.142    
#Compartment:Label      1    0.0469  0.0469   2.208 0.00679  0.113    
#Compartment:Line       1    0.0659  0.0659   3.102 0.00954  0.069 .  
#Compartment:Treatment  2    0.0484  0.0242   1.140 0.00701  0.332    
#Residuals             51    1.0834  0.0212         0.15690           
#Total                 60    6.9049                 1.00000        

# root only
dist.root <- distance(df.r, "bray-curtis")

adonis(dist.root~Label + Line + Treatment + Label*Line + Label*Treatment + Line*Treatment, df.stat.r)
#Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)
#Label            1   0.00610 0.0060988 0.34669 0.01251  0.910
#Line             1   0.01727 0.0172672 0.98157 0.03541  0.390
#Treatment        2   0.02675 0.0133731 0.76021 0.05485  0.624
#Label:Line       1   0.01600 0.0160014 0.90962 0.03282  0.436
#Label:Treatment  2   0.02060 0.0103020 0.58563 0.04226  0.813
#Line:Treatment   2   0.04904 0.0245214 1.39395 0.10058  0.179
#Residuals       20   0.35183 0.0175913         0.72157       
#Total           29   0.48759                   1.00000       

#soil only
dist.soil <- distance(df.s, "bray-curtis")

adonis(dist.soil~Label + Line + Treatment + Label*Line + Label*Treatment + Line*Treatment, df.stat.s)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#Label            1   0.04282 0.042816  1.5411 0.04528  0.138  
#Line             1   0.04907 0.049067  1.7661 0.05189  0.085 .
#Treatment        2   0.06387 0.031937  1.1495 0.06755  0.311  
#Label:Line       1   0.08926 0.089260  3.2128 0.09440  0.013 *
#  Label:Treatment  2   0.06164 0.030820  1.1093 0.06519  0.353  
#Line:Treatment   2   0.05547 0.027734  0.9983 0.05866  0.445  
#Residuals       21   0.58344 0.027783         0.61703         
#Total           30   0.94557                  1.00000             
