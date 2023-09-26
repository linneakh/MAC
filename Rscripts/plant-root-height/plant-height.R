##Linnea Honeker
##plot plant height
##9/26/23

library(tidyverse)
library(dplyr)
library(ggplot2)


list_of_colors <- c("hotpink3", "lightgreen", "cyan3")



###inport data
df.0 = read.delim("./Data/plant-root-height/plantheights_corrected20211209.txt", header = TRUE, sep=",")

df <- df.0 %>%
  select(sprout_name,plot_number, height, prod_des, c_date) %>%
  separate(prod_des, c("line", "extra"), sep = " ") %>%
  select(-extra) %>% #remove outliers
  filter(height<250) %>%
  filter(height>25) %>%
  group_by(plot_number, c_date, line) %>%
  summarize(mean = mean(height)) %>%
  filter(c_date != "2021-09-28")

#change factor levels
df$line <- factor(df$line, c("Low", "Medium", "High"))
  

##plot height data for each genotype

p <- ggplot(df, aes(x=c_date, y=mean, color = line)) +
  geom_point() +
  geom_line(aes(group = plot_number)) +
  labs(y="height (cm)") +
  theme_bw() +
  theme(text = element_text(size = 18))
