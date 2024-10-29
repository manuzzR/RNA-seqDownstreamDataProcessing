#Load necessary packages
library(dplyr)
library(ggplot2)

#Read in the files
trinity_cog <- read.table("cog_annotations.txt", sep = "\t", header = TRUE)
fun_20 <- read.table("fun-20.tab",sep = "\t",  header = FALSE) %>% 
  select(V1,V3)

cog_20_def <- read.table("cog-20.def.tab", sep = "\t",  header = FALSE) %>%
  select(V1, V2)

# Join the data
data <- left_join(cog_20_def, fun_20, by = c("V2" = "V1"))
data <- left_join(data, trinity_cog, by = c("V1" = "COG_number"))

# Group the data by COG category and count the number of genes in each category
data_count <- data %>% group_by(V3) %>% summarize(count = n())
final_data <- left_join(data_count, fun_20, by = c("V3" = "V3")) %>% na.omit()
final_data$legend <- paste(final_data$V1,":" , final_data$V3)

# Create the graph using ggplot2

cog_analysis_figure <- ggplot(final_data, aes(x = V1, y = count, fill = legend), plot.width = 12, plot.height = 8) +
  geom_bar(stat = "identity") +
  labs(title= "COG Category Distribution",
       x ="COG Category",
       y = "Number of genes") +
  #ggtitle("COG Category Distribution") +
  #xlab("COG Category") +
  #ylab("Number of Genes") +
  theme (axis.text.x = element_text(hjust = 0.5, vjust = 0.5, face = "bold", size = 12),
        axis.text.y = element_text(face="bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
        plot.title = element_text(face = "bold", size = 20, colour = "black", hjust = 0.5),
    #legend.key.size = unit(0.5, 'cm'),
        #legend.key.height = unit(0.5, 'cm'),
        #legend.key.width = unit(0.5, 'cm'),
        #legend.title = element_text(size=12),
        #legend.text = element_text(size=10),
        legend.box = "vertical") +
  guides(fill = guide_legend(nrow = 36,
                             label.vjust = 0.5,
                             position = "right",
                             direction = "vertical",
                             size = 16))



        
ggsave(paste0("cog_analysis_figure", Sys.Date(), ".jpg"), height = 8, width = 8, units = "in")


