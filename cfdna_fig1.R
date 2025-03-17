# Cargar datos
library(ggplot2)
library(ggbreak)
cfDNA <- read.table("~/Downloads/cfdna_data.txt", header = TRUE, sep = "\t")

# Mock dataset with specified ranges
range_data <- data.frame(
  Range = c("20-150 bp", "160-180 bp", "180-220 bp", "250-320 bp", "320-700 bp", "700-1300 bp"),
  Start = c(20, 160, 180, 250, 320, 700),
  End = c(150, 180, 220, 320, 700, 1300),
  Label = c("20-150 bp", "160-180 bp", "180-220 bp", "250-320 bp", "320-700 bp", "700-1300 bp")
)

# Print the dataset
print(range_data)

# Add a midpoint for each range to place the vectors
range_data <- range_data %>%
  mutate(Midpoint = (Start + End) / 2)
range_data$Range <- factor(range_data$Range, levels = c("20-150 bp", "160-180 bp", "180-220 bp", 
                                                        "250-320 bp", "320-700 bp", "700-1300 bp"))                                                        
                                                    
# Updated plot code with elegant segments. Iteration = 42
p4 <- ggplot(data = cfDNA, aes(x = Size, y = density)) +
  geom_line(aes(color = 'cornflowerblue'), size=1.75, show.legend = FALSE) +
  labs(title = "Example Size Profile",
       x = "Fragment size in bp",
       y = "Percentage of fragments") +
  xlim(1, 1312) +
  scale_x_break(c(330, 680), space = 0.25) +
  scale_x_break(c(720, 1275), space = 0.25) +
  geom_vline(xintercept = c(20, 150, 160, 180, 220, 250, 320, 700, 1300), linetype = c(5, 5, 5, 5, 5, 5, 5, 5, 5), col = 'grey') +
  #geom_vline(xintercept = c(346,678.5), linetype = c(1,1), col = 'black') +
  annotate("text", x = c(20, 150, 160, 180, 220, 250, 320, 700, 1300), y = -0.1,
           label = c("20", "150", "\n160", "180", "220", "250", "320", "700", "1300"), size=7) +
  # Add segments for each range
  geom_segment(data = range_data,
               aes(x = Start, xend = End, y = -0.025, yend = -0.025, color = Range),
               size = 4, inherit.aes = FALSE) + # Elegant lines below the plot
  scale_color_manual(values = c("20-150 bp" = "#FF0000", "160-180 bp" = "#00A08A",
                                "180-220 bp" = "#F2AD00", "250-320 bp" = "#5BBCD6",
                                "320-700 bp" = "#F98400", "700-1300 bp" = "#9B59B6")) + # Custom colors for ranges
  theme_linedraw() +
  labs(color = "Ranges") + # Legend for segments
  theme(axis.text.x = element_blank(),        # Remove x-axis labels
        axis.ticks.x = element_blank(),        # Remove x-axis ticks
        axis.text.y = element_text(size = 18),  # Adjust y-axis text size for visibility
        axis.title.y = element_text(size = 24), # Adjust y-axis title size for visibility
        axis.title.x = element_text(size = 24), 
        panel.grid = element_blank(),
        legend.text = element_text(size = 20),    # Increase legend text size
        legend.title = element_text(size = 22),   # Increase legend title size
        legend.key.size = unit(1.5, "cm"))        # Increase legend key size  
p4

# Save the plot with high resolution
ggsave("/home/jlvillanueva/Downloads/cfDNA_plot_high_res.png", width = 16, height = 9, dpi = 300, units = "in") # High-resolution settings
