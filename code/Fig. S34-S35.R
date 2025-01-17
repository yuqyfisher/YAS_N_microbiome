####Fig.S34-S35
#R3.6.2

####Fig.S34
library(ggplot2)
setwd("mine/plot")

phylum1 <- read.delim('data/virus_bacteria.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

p1 <- ggplot(phylum1, aes(x = '', y = proportion, fill = taxa)) +
  geom_bar(stat = 'identity', width = 1) +
  coord_polar(theta = 'y') +
  facet_wrap(~sample,nrow=6) +
  geom_text(aes(label = paste0(proportion, "%")), position = position_stack(vjust = 0.5), size = 4) + 
  scale_fill_manual(values = rev(c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#B83945', '#CFE7C4'))) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(face = 'bold'), legend.title = element_blank()) +
  labs(x = '', y = '', title = 'All samples', fill = 'taxa')

p1
ggsave("Fig.S34.pdf", p1,height=20,width =30,limitsize = FALSE )


###Fig.S35
phylum1 <- read.delim('data/prokaryotic_eukaryotic.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

p2 <- ggplot(phylum1, aes(x = '', y = proportion, fill = taxa)) +
  geom_bar(stat = 'identity', width = 1) +
  coord_polar(theta = 'y') +
  facet_wrap(~sample,nrow=6) +
  geom_text(aes(label = paste0(proportion, "%")), position = position_stack(vjust = 0.5), size = 4) + 
  scale_fill_manual(values = rev(c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#BEBADA', '#E47159'))) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(face = 'bold'), legend.title = element_blank()) +
  labs(x = '', y = '', title = 'All samples', fill = 'taxa')

p2
ggsave("Fig.S35.pdf", p1,height=20,width =30,limitsize = FALSE )

