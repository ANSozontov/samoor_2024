M2 %>% 
    filter(taxa == "Oribatida" & type == "Binary data (Jaccard)") %>% 
    filter(plants.d == "Juncus" | plants.d == "Convolvulus") %>% 
ggplot(aes(x = axis1, y = axis2, color = plants.d)) + 
    geom_point() + 
    stat_ellipse() 
    # geom_text(aes(label = axis1, x = 0, y = -0.77), color = "black", data = eig, alpha = 0.68) +
    # geom_text(aes(label = axis2, x = -0.99, y = 0), color = "black", data = eig, alpha = 0.68, angle = 90) +
    # facet_grid(cols = vars(type), rows = vars(taxa)) + 
    labs(x = NULL, y = NULL, color = NULL, #subtitle = "Б. Доминантные виды растений") + 
         subtitle = "A. Dominant plant species") + 
    theme(legend.text = element_text(face = "italic"))
