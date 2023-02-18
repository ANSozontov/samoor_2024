# initial -----------------------------------------------------------------
library(tidyverse)
theme_set(theme_bw() + theme(legend.position = "bottom"))
load("17.02.2023.RData")
div <- div %>% 
    # filter(substrate != "debris") %>% 
    select(-Astigmata_total) %>% 
    mutate(
        coast = factor(coast, levels = c("pebbly", "sandy beach", "reeds")), 
        vegetation = case_when(vegetation == "no" ~ "0", TRUE ~ vegetation), 
        vegetation = factor(as.numeric(vegetation), ordered = TRUE)
)

# general -----------------------------------------------------------------
# coast
div %>% 
    pivot_longer(names_to = "taxa", values_to = "abu", -c(1:16)) %>% 
    ggplot(aes(x = coast, y = abu, fill = taxa)) + 
    geom_col(position = "dodge") + 
    labs(x = "Coast", y = "Abundance raw")
div %>% 
    pivot_longer(names_to = "taxa", values_to = "abu", -c(1:16)) %>% 
    ggplot(aes(x = coast, y = abu, fill = taxa)) + 
    geom_col(position = "fill") + 
    labs(x = "Coast", y = "Abundance relative")

# skew
div %>% 
    pivot_longer(names_to = "taxa", values_to = "abu", -c(1:16)) %>% 
    ggplot(aes(x = skew, y = abu, fill = taxa)) + 
    geom_col(position = "dodge") + 
    labs(x = "skew", y = "Abundance raw")
div %>% 
    pivot_longer(names_to = "taxa", values_to = "abu", -c(1:16)) %>% 
    ggplot(aes(x = skew, y = abu, fill = taxa)) + 
    geom_col(position = "fill") + 
    labs(x = "skew", y = "Abundance relative")

# soil
div %>% 
    pivot_longer(names_to = "taxa", values_to = "abu", -c(1:16)) %>% 
    ggplot(aes(x = soil, y = abu, fill = taxa)) + 
    geom_col(position = "dodge") + 
    labs(x = "soil", y = "Abundance raw")
div %>% 
    pivot_longer(names_to = "taxa", values_to = "abu", -c(1:16)) %>% 
    ggplot(aes(x = soil, y = abu, fill = taxa)) + 
    geom_col(position = "fill") + 
    labs(x = "soil", y = "Abundance relative")

# vegetation
div %>% 
    pivot_longer(names_to = "taxa", values_to = "abu", -c(1:16)) %>% 
    ggplot(aes(x = vegetation, y = abu, fill = taxa)) + 
    geom_col(position = "dodge") + 
    labs(x = "Vegetation cover, %", y = "Abundance raw")
div %>% 
    pivot_longer(names_to = "taxa", values_to = "abu", -c(1:16)) %>% 
    ggplot(aes(x = vegetation, y = abu, fill = taxa)) + 
    geom_col(position = "fill") + 
    labs(x = "Vegetation cover, %", y = "Abundance relative")

# Dominant plants
div %>% 
    pivot_longer(names_to = "taxa", values_to = "abu", -c(1:16)) %>% 
    ggplot(aes(x = plants.d, y = abu, fill = taxa)) + 
    geom_col(position = "dodge") + 
    labs(x = "Dominant plants", y = "Abundance raw") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
div %>% 
    pivot_longer(names_to = "taxa", values_to = "abu", -c(1:16)) %>% 
    ggplot(aes(x = plants.d, y = abu, fill = taxa)) + 
    geom_col(position = "fill") + 
    labs(x = "Dominant plants", y = "Abundance relative")+ 
    theme(axis.text.x = element_text(angle = 45, hjust = 0.9))

# Diversity ---------------------------------------------------------------

# Coast
div %>% 
    select(coast:mst_iH) %>% 
    pivot_longer(names_to = "v1", values_to = "vr", -c(1:7)) %>% 
    filter(!is.na(vr)) %>% 
    separate(v1, into = c("taxa", "diversity"), sep = "_", extra = "merge") %>% 
    mutate(diversity2 = case_when(
        diversity == "obs_m"  ~"1. Abundance", 
        diversity == "obs_qD" ~ "2. Observed N of sp.",
        diversity == "iH"     ~ "4. Shannon index", 
        TRUE ~ "3. Rarefied N of sp."
        ), taxa = case_when(taxa == "orb" ~ "Oribatida", TRUE ~ "Mesostigmata")) %>% 
    ggplot(aes(x = coast, y = vr, fill = coast)) + 
    geom_boxplot() + 
    facet_grid(rows = vars(diversity2), cols = vars(taxa), scales = "free") + 
    labs(subtitle = "Oribatida rarefied to 20, Mesostigmata to 10 individuals", 
         y = NULL, x = NULL)

# skew
div %>% 
    select(coast:mst_iH) %>% 
    pivot_longer(names_to = "v1", values_to = "vr", -c(1:7)) %>% 
    filter(!is.na(vr)) %>% 
    separate(v1, into = c("taxa", "diversity"), sep = "_", extra = "merge") %>% 
    mutate(diversity2 = case_when(
        diversity == "obs_m"  ~"1. Abundance", 
        diversity == "obs_qD" ~ "2. Observed N of sp.",
        diversity == "iH"     ~ "4. Shannon index", 
        TRUE ~ "3. Rarefied N of sp."
    ), taxa = case_when(taxa == "orb" ~ "Oribatida", TRUE ~ "Mesostigmata")) %>% 
    ggplot(aes(x = skew, y = vr, fill = skew)) + 
    geom_boxplot() + 
    facet_grid(rows = vars(diversity2), cols = vars(taxa), scales = "free") + 
    labs(subtitle = "Oribatida rarefied to 20, Mesostigmata to 10 individuals", 
         y = NULL, x = NULL)

# soil
div %>% 
    select(coast:mst_iH) %>% 
    pivot_longer(names_to = "v1", values_to = "vr", -c(1:7)) %>% 
    filter(!is.na(vr)) %>% 
    separate(v1, into = c("taxa", "diversity"), sep = "_", extra = "merge") %>% 
    mutate(diversity2 = case_when(
        diversity == "obs_m"  ~"1. Abundance", 
        diversity == "obs_qD" ~ "2. Observed N of sp.",
        diversity == "iH"     ~ "4. Shannon index", 
        TRUE ~ "3. Rarefied N of sp."
    ), taxa = case_when(taxa == "orb" ~ "Oribatida", TRUE ~ "Mesostigmata")) %>% 
    ggplot(aes(x = soil, y = vr, fill = soil)) + 
    geom_boxplot() + 
    facet_grid(rows = vars(diversity2), cols = vars(taxa), scales = "free") + 
    labs(subtitle = "Oribatida rarefied to 20, Mesostigmata to 10 individuals", 
         y = NULL, x = NULL)

# vegetation
div %>% 
    select(coast:mst_iH) %>% 
    pivot_longer(names_to = "v1", values_to = "vr", -c(1:7)) %>% 
    filter(!is.na(vr)) %>% 
    separate(v1, into = c("taxa", "diversity"), sep = "_", extra = "merge") %>% 
    mutate(
        diversity2 = case_when(
        diversity == "obs_m"  ~"1. Abundance", 
        diversity == "obs_qD" ~ "2. Observed N of sp.",
        diversity == "iH"     ~ "4. Shannon index", 
        TRUE ~ "3. Rarefied N of sp."), 
        taxa = case_when(taxa == "orb" ~ "Oribatida", TRUE ~ "Mesostigmata"))  %>% 
    ggplot(aes(x = vegetation, y = vr, fill = taxa)) + 
    geom_boxplot() + 
    facet_grid(rows = vars(diversity2), scales = "free") + 
    labs(subtitle = "Oribatida rarefied to 20, Mesostigmata to 10 individuals", 
         y = NULL, x = "Vegetation cover, %")

# d.plants
div %>% 
    select(coast:mst_iH) %>% 
    pivot_longer(names_to = "v1", values_to = "vr", -c(1:7)) %>% 
    filter(!is.na(vr)) %>% 
    separate(v1, into = c("taxa", "diversity"), sep = "_", extra = "merge") %>% 
    mutate(diversity2 = case_when(
        diversity == "obs_m"  ~"1. Abundance", 
        diversity == "obs_qD" ~ "2. Observed N of sp.",
        diversity == "iH"     ~ "4. Shannon index", 
        TRUE ~ "3. Rarefied N of sp."
    ), taxa = case_when(taxa == "orb" ~ "Oribatida", TRUE ~ "Mesostigmata")) %>% 
    ggplot(aes(x = plants.d, y = vr, fill = taxa)) + 
    geom_boxplot() + 
    facet_grid(rows = vars(diversity2),scales = "free") + 
    labs(subtitle = "Oribatida rarefied to 20, Mesostigmata to 10 individuals", 
         y = NULL, x = "Dominant plant species") + 
    theme(axis.text.x = element_text(angle = 30, hjust = 0.9))

# Multidimensional --------------------------------------------------------
M <- PCOA %>% 
    pluck("pc") %>% 
    map_df(rbind, .id = "D") %>% 
    separate(D, into = c("taxa", "type")) %>% 
    left_join(div, by = "id") %>% 
    select(taxa:plants.d) %>% 
    mutate(
        axis1 = case_when(taxa == "ms" & type == "bin" ~ axis1*-1, TRUE ~ axis1),
        taxa = case_when(taxa == "or" ~ "Oribatida", TRUE ~ "Mesostigmata"), 
        type = case_when(type == "bin" ~ "Binary data", 
                            TRUE ~ "Numeric data")) %>% 
    filter(!is.na(coast))
# general
plotly::ggplotly(
ggplot(M, aes(x = axis1, y = axis2, color = id)) + 
    geom_point() + 
    facet_grid(rows = vars(type), cols = vars(taxa)) + 
    labs(x = NULL, y = NULL, title = "General topology") + 
    theme(legend.position = "none")
)

#coast
ggplot(M, aes(x = axis1, y = axis2, color = coast)) + 
    geom_point() + 
    stat_ellipse() +
    facet_grid(rows = vars(type), cols = vars(taxa)) + 
    labs(x = NULL, y = NULL) 

# skew
ggplot(M, aes(x = axis1, y = axis2, color = skew)) + 
    geom_point() + 
    stat_ellipse() +
    facet_grid(rows = vars(type), cols = vars(taxa)) + 
    labs(x = NULL, y = NULL) 

# soil
ggplot(M, aes(x = axis1, y = axis2, color = soil)) + 
    geom_point() + 
    stat_ellipse() +
    facet_grid(rows = vars(type), cols = vars(taxa)) + 
    labs(x = NULL, y = NULL) 

# plants.d
M %>% 
    # filter(taxa != "Mesostigmata" & type == "Numeric data") %>% 
ggplot(aes(x = axis1, y = axis2, color = plants.d)) + 
    geom_point() + 
    stat_ellipse() +
    facet_grid(rows = vars(type), cols = vars(taxa)) + 
    labs(x = NULL, y = NULL, color = "Dominant plant species")

# dissimilarity -----------------------------------------------------------
dis <- list()
dis$or.bin <- or.w %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "jaccard", binary = TRUE)
dis$or.num <- or.w %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "bray", binary = FALSE)
dis$ms.bin <- ms.w %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "jaccard", binary = TRUE)
dis$ms.num <- ms.w %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "bray", binary = FALSE)

# permanova ---------------------------------------------------------------

vegan::adonis2(dis$or.bin ~ coast + soil + vegetation + plants.d, data = div)


# For future  -------------------------------------------------------------
# hill profiles 
# cluster analysis 
