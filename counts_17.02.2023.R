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


# dissimilarity -----------------------------------------------------------
dis2 <- list()
dis2$or.bin <- or.w %>% 
    select(!starts_with("Sw")) %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "jaccard", binary = TRUE)
dis2$or.num <- or.w %>% 
    select(!starts_with("Sw")) %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "bray", binary = FALSE)
dis2$ms.bin <- ms.w %>% 
    select(!starts_with("Sw")) %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "jaccard", binary = TRUE)
dis2$ms.num <- ms.w %>% 
    select(!starts_with("Sw")) %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "bray", binary = FALSE)

# Multidimensional --------------------------------------------------------
PCOA2 <- dis2 %>% 
    lapply(function(a){
        p <- ape::pcoa(a)
        e <- p$values$Eigenvalues
        if(min(e) < 0){
            e <- e + abs(min(e))
            e <- round(e/sum(e)*100, 1)
        } else { 
            e <- round(e/sum(e)*100, 1)
        }
        p <- tibble::tibble(id = rownames(p$vectors), 
                            axis1 = p$vectors[,1], 
                            axis2 = p$vectors[,2]) 
        list(eig = e, pc = p)
    }) %>% 
    purrr::transpose()

M2 <- PCOA2 %>% 
    pluck("pc") %>% 
    map_df(rbind, .id = "D") %>% 
    separate(D, into = c("taxa", "type")) %>% 
    left_join(div, by = "id") %>% 
    select(taxa:plants.d) %>% 
    mutate(
        axis1 = case_when(taxa == "ms" & type == "bin" ~ axis1*-1, TRUE ~ axis1),
        axis1 = case_when(taxa == "or" & type == "bin" ~ axis1*-1, TRUE ~ axis1),
        taxa = case_when(taxa == "or" ~ "Oribatida", TRUE ~ "Mesostigmata"), 
        type = case_when(type == "bin" ~ "Binary data", 
                         TRUE ~ "Numeric data")) %>% 
    filter(!is.na(coast))

M <- PCOA %>% 
    pluck("pc") %>% 
    map_df(rbind, .id = "D") %>% 
    separate(D, into = c("taxa", "type")) %>% 
    left_join(div, by = "id") %>% 
    select(taxa:plants.d) %>% 
    mutate(
        axis1 = case_when(taxa == "ms" & type == "bin" ~ axis1*-1, TRUE ~ axis1),
        # axis1 = case_when(taxa == "or" & type == "bin" ~ axis1*-1, TRUE ~ axis1),
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
ggplot(M2, aes(x = axis1, y = axis2, color = coast)) + 
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
ggplot(M2, aes(x = axis1, y = axis2, color = skew)) + 
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
M2 %>% 
    # filter(taxa != "Mesostigmata" & type == "Numeric data") %>% 
    ggplot(aes(x = axis1, y = axis2, color = plants.d)) + 
    geom_point() + 
    stat_ellipse() +
    facet_grid(rows = vars(type), cols = vars(taxa)) + 
    labs(x = NULL, y = NULL, color = "Dominant plant species")

# permanova ---------------------------------------------------------------
vegan::adonis2(dis2$or.bin ~ coast + soil + vegetation + plants.d, 
               data = filter(labs, id %in% attr(dis2$or.bin, "Labels")))
vegan::adonis2(dis2$or.num ~ coast + soil + vegetation + plants.d, 
               data = filter(labs, id %in% attr(dis2$or.num, "Labels")))
vegan::adonis2(dis2$ms.bin ~ coast + soil + vegetation + plants.d, 
               data = filter(labs, id %in% attr(dis2$ms.bin, "Labels")))
vegan::adonis2(dis2$ms.num ~ coast + soil + vegetation + plants.d, 
               data = filter(labs, id %in% attr(dis2$ms.num, "Labels")))


# For future  -------------------------------------------------------------
# hill profiles 
# cluster analysis 

# map ---------------------------------------------------------------------
# library(leaflet)
# 
# readxl::read_excel("Caspian data_03.02.2023_SA.xlsx", sheet = "samples") %>% 
#     filter(distr == "Samoor") %>% 
#     transmute(id = substr(id, 3, 7), coast, skew, 
#               soil, substrate, veg, plants.d, N, E)




