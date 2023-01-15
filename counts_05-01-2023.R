# 1. Load data ------------------------------------------------------------
suppressMessages(library(tidyverse))
library(parallel) # are you sure
# library(sf)
# library(vegan)
source("indval.R")

taxa <- readxl::read_xlsx(
        "Caspian data_15.09.2022_SA.xlsx", 
        sheet = "taxa") %>% 
    select(sp, order, area)
labs <- readxl::read_xlsx(
        "Caspian data_15.09.2022_SA.xlsx", 
        sheet = "samples") %>% 
    filter(distr == "Samoor") %>% 
    transmute(id, code, 
              coast, skew, soil, substrate, zone, RH,
              veg = factor(veg, ordered = TRUE))
dfl <- readxl::read_xlsx("Caspian data_15.09.2022_SA.xlsx", sheet = "main") %>% 
    mutate_all(as.character) %>% #equalize all columns
    pivot_longer(names_to = "id", values_to = "abu", -sp) %>% 
    filter(#abu != "0",
           id %in% labs$id, 
           sp != "Oribatida Juvenile instars") %>% 
    left_join(taxa, by = "sp") %>% 
    filter(order == "Oribatida") %>% 
    separate(col = abu, into = c("adu", "juv"), sep = "\\+", fill = "right") %>% 
    mutate(adu = as.numeric(adu), juv = as.numeric(juv), 
           juv = case_when(is.na(juv) ~ 0, TRUE ~ juv)) %>% 
    transmute(sp, id, abu = adu + juv, 
              abu = case_when(is.na(abu) ~ 0, TRUE ~ abu)) 
dfw <- dfl %>% 
    pivot_wider(names_from = id, values_from = abu, values_fill = 0)

# 3. rarefication ----------------------------------------------------------
cl <- makeCluster(detectCores()-1)
rar1 <- dfw %>% 
    select(-sp) %>% 
    lapply(function(a){sort(a[a>0], decreasing = TRUE)}) %>% 
    discard(~ length(.x) < 2) %>% 
    parLapply(cl, ., function(a){
        b <- iNEXT::iNEXT(a, q = 0, size = c(25), 
                 datatype = "abundance", nboot = 999) # 999
        b <- b[["iNextEst"]][["size_based"]]
        b[b$m == 25,c("qD")]
    }) %>% 
    map_dbl(c) %>% 
    tibble(id = names(.), rar25 = .)
stopCluster(cl)

# diversity & PCoA --------------------------------------------------------
div <- dfw %>% 
    column_to_rownames("sp") %>% 
    t() %>% 
    as.data.frame() %>% 
    transmute(
        abu  = apply(., 1, function(a){sum(a)}), 
        nsp  = apply(., 1, function(a){length(a[a>0])}),
        shan = apply(., 1, function(a){vegan::diversity(a, "shannon")})
    ) %>% 
    rownames_to_column("id") %>% 
    full_join(labs, ., by = "id") %>% 
    left_join(rar1, by = "id") %>% 
    lapply(function(a){if(!is.factor(a)){attributes(a) <- NULL}; a}) %>%
    as_tibble()
    
dist_n <- dfw %>% 
    select_at(div$id[div$nsp > 1]) %>% 
    t() %>% 
    vegan::vegdist(method = "bray", binary = FALSE) 
pcoa_n <- ape::pcoa(dist_n)
pcoa_n <- pcoa_n$vectors %>% 
    as.data.frame() %>% 
    select(axis1n = 1, axis2n = 2) %>% 
    rownames_to_column("id") %>% 
    mutate(eig_n = c(
        (pcoa_n$values$Eigenvalues/sum(pcoa_n$values$Eigenvalues)*100)[1:2], 
        rep(NA, nrow(.)-2)))

dist_q <- dfw %>% 
    select_at(div$id[div$nsp > 1]) %>% 
    t() %>% 
    vegan::vegdist(method = "jaccard", binary = TRUE) 
pcoa_q <- ape::pcoa(dist_q)
pcoa_q <- pcoa_q$vectors %>% 
    as.data.frame() %>% 
    select(axis1q = 1, axis2q = 2) %>% 
    rownames_to_column("id") %>% 
    mutate(eig_q = c(
        (pcoa_q$values$Eigenvalues/sum(pcoa_q$values$Eigenvalues)*100)[1:2], 
        rep(NA, nrow(.)-2)))

div <- div %>% 
    left_join(pcoa_n, by = "id") %>% 
    left_join(pcoa_q, by = "id")

# points <- readxl::read_xlsx(
#     "Caspian data_15.09.2022_SA.xlsx", 
#     sheet = "samples") %>% 
#     filter(distr == "Samoor") 

# export data -------------------------------------------------------------
rm(list = ls()[!(ls() %in% c("dfl", "dfw", "labs", "taxa", "div", "dist_n", "dist_q"))])
save.image(paste0(Sys.Date(), "_results.RData"))


