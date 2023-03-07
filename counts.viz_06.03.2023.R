# loading -----------------------------------------------------------------
library(tidyverse)
library(parallel)
cl <- makeCluster(detectCores()-1)

long <- readxl::read_excel("Caspian data_01.03.2023_SA.xlsx", sheet = "main") %>% 
    select(O:SmSw__5) %>% 
    mutate_all(as.character) %>% 
    pivot_longer(names_to = "id", values_to = "v", -c("O", "sp")) %>% 
    separate(v, into = c("ad", "jv"), sep = "\\+", convert = TRUE) %>% # 
    pivot_longer(names_to = "age", values_to = "abu", -c(1:3)) %>% 
    filter(!is.na(abu)) %>% 
    mutate(id = substr(id, 3, 7)) %>%  #id = str_replace_all(id, "_", "")) %>% 
    arrange(id)
or.w <- long %>% 
    filter(O == "Oribatida", sp != "Oribatida_juvenile_indet") %>% 
    select(-O, -age) %>% 
    pivot_wider(names_from = id, values_from = abu, values_fn = sum) #id or taxa???
ms.w <- long %>% 
    filter(O == "Mesostigmata") %>% 
    select(-O, -age) %>%
    pivot_wider(names_from = id, values_from = abu, values_fn = sum)

labs <- readxl::read_excel("Caspian data_01.03.2023_SA.xlsx", sheet = "samples") %>% 
    filter(distr == "Samoor") %>% 
    transmute(id = substr(id, 3, 7), coast, plants.d) %>% 
    separate(plants.d, sep = " ", into = "plants.d", extra = "drop")

# rarefication ------------------------------------------------------------
# Oribatida - d20
div <- or.w %>% 
    select(-sp) %>% 
    lapply(function(a)a <- a[a>0]) %>% 
    keep(~length(.x)>0) %>% # >0 !!!
    parLapply(cl = cl, ., function(a){
        a |> 
            iNEXT::iNEXT(size = 20, nboot = 0) |>
            purrr::pluck("iNextEst", "size_based") |> 
            dplyr::select(m, Method, qD) |> 
            dplyr::filter(m == 20 | Method == "Observed") |> 
            tidyr::pivot_longer(names_to = "mm", values_to = "vv", -Method) |>
            dplyr::filter(mm == "qD" | Method == "Observed") |> 
            dplyr::mutate(Method = dplyr::case_when(Method == "Observed" ~ "orb_obs", TRUE ~ "orb_d20")) |> 
            tidyr::pivot_wider(names_from = c("Method", "mm"), values_from = vv) |> 
            dplyr::mutate(orb_iH = vegan::diversity(a))
        }) %>% 
    map_df(rbind, .id = "id") %>% 
    left_join(labs, ., by = "id")
# Mesostigmata d10
div <- ms.w %>% 
    select(-sp) %>% 
    lapply(function(a)a <- a[a>0]) %>% 
    keep(~length(.x)>0) %>% # >0 !!!
    parLapply(cl = cl, ., function(a){
        a |> 
            iNEXT::iNEXT(size = 10, nboot = 0) |>
            purrr::pluck("iNextEst", "size_based") |> 
            dplyr::select(m, Method, qD) |> 
            dplyr::filter(m == 10 | Method == "Observed") |> 
            tidyr::pivot_longer(names_to = "mm", values_to = "vv", -Method) |>
            dplyr::filter(mm == "qD" | Method == "Observed") |> 
            dplyr::mutate(Method = dplyr::case_when(Method == "Observed" ~ "mst_obs", TRUE ~ "mst_d10")) |> 
            tidyr::pivot_wider(names_from = c("Method", "mm"), values_from = vv) |> 
            dplyr::mutate(mst_iH = vegan::diversity(a))
    }) %>% 
    map_df(rbind, .id = "id") %>% 
    left_join(div, ., by = "id")
# general
div <- long %>% 
    filter(O == "total", 
           !(sp %in% c("Oribatida_adult", 
                       "Oribatida_juvenile_det", 
                       "Oribatida_juvenile_indet"))) %>% 
    pivot_wider(names_from = sp, values_from = abu) %>% 
    select(-O, -age) %>% 
    left_join(div, ., by = "id") 

# multivariate ------------------------------------------------------------
dis <- list()
# dissimilarity
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

# pcoa
PCOA <- dis %>% 
    parLapply(cl = cl, ., function(a){
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
rm(cl)

# Viz me ---------------------------------------------------------------



# save.image("01.03.2023.RData")

