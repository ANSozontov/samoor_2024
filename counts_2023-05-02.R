# loading -----------------------------------------------------------------
library(tidyverse)
theme_set(theme_bw() + theme(legend.position = "bottom"))

total <- readxl::read_excel("Caspian data_01.03.2023_SA.xlsx", sheet = "main") %>% 
    select(O:SmSdTu5) %>% 
    mutate_all(as.character) %>% 
    pivot_longer(names_to = "id", values_to = "v", -c("O", "sp")) %>% 
    separate(v, into = c("ad", "jv"), sep = "\\+", convert = TRUE) %>% # 
    mutate_if(is.numeric, function(a){a[is.na(a)] <- 0; a}) %>% 
    transmute(id = substr(id, 3, 7), seria = substr(id, 1, 4),
              O, sp, ad, jv, abu = ad+jv) %>% 
    # filter(!is.na(abu), abu > 0) %>% 
    arrange(id)
long <- total %>% 
    filter(O != "total", O != "Astigmata") %>% 
    select(-ad, -jv)
orwi <- long %>% 
    filter(O == "Oribatida", sp != "Oribatida_juvenile_indet") %>% 
    select(-O, -seria) %>% 
    pivot_wider(names_from = id, values_from = abu, values_fn = sum, values_fill = 0)
orws <- long %>% 
    filter(O == "Oribatida", sp != "Oribatida_juvenile_indet") %>% 
    select(-O, -id) %>% 
    pivot_wider(names_from = seria, values_from = abu, values_fn = sum, values_fill = 0)
mswi <- long %>% 
    filter(O == "Mesostigmata") %>% 
    select(-O, -seria) %>%
    pivot_wider(names_from = id, values_from = abu, values_fn = sum, values_fill = 0)
msws <- long %>% 
    filter(O == "Mesostigmata") %>% 
    select(-O, -id) %>%
    pivot_wider(names_from = seria, values_from = abu, values_fn = sum, values_fill = 0)

labs <- readxl::read_excel("Caspian data_01.03.2023_SA.xlsx", sheet = "samples") %>% 
    filter(distr == "Samoor", code != "SmSw") %>% 
    transmute(id = substr(id, 3, 7), seria = substr(id, 1, 4), plants.d, 
              coast = case_when(coast == "sandy dunes" ~ "dunes", TRUE ~ coast)) %>% 
    separate(plants.d, sep = " ", into = "plants.d", extra = "drop") %>% 
    arrange(id)
# taxa <- readxl::read_excel("Caspian data_01.03.2023_SA.xlsx", sheet = "taxa") %>% 
#     select(sp, area2, order)
# 
# xaxis <- c(
#     "SdJj", #sandy beach
#     "SdEq", #sandy beach
#     "SdTu", #sandy beach
#     "SdEc", #sandy beach
#     "SdTa", #sandy beach
#     "SdJm", #sandy beach
#     "SdFn", #sandy beach
#     "PbAe", #pebbly     
#     "PbDe", #pebbly     
#     "PbPo", #pebbly     
#     "PbTu", #pebbly     
#     "PbTl", #pebbly     
#     "SdCJ", #dunes      
#     "SdCS", #dunes   
#     "RsFd") #reeds      

# rarefication ------------------------------------------------------------
library(parallel)
cl <- makeCluster(detectCores()-1)
# Oribatida - d20
rar <- orwi %>% 
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
rar <- mswi %>% 
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
    left_join(rar, ., by = "id") %>% 
    mutate_at(c("orb_obs_m", "orb_obs_qD", "mst_obs_m", "mst_obs_qD"), 
              function(a){a[is.na(a)] <- 0; a})

# comparison --------------------------------------------------------------
rar %>% 
    select(coast:mst_iH) %>% 
    pivot_longer(names_to = "v1", values_to = "vr", -1) %>% 
    filter(!is.na(vr)) %>% 
    separate(v1, into = c("taxa", "diversity"), sep = "_", extra = "merge") %>% 
    mutate(
        diversity2 = case_when(
            diversity == "obs_m"  ~"1. Abundance", 
            diversity == "obs_qD" ~ "2. Observed N of sp.",
            diversity == "iH"     ~ "4. Shannon index", 
            TRUE ~ "3. Rarefied N of sp."), 
        taxa = case_when(
            taxa == "orb" ~ "Oribatida", 
            TRUE ~ "Mesostigmata")) %>% 
    ggplot(aes(x = coast, y = vr, fill = taxa)) + 
    geom_boxplot() + 
    facet_wrap(~diversity2, scales = "free") +
    labs(subtitle = "Oribatida rarefied to 20, Mesostigmata to 10 individuals", 
         y = NULL, x = NULL, fill = "")

#oribatida
fit_or.plant <- lm(orb_obs_m ~ plants.d, data = rar) 
fit_or.coast <- lm(orb_obs_m ~ coast, data = rar) 

ggplot() + 
    geom_density(aes(x = residuals(fit_or.coast)), color = "blue") +
    geom_density(aes(x = residuals(fit_or.plant)), color = "green") 
shapiro.test(residuals(fit_or.coast))
shapiro.test(residuals(fit_or.plant))
AIC(fit_or.coast)
AIC(fit_or.plant)
summary(fit_or.coast)
summary(fit_or.plant)

# mesostigmata
fit_ms.plant <- glm(mst_obs_m ~ plants.d, data = rar, family = binomial(link = "logit")) 
fit_ms.coast <- glm(mst_obs_m ~ coast, data = rar, family = binomial(link = "logit")) 
fit <- lme4::lmer(mst_obs_m ~ coast + (1|seria), data = rar) 

ggplot() + 
    geom_density(aes(x = residuals(fit_ms.coast)), color = "blue") +
    geom_density(aes(x = residuals(fit_ms.plant)), color = "green") 
shapiro.test(residuals(fit_ms.coast))
shapiro.test(residuals(fit_ms.plant))
AIC(fit_ms.coast)
AIC(fit_ms.plant)
summary(fit_ms.coast)
summary(fit_ms.plant)




# multivariate ------------------------------------------------------------
dis <- list()
# dissimilarity

dis <- list(orws, orws, msws,  msws #, rbind(orws, msws), rbind(orws, msws)
            ) %>% 
    map(~.x %>% 
            column_to_rownames("sp") %>% 
            select_if(function(a){sum(a)>0}) %>% 
            t %>% 
            as.data.frame)
names(dis) <- c("or.bin", "or.num", "ms.bin", "ms.num") #, "all.bin", "all.num")
dis$or.bin <- vegan::vegdist(dis$or.bin, method = "jaccard", binary = TRUE)
dis$or.num <- vegan::vegdist(dis$or.num, method = "bray", binary = FALSE)
dis$ms.bin <- vegan::vegdist(dis$ms.bin, method = "jaccard", binary = TRUE)
dis$ms.num <- vegan::vegdist(dis$ms.num, method = "bray", binary = FALSE)
# dis$all.bin<- vegan::vegdist(dis$all.bin, method = "jaccard", binary = TRUE)
# dis$all.num<- vegan::vegdist(dis$all.num, method = "bray", binary = FALSE)

# pcoa
PCOA <- dis %>% 
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

eig <- PCOA %>%
    pluck("eig") %>%
    map(~data.frame(axis1 = .x[1], axis2 = .x[2])) %>%
    map_df(rbind, .id = "a") %>%
    mutate_if(is.numeric, function(a){paste0(a, " %")}) %>%
    separate(a, into = c("taxa", "type")) %>% 
    mutate(taxa = case_when(taxa == "or" ~ "Oribatida", 
                            taxa == "ms" ~ "Mesostigmata", 
                            TRUE ~ "Both orders"), 
           type = case_when(type == "bin" ~ "Binary data (Jaccard)", 
                            TRUE ~ "Numeric data (Bray-Curtis)"))
M2 <- PCOA %>% 
    pluck("pc") %>% 
    map_df(rbind, .id = "D") %>% 
    separate(D, into = c("taxa", "type")) %>% 
    left_join(distinct(select(labs, id = seria, coast)), by = "id") %>% 
    mutate(
        axis1 = case_when(taxa == "all" & type == "bin" ~ axis1*-1, TRUE ~ axis1),
        # axis1 = case_when(taxa == "or" & type == "bin" ~ axis1*-1, TRUE ~ axis1),
        # axis2 = case_when(taxa == "or" & type == "num" ~ axis2*-1, TRUE ~ axis2),
        taxa = case_when(taxa == "or" ~ "Oribatida", 
                         taxa == "ms" ~ "Mesostigmata", 
                         TRUE ~ "Both orders"), 
        type = case_when(type == "bin" ~ "Binary data (Jaccard)", 
                         TRUE ~ "Numeric data (Bray-Curtis)")) 
M3 <- M2 %>% 
    group_by(taxa, type, coast) %>% 
    slice(chull(axis1, axis2))

ggplot() + 
    geom_polygon(aes(axis1, axis2, color = coast, fill = coast), data = M3, alpha = 0.2) +
    geom_point(aes(axis1, axis2, color = coast, fill = coast), 
               data = M2, 
               color = "black", shape = 21) + 
    geom_text(aes(x = 0, y = -0.5, label = axis1), 
              data = eig, size = 3)+
    geom_text(aes(x = -0.4, y = 0, label = axis2), 
              data = eig, size = 3, angle = 90)+
    facet_grid(cols = vars(type), rows = vars(taxa)) +
    theme(panel.grid = element_blank())

# indval ------------------------------------------------------------------
library(indicspecies)

set.seed(2)
iv <- orws %>% 
    rbind(msws) %>%
    mutate(total = apply(.[,-1], 1, sum), .before = 2) %>% 
    filter(total >= 5) %>% 
    select(-total) %>% 
    column_to_rownames("sp") %>% 
    t %>%
    indicspecies::multipatt(
        pull(arrange(distinct(select(labs, seria, coast)), seria) , coast), 
        control = how(nperm=999),
        max.order = 4,
        func = "indval",
        duleg = FALSE
    )
summary(iv)
iv <- iv$str %>% 
    as.data.frame() %>% 
    rownames_to_column("sp") %>% 
    as_tibble() %>% 
    pivot_longer(names_to = "biotop", values_to = "iv", -sp) %>% 
    group_by(sp) %>% 
    filter(iv == max(iv)) %>% 
    ungroup() %>% 
    left_join(select(rownames_to_column(iv$sign, "sp"), sp, p.value), 
              by = "sp") %>% 
    mutate(`iv, %` = round(iv*100), 
           p.value = round(p.value, 4),
           sign = case_when(p.value <= 0.05 ~ "*", TRUE ~ NA), 
           .keep = "unused", 
           .after = 2)




