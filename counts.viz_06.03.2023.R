# loading -----------------------------------------------------------------
library(tidyverse)
theme_set(theme_bw() + theme(legend.position = "bottom"))
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
taxa <- readxl::read_excel("Caspian data_01.03.2023_SA.xlsx", sheet = "taxa") %>% 
    select(sp, area2, order)

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
div
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

# Viz me 1 ---------------------------------------------------------------
abundance <- div %>% select(id, coast, plants.d, Collembola:Oribatida_total) %>% 
    pivot_longer(names_to = "taxa", values_to = "abu", -c(id:plants.d)) %>% 
    # pivot_longer(names_to = "VAR", values_to = "VAL", -c("id", "taxa", "abu")) %>% 
    mutate(seria = substr(id, 1, 4), .before = 2, 
           taxa = str_replace_all(taxa, "_total", ""),
           taxa = factor(taxa, levels = c("Collembola", "Astigmata", 
                "Mesostigmata", "Oribatida", "Prostigmata")))

p1a <- abundance %>% 
    filter(taxa %in% c("Astigmata", "Oribatida", "Mesostigmata")) %>% 
    group_by(seria, coast, taxa) %>%  
    summarise(abu = sum(abu), 
              plants.e = paste0(unique(plants.d), collapse = " / "), 
              .groups = "drop") %>% 
    ggplot(aes(x = seria, y = abu, fill = taxa)) + 
    geom_col(width = 0.68) + 
    geom_text(mapping = aes(y = 600, label = plants.e), angle = 90, color = "black") +
    geom_text(mapping = aes(y = 900, label = coast), angle = 90, color = "black") +
    # scale_y_reverse() +
    scale_fill_manual(values = c("#C77CFF", "#00B8E7", "#F8766D")) + 
    guides(fill="none") +
    labs(y = NULL, #x = "Серия", fill = NULL) + 
         x = NULL)+ 
    theme(axis.text.x = element_text(angle = 90, hjust = -0.5))
p1b <- abundance %>% 
    group_by(seria, coast, taxa) %>% 
    summarise(abu = sum(abu), 
              plants.e = paste0(unique(plants.d), collapse = " / "), 
              .groups = "drop") %>% 
    ggplot(aes(x = seria, y = abu, fill = taxa)) + 
    geom_col(width = 0.68) + 
    geom_text(mapping = aes(y = 3500, label = plants.e), angle = 90, color = "black") +
    geom_text(mapping = aes(y = 4900, label = coast), angle = 90, color = "black") +
    scale_fill_manual(values = c("#ABA300", "#C77CFF", "#00B8E7", "#F8766D", "#00c19A"), drop = TRUE) + 
    scale_y_reverse() +
    scale_x_discrete(position = "top") +
    labs(y = NULL, x = NULL, fill = NULL) + 
    theme(axis.text.x = element_blank())

p1 <- gridExtra::grid.arrange(p1a, p1b, ncol = 1, #left = "Суммарное обилие в серии, экз.")
                              left = "Total abundance, individuals")
ggsave("plot1.pdf", plot = p1, width = 297/25, height = 210/25)    

# Viz me 2 ----------------------------------------------------------------

p2 <- div %>% 
    select(id, coast, plants.d, Oribatida = orb_obs_qD, Mesostigmata = mst_obs_qD) %>% 
    mutate_if(is.numeric, function(a){a[is.na(a)] <- 0; a}) %>% 
    mutate(seria = substr(id, 1, 4), .before = 2) %>% 
    pivot_longer(names_to = "taxa", values_to = "nsp", -c(1:4)) %>% 
    group_by(seria, taxa) %>% 
    summarise(
        coast = paste0(unique(coast), collapse = " / "), 
        plants.e = paste0(unique(plants.d), collapse = " / "), 
        nsp_mean = mean(nsp),
        nsp_sd   = sd(nsp), 
        .groups = "drop") %>% 
    mutate(ymax = nsp_mean + nsp_sd, 
           ymin = nsp_mean - nsp_sd, 
           ymin = case_when(ymin <= 0 ~ 0, TRUE ~ ymin)) %>% 
    ggplot(aes(x = seria, y = nsp_mean, 
               fill = taxa, #color = taxa,
               ymin = ymin, ymax = ymax)) + 
        geom_col(position = "dodge", width = 0.68) + 
        geom_errorbar(position = "dodge", color = "black", alpha = 0.5) + 
    geom_text(mapping = aes(y = 10, label = plants.e), angle = 90, color = "black") +
    geom_text(mapping = aes(y = 4, label = coast), angle = 90, color = "black") +
        scale_fill_manual(values = c("#00B8E7", "#F8766D")) + 
        theme(axis.text.x = element_text(angle = 90)) +
        labs(x = NULL, fill = NULL, #y = "Среднее количество видов в серии ± SD")
             y = "Average number of species in seria ± SD")
ggsave("plot2.pdf", p2, width = 297/25, height = 210/25)

# Viz me 3 ----------------------------------------------------------------
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
        type = case_when(type == "bin" ~ "Binary data (Jaccard)", 
                         TRUE ~ "Numeric data (Bray-Curtis)")) %>% 
    filter(!is.na(coast))
eig <- PCOA2 %>%
    pluck("eig") %>%
    map(~data.frame(axis1 = .x[1], axis2 = .x[2])) %>%
    map_df(rbind, .id = "a") %>%
    mutate_if(is.numeric, function(a){paste0(a, " %")}) %>%
    separate(a, into = c("taxa", "type")) %>% 
    mutate(taxa = case_when(taxa == "or" ~ "Oribatida", TRUE ~ "Mesostigmata"), 
           type = case_when(type == "bin" ~ "Binary data (Jaccard)", 
                            TRUE ~ "Numeric data (Bray-Curtis)"))

ggplot(M2, aes(x = axis1, y = axis2, color = coast)) + 
    geom_point() + 
    stat_ellipse() +
    geom_text(aes(label = axis1, x = 0, y = -0.63), color = "black", data = eig, alpha = 0.68) +
    geom_text(aes(label = axis2, x = -0.9, y = 0), color = "black", data = eig, alpha = 0.68, angle = 90) +
    scale_color_manual(values = c("#00A9FF", "#0CB720", "#CD9600", "#F8766D"))+
    facet_grid(cols = vars(type), rows = vars(taxa)) + 
    labs(x = NULL, y = NULL, color  = NULL, #subtitle = "A. Тип берега")
         subtitle = "A. Coast type") 
p3b <- ggplot(M2, aes(x = axis1, y = axis2, color = plants.d)) + 
    geom_point() + 
    stat_ellipse() +
    geom_text(aes(label = axis1, x = 0, y = -0.77), color = "black", data = eig, alpha = 0.68) +
    geom_text(aes(label = axis2, x = -0.99, y = 0), color = "black", data = eig, alpha = 0.68, angle = 90) +
    facet_grid(cols = vars(type), rows = vars(taxa)) + 
    labs(x = NULL, y = NULL, color = NULL, #subtitle = "Б. Доминантные виды растений") + 
         subtitle = "B. Dominant plant species") + 
    theme(legend.text = element_text(face = "italic"))
p3 <- gridExtra::grid.arrange(p3a, p3b, ncol = 1) #, left = "Суммарное обилие в серии, экз.")
ggsave("plot3.pdf", p3, width = 210/25, height = 297/25)

# Viz me 4 ----------------------------------------------------------------
p4 <- rbind(ms.w, or.w) %>% 
    pivot_longer(names_to = "id", values_to = "abu", -sp) %>% 
    left_join(taxa, by = "sp") %>% 
    rbind(., .) %>% 
    mutate(type = rep(c("diff", "Oribatida & Mesostigmata"), each = 14160)) %>% 
    filter(area2 != "unknown" | type != "Oribatida & Mesostigmata") %>% 
    transmute(id, seria = substr(id, 1, 4), sp, 
              area2 = factor(area2, ordered = TRUE, levels = 
                    c("(Semi)Cosmopolitan", "Holarctic", "Palaearctic", 
                      "European-Caucasian", "Mediterranean-Caucasian",
                      "Caspian and/or Caucasian")), 
              order = case_when(type == "diff" ~ order, TRUE ~ type),
              abu, fau = case_when(abu > 0 ~ 1, TRUE ~ 0)) %>% 
    group_by(seria, order, area2) %>% 
    summarise_if(is.numeric, sum) %>% 
    mutate(Community = abu/sum(abu)*100, Fauna = fau/sum(fau)*100, .keep = "unused") %>% 
    ungroup() %>% 
    pivot_longer(names_to = "VAR", values_to = "VAL", -1:-3) %>% 
ggplot(aes(x = seria, y = VAL, fill = area2)) + 
    geom_col(width = 0.68) + 
    facet_grid(cols = vars(VAR), rows = vars(order)) + 
    scale_fill_viridis_d() + 
    scale_y_reverse() +
    labs(y = "Ratio, %", x = NULL, fill = NULL) + 
    theme(axis.text.x = element_text(angle = 90))
# p4
ggsave("plot4.pdf", p4, height = 210/25, width = 297/25)
