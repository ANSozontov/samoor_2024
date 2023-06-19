# loading -----------------------------------------------------------------
library(tidyverse)
theme_set(theme_bw() + theme(legend.position = "bottom"))
library(parallel)
cl <- makeCluster(detectCores()-1)

long <- readxl::read_excel("Caspian data_01.03.2023_SA.xlsx", sheet = "main") %>% 
    select(O:SmSdTu5) %>% 
    mutate_all(as.character) %>% 
    pivot_longer(names_to = "id", values_to = "v", -c("O", "sp")) %>% 
    separate(v, into = c("ad", "jv"), sep = "\\+", convert = TRUE) %>% # 
    pivot_longer(names_to = "age", values_to = "abu", -c(1:3)) %>% 
    mutate(id = substr(id, 3, 7)) %>%  # seria = substr(id, 1, 4)) %>%  #id = str_replace_all(id, "_", "")) %>% 
    filter(!is.na(abu)) %>% # seria != "Sw__") %>% 
    arrange(id)
or.w <- long %>% 
    filter(O == "Oribatida", sp != "Oribatida_juvenile_indet") %>% 
    select(-O, -age) %>% 
    pivot_wider(names_from = id, values_from = abu, values_fn = sum, values_fill = 0) #id or taxa???
ms.w <- long %>% 
    filter(O == "Mesostigmata") %>% 
    select(-O, -age) %>%
    pivot_wider(names_from = id, values_from = abu, values_fn = sum, values_fill = 0)

labs <- readxl::read_excel("Caspian data_01.03.2023_SA.xlsx", sheet = "samples") %>% 
    filter(distr == "Samoor", code != "SmSw") %>% 
    transmute(id = substr(id, 3, 7), seria = substr(id, 1, 4), plants.d, 
              coast = case_when(coast == "sandy dunes" ~ "dunes", TRUE ~ coast)) %>% 
    separate(plants.d, sep = " ", into = "plants.d", extra = "drop")
taxa <- readxl::read_excel("Caspian data_01.03.2023_SA.xlsx", sheet = "taxa") %>% 
    select(sp, area2, order)

# Рис. 4. Обилие микроартропод по типам берега и парцеллам  ---------------
abundance <- long %>% 
    filter(O == "total", 
           !(sp %in% c("Oribatida_adult", 
                       "Oribatida_juvenile_det", 
                       "Oribatida_juvenile_indet"))) %>% 
    pivot_wider(names_from = sp, values_from = abu) %>% 
    select(-O, -age) %>% 
    left_join(labs, ., by = "id") %>% 
    pivot_longer(names_to = "taxa", values_to = "abu", -c(id:coast)) %>% 
    mutate(
        seria = substr(id, 1, 4), 
        taxa = str_replace_all(taxa, "_total", ""),
        taxa = factor(taxa, levels = c("Collembola", "Astigmata", 
            "Mesostigmata", "Oribatida", "Prostigmata")), 
        .before = 2)
xaxis <- c("SdJj", "SdEq", "SdTu", "SdEc", "SdTa", "SdJm", "SdFn", 
    "PbAe", "PbDe", "PbPo", "PbTu", "PbTl", "SdCJ", "SdCS", "RsFd") 

p4a <- abundance %>% 
    group_by(seria, coast, taxa) %>% 
    summarise(abu = sum(abu), 
              plants.e = paste0(unique(plants.d), collapse = " / "), 
              .groups = "drop") %>% 
    ggplot(aes(x = seria, y = abu, fill = taxa)) + 
    geom_col(width = 0.68) + 
    geom_text(mapping = aes(y = 3500, label = plants.e), angle = 90, color = "black", fontface = "italic") +
    geom_text(mapping = aes(y = 1900, label = coast), angle = 90, color = "black") +
    scale_fill_manual(values = c("#ABA300", "#C77CFF", "#00B8E7", "#F8766D", "#00c19A"), drop = TRUE) + 
    scale_x_discrete(limits = xaxis) +
    labs(y = NULL, x = NULL, fill = NULL) + 
    theme(axis.text.x = element_text(angle = 90, hjust = -0.5))
p4b <- abundance %>% 
    filter(taxa %in% c("Oribatida", "Mesostigmata")) %>% 
    group_by(seria, coast, taxa) %>%  
    summarise(abu = sum(abu), 
              plants.e = paste0(unique(plants.d), collapse = " / "), 
              .groups = "drop") %>% 
    ggplot(aes(x = seria, y = abu, fill = taxa)) + 
    geom_col(width = 0.68) + 
    geom_text(mapping = aes(y = 600, label = plants.e), angle = 90, color = "black", fontface = "italic") +
    geom_text(mapping = aes(y = 300, label = coast), angle = 90, color = "black") +
    scale_x_discrete(limits = xaxis) +
    scale_fill_manual(values = c("#00B8E7", "#F8766D")) + 
    guides(fill="none") +
    labs(y = NULL, x = NULL)+ 
    theme(axis.text.x = element_blank())
p4 <- gridExtra::grid.arrange(p4a, p4b, ncol = 1, 
                              left = "Total abundance, individuals")
ggsave("plot_4.pdf", plot = p4, width = 297/25, height = 210/25)   

# Табл. 2. Корреляция обилия групп друг с другом и RH%, C%. ---------------

# Табл. 3. Таксономический состав панцирных и мезостигматических к --------
## НИПОНЕЛ, обсудить

# Рис. 5. Куммуляты по типам берега (АС). ---------------------------------

# rar <- 
long %>% 
    filter(O == "Oribatida") %>% 
    left_join(labs, by = "id") %>% 
    group_by(coast, sp) %>% 
    summarise(abu = mean(abu), .groups = "drop") %>% 
    filter(abu > 0) %>% 
    pivot_wider(names_from = coast, values_from = abu, values_fill = 0) %>% 
    select(-sp) %>% 
    lapply(function(a)a <- a[a>0]) %>% 
    keep(~length(.x)>0) %>% # >0 !!!
    parLapply(cl = cl, ., function(a){
        a |> 
            iNEXT::iNEXT(size = seq(0, 100, by = 5), nboot = 0) |>
            purrr::pluck("iNextEst", "size_based") |> 
            dplyr::select(m, Method, qD) |> 
            dplyr::filter(m %in% seq(0, 100, by = 5) | Method == "Observed") |> 
            dplyr::mutate(taxa = "Oribatida")
    }) 
    map_df(rbind, .id = "id") %>% 
    as_tibble() %>% 
    left_join(labs, by = "id")


# Oribatida - d20
# rar <- 
or.w %>% 
    select(-sp) %>% 
    lapply(function(a)a <- a[a>0]) %>% 
    keep(~length(.x)>0) %>% # >0 !!!
    parLapply(cl = cl, ., function(a){
        a |> 
            iNEXT::iNEXT(size = seq(0, 50, by = 2), nboot = 0) |>
            purrr::pluck("iNextEst", "size_based") |> 
            dplyr::select(m, Method, qD) |> 
            dplyr::filter(m %in% seq(0, 50, by = 2) | Method == "Observed") |> 
            dplyr::mutate(taxa = "Oribatida")
            # tidyr::pivot_longer(names_to = "mm", values_to = "vv", -Method) |>
            # dplyr::filter(mm == "qD" | Method == "Observed") |> 
            # dplyr::mutate(Method = dplyr::case_when(Method == "Observed" ~ "orb_obs", TRUE ~ "orb_d20")) |> 
            # tidyr::pivot_wider(names_from = c("Method", "mm"), values_from = vv) #|> 
            #dplyr::mutate(orb_iH = vegan::diversity(a))
    }) %>% 
    map_df(rbind, .id = "id") %>% 
    as_tibble() %>% 
    left_join(labs, by = "id")
# Mesostigmata d10
rar <- ms.w %>% 
    select(-sp) %>% 
    lapply(function(a)a <- a[a>0]) %>% 
    keep(~length(.x)>0) %>% # >0 !!!
    parLapply(cl = cl, ., function(a){
        a |> 
            iNEXT::iNEXT(size = seq(0, 25, by = 1), nboot = 0) |>
            purrr::pluck("iNextEst", "size_based") |> 
            dplyr::select(m, Method, qD) |> 
            dplyr::filter(m %in% seq(0, 25, by = 1) | Method == "Observed") |> 
            dplyr::mutate(taxa = "Mesostigmata")
            # tidyr::pivot_longer(names_to = "mm", values_to = "vv", -Method) |>
            # dplyr::filter(mm == "qD" | Method == "Observed") |> 
            # dplyr::mutate(Method = dplyr::case_when(Method == "Observed" ~ "mst_obs", TRUE ~ "mst_d10")) |> 
            # tidyr::pivot_wider(names_from = c("Method", "mm"), values_from = vv) |> 
            # dplyr::mutate(mst_iH = vegan::diversity(a))
    }) %>% 
    map_df(rbind, .id = "id") %>% 
    as_tibble()
    left_join(rar, ., by = "id")
# general
div <- long %>% 
    filter(O == "total", 
           !(sp %in% c("Oribatida_adult", 
                       "Oribatida_juvenile_det", 
                       "Oribatida_juvenile_indet"))) %>% 
    pivot_wider(names_from = sp, values_from = abu) %>% 
    select(-O, -age) %>% 
    left_join(rar, ., by = "id") %>% 
    mutate(seria = substr(id, 1, 4), .before = 2) 


# Рис. 6. Кладограмма по фаунистическим спискам отдельных мероцено --------


# Рис. 7. Видовое богатство по сериям  ------------------------------------

div <- or.w %>% 
    select(-sp) %>% 
    as.list() %>% 
    lapply(function(a){data.frame(Oribatida = length(a[a>0]))}) %>% 
    map_df(rbind, .id = "id")
div <- ms.w %>% 
    select(-sp) %>% 
    as.list() %>% 
    lapply(function(a){data.frame(Mesostigmata = length(a[a>0]))}) %>% 
    map_df(rbind, .id = "id") %>% 
    left_join(div, by = "id") %>% 
    left_join(labs, ., by = "id") %>% 
    as_tibble()

p7 <- div %>% 
    # select(id, seria, coast, plants.d, Oribatida = orb_obs_qD, Mesostigmata = mst_obs_qD) %>% 
    # mutate_if(is.numeric, function(a){a[is.na(a)] <- 0; a}) %>% 
    # mutate(seria = substr(id, 1, 4), .before = 2) %>% 
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
    scale_x_discrete(limits = xaxis) +
    scale_fill_manual(values = c("#00B8E7", "#F8766D")) + 
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = NULL, fill = NULL, #y = "Среднее количество видов в серии ± SD")
         y = "Average number of species in seria ± SD")
ggsave("plot_7.pdf", p2, width = 297/25, height = 210/25)

# Рис. 8. Ареалогический состав по сериям: качественный и количест --------

p8 <- rbind(mutate(ms.w, Order = "Mesostigmata", .before = 1), 
            mutate(or.w, Order = "Oribatida", .before = 1), 
            mutate(rbind(ms.w, or.w), .before = 1,
                   Order = "Oribatida & Mesostigmata")) %>% 
    pivot_longer(names_to = "id", values_to = "abu", -Order:-sp) %>% 
    left_join(taxa, by = "sp") %>% 
    filter(area2 != "unknown" | Order != "Oribatida & Mesostigmata", 
           abu > 0) %>%
    transmute(
              seria = substr(id, 1, 4), 
              area2 = factor(area2, ordered = TRUE, levels = 
                                 c("(Semi)Cosmopolitan", "Holarctic", "Palaearctic", 
                                   "European-Caucasian", "Mediterranean-Caucasian",
                                   "Caspian and/or Caucasian")), 
              id, Order, sp, abu, 
              fau = case_when(abu > 0 ~ 1, TRUE ~ 0)) %>% 
    group_by(seria, Order, area2) %>% 
    summarise_if(is.numeric, sum) %>% 
    mutate(Community = abu/sum(abu)*100, Fauna = fau/sum(fau)*100, .keep = "unused") %>% 
    ungroup() %>% 
    pivot_longer(names_to = "VAR", values_to = "VAL", -1:-3) %>% 
    # filter(order  == "Oribatida & Mesostigmata") %>% 
    ggplot(aes(x = seria, y = VAL, fill = area2)) + 
        geom_col(width = 0.68) + 
        facet_grid(cols = vars(VAR), rows = vars(Order)) + 
        scale_fill_manual(values = c(
            "#0D0786", 
            "#3FA9F5", 
            "#79D151", 
            "#8F2773", 
            "#FB9009", 
            "#FCE724", 
            NA)) +
        scale_y_reverse() +
        scale_x_discrete(limits = xaxis) +
        labs(y = "Ratio, %", x = NULL, fill = NULL) + 
        theme(axis.text.x = element_text(angle = 90))
p8
ggsave("plot_8.pdf", p4, height = 210/40, width = 297/40)

# Табл. 5. Виды-специалисты (индикаторы) ----------------------------------


# Рис. 9. ССА распределения видов -----------------------------------------
# Факторы: тип берега (? 1-4), 
# тип растительности (рогозы, злаки, ситники, хвощ, лох, тростники), 
# расстояние до моря (м), RH%, C%  


# Рис. 10. Ординация мероценозов Mesostigmata и Oribatida  ----------------
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

M2 <- PCOA %>% 
    pluck("pc") %>% 
    map_df(rbind, .id = "D") %>% 
    separate(D, into = c("taxa", "type")) %>% 
    left_join(labs, by = "id") %>% 
    select(taxa:coast) %>% 
    mutate(
        axis1 = case_when(taxa == "ms" & type == "bin" ~ axis1*-1, TRUE ~ axis1),
        axis1 = case_when(taxa == "or" & type == "bin" ~ axis1*-1, TRUE ~ axis1),
        axis2 = case_when(taxa == "or" & type == "num" ~ axis2*-1, TRUE ~ axis2),
        taxa = case_when(taxa == "or" ~ "Oribatida", TRUE ~ "Mesostigmata"), 
        type = case_when(type == "bin" ~ "Binary data (Jaccard)", 
                         TRUE ~ "Numeric data (Bray-Curtis)")) %>% 
    filter(!is.na(coast))
eig <- PCOA %>%
    pluck("eig") %>%
    map(~data.frame(axis1 = .x[1], axis2 = .x[2])) %>%
    map_df(rbind, .id = "a") %>%
    mutate_if(is.numeric, function(a){paste0(a, " %")}) %>%
    separate(a, into = c("taxa", "type")) %>% 
    mutate(taxa = case_when(taxa == "or" ~ "Oribatida", TRUE ~ "Mesostigmata"), 
           type = case_when(type == "bin" ~ "Binary data (Jaccard)", 
                            TRUE ~ "Numeric data (Bray-Curtis)"))


p3a <- ggplot(M2, aes(x = axis1, y = axis2, color = plants.d)) + 
    geom_point() + 
    stat_ellipse() +
    geom_text(aes(label = axis1, x = 0, y = -0.77), color = "black", data = eig, alpha = 0.68) +
    geom_text(aes(label = axis2, x = -0.99, y = 0), color = "black", data = eig, alpha = 0.68, angle = 90) +
    facet_grid(cols = vars(type), rows = vars(taxa)) + 
    labs(x = NULL, y = NULL, color = NULL, #subtitle = "Б. Доминантные виды растений") + 
         subtitle = "A. Dominant plant species") + 
    theme(legend.text = element_text(face = "italic"))
p3b <- ggplot(M2, aes(x = axis1, y = axis2, color = coast)) + 
    geom_point() + 
    stat_ellipse() +
    geom_text(aes(label = axis1, x = 0, y = -0.63), color = "black", data = eig, alpha = 0.68) +
    geom_text(aes(label = axis2, x = -0.9, y = 0), color = "black", data = eig, alpha = 0.68, angle = 90) +
    scale_color_manual(values = c("#F8766D", "#00A9FF", "#0CB720", "#CD9600"))+
    facet_grid(cols = vars(type), rows = vars(taxa)) + 
    labs(x = NULL, y = NULL, color  = NULL, #subtitle = "A. Тип берега")
         subtitle = "B. Coast type") 
p3 <- gridExtra::grid.arrange(p3a, p3b, ncol = 1) #, left = "Суммарное обилие в серии, экз.")
ggsave("plot_3.pdf", p3, width = 210/25, height = 297/25)




