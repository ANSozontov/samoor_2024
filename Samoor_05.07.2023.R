# loading -----------------------------------------------------------------
library(indicspecies)
library(tidyverse)
library(dendextend)
library(parallel)
cl <- makeCluster(detectCores()-1)
theme_set(theme_bw() + theme(legend.position = "bottom"))

version
cat("The package's used versions\n"); sessionInfo() %>% 
    pluck("otherPkgs") %>% 
    map_dfr(~tibble(Package = .x$Package, Version = .x$Version)) %>% 
    arrange(Package)

long <- readxl::read_excel("Caspian data_01.03.2023_SA.xlsx", sheet = "main") %>% 
    select(O:SmSdTu5) %>% 
    mutate_all(as.character) %>% 
    pivot_longer(names_to = "id", values_to = "v", -c("O", "sp")) %>% 
    separate(v, into = c("ad", "jv"), sep = "\\+", convert = TRUE) %>% # 
    pivot_longer(names_to = "age", values_to = "abu", -c(1:3)) %>% 
    mutate(id = substr(id, 3, 7)) %>%  # seria = substr(id, 1, 4)) %>%  #id = str_replace_all(id, "_", "")) %>% 
    filter(!is.na(abu)) %>% # seria != "Sw__") %>% 
    arrange(id)

orwi <- long %>% 
    filter(O == "Oribatida", sp != "Oribatida_juvenile_indet") %>% 
    select(-O, -age) %>% 
    pivot_wider(names_from = id, values_from = abu, values_fn = sum, values_fill = 0)
orws <- long %>% 
    filter(O == "Oribatida", sp != "Oribatida_juvenile_indet") %>% 
    mutate(seria = substr(id, 1, 4)) %>% 
    select(-O, -id, -age) %>% 
    pivot_wider(names_from = seria, values_from = abu, values_fn = sum, values_fill = 0)
mswi <- long %>% 
    filter(O == "Mesostigmata") %>% 
    select(-O, -age) %>%
    pivot_wider(names_from = id, values_from = abu, values_fn = sum, values_fill = 0)
msws <- long %>% 
    filter(O == "Mesostigmata") %>% 
    mutate(seria = substr(id, 1, 4)) %>% 
    select(-O, -id, -age) %>%
    pivot_wider(names_from = seria, values_from = abu, values_fn = sum, values_fill = 0)


labs <- readxl::read_excel("Caspian data_01.03.2023_SA.xlsx", sheet = "samples") %>% 
    filter(distr == "Samoor", code != "SmSw") %>% 
    transmute(id = substr(id, 3, 7), 
              seria = substr(id, 1, 4), 
              plants.d, 
              coast = case_when(coast == "sandy dunes" ~ "dunes", TRUE ~ coast),
              # C, %
              RH) %>% 
    separate(plants.d, sep = " ", into = "plants.d", extra = "drop")
taxa <- readxl::read_excel("Caspian data_01.03.2023_SA.xlsx", sheet = "taxa") %>% 
    select(sp, area2, order, family)

# + Рис. 4. Обилие микроартропод по типам берега и парцеллам  ---------------
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

# + Табл. 2. Корреляция обилия групп друг с другом и RH%, C%. ---------------
# By samples 
cor.data1 <- long %>% 
    filter(sp %in% c("Astigmata_total", "Collembola", "Mesostigmata_total", "Oribatida_total", "Prostigmata")) %>% 
    transmute(sp = str_replace(sp, "_total", ""), 
              id, 
              abu) %>% 
    pivot_wider(names_from = sp, values_from = abu) %>% 
    left_join(select(labs, id, RH), by = "id")

cor.val1 <- cor(cor.data1[,2:ncol(cor.data1)], method = "spearman")
cor.pval1 <- expand_grid(v1 = 2:ncol(cor.data1), v2 = 2:ncol(cor.data1)) %>% 
    split(1:nrow(.)) %>% 
    lapply(function(a){
        p = cor.test(as_vector(cor.data1[,a$v1]), as_vector(cor.data1[,a$v2]), 
                     method = "spearman")
        data.frame(p = p$p.value, 
                   v1 = colnames(cor.data1)[a$v1], 
                   v2 = colnames(cor.data1)[a$v2])
    }) %>% 
    map_df(tibble) %>% 
    mutate(p = p.adjust(p, method = "BY")) %>% 
    pivot_wider(names_from = v2, values_from = p) %>% 
    column_to_rownames("v1") %>% 
    as.matrix()

corrplot::corrplot(
    corr = cor.val1, 
    p.mat = cor.pval1, 
    type="upper", 
    order = "original",
    diag = FALSE,
    col =  corrplot::COL2('RdYlBu', 10)[10:1], 
    sig.level = 0.05)

paste("ρ=", round(cor.val1, 2), "; p=", round(cor.pval1, 2), sep = "") %>% 
    matrix(ncol = 6, byrow = TRUE) %>% 
    `colnames<-`(colnames(cor.val1)) %>% 
    `rownames<-`(rownames(cor.val1)) %>% 
    as.data.frame() %>% 
    rownames_to_column("taxa") 

# + Табл. 3. Таксономический состав панцирных и мезостигматических к --------
long %>% 
    left_join(labs, by = "id") %>% 
    group_by(coast, sp) %>% 
    summarise(abu = sum(abu), .groups = "drop") %>% 
    pivot_wider(names_from = coast, values_from = abu, values_fill = 0) %>% 
    left_join(taxa, by = "sp") %>% 
    select(order, family, species = sp, range = area2, `sandy beach`, pebbly, dunes, reeds)


# + Рис. 5. Куммуляты по типам берега (АС). ---------------------------------
s <- Sys.time()
rar.or.c <- long %>% 
    filter(O == "Oribatida") %>% 
    left_join(labs, by = "id") %>% 
    select(coast, abu, sp) %>% 
    pivot_wider(names_from = coast, values_from = abu, values_fill = 0, values_fn = sum) %>% 
    select(-sp) %>% 
    lapply(function(a){sort(a[a>0], decreasing = TRUE)}) %>% 
    parLapply(cl = cl, ., function(a){
        a |>
            iNEXT::iNEXT( 
                q=0, 
                datatype="abundance", 
                se = TRUE,
                size = seq(0, 2000, by = 5))
    }) %>%
    map_dfr(~.x %>% 
                pluck("iNextEst", "size_based") %>% 
                select(m, Method, qD, qD.LCL, qD.UCL), 
            .id = "coast") %>% 
    mutate(taxa = "Oribatida")
Sys.time() - s

s <- Sys.time()
rar.ms.c <- long %>% 
    filter(O == "Mesostigmata") %>% 
    left_join(labs, by = "id") %>% 
    select(coast, abu, sp) %>% 
    pivot_wider(names_from = coast, values_from = abu, values_fill = 0, values_fn = sum) %>% 
    select(-sp) %>% 
    lapply(function(a){sort(a[a>0], decreasing = TRUE)}) %>% 
    parLapply(cl = cl, ., function(a){
        a |>
            iNEXT::iNEXT( 
                q=0, 
                datatype="abundance", 
                se = TRUE,
                size = seq(0, 800, by = 5))
    }) %>%
    map_dfr(~.x %>% 
                pluck("iNextEst", "size_based") %>% 
                select(m, Method, qD, qD.LCL, qD.UCL), 
            .id = "coast") %>% 
    mutate(taxa = "Mesostigmata")
Sys.time() - s

rar <- rbind(rar.ms.c, rar.or.c) %>% as_tibble()

p5o <- rar %>% 
    ggplot(aes(x = m, y = qD, group = coast, fill = coast, color = coast)) + 
    labs(x = "individuals", y = "Number of species") + 
    facet_wrap(~taxa, scales = "free")
p5b <- p5o + 
    geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.3, color = "transparent")
p5o <- p5o +
    geom_line(data = filter(rar, Method != "Extrapolation")) +
    geom_line(data = filter(rar, Method == "Extrapolation"), linetype = "dashed") + 
    geom_point(data = filter(rar, Method == "Observed"), 
               shape = 15, size = 2)

p5b + 
    geom_line(data = filter(rar, Method != "Extrapolation")) +
    geom_line(data = filter(rar, Method == "Extrapolation"), linetype = "dashed") + 
    geom_point(data = filter(rar, Method == "Observed"), 
               shape = 15, size = 2)
p5o

# + Рис. 6. Кладограмма по фаунистическим спискам отдельных мероценозов --------
dis <- list()
# dissimilarity
dis$or.bin <- orwi %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "jaccard", binary = TRUE)
dis$or.num <- orwi %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "bray", binary = FALSE)
dis$ms.bin <- mswi %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "jaccard", binary = TRUE)
dis$ms.num <- mswi %>% 
    column_to_rownames("sp") %>% 
    select_if(function(a){sum(a)>0}) %>% 
    t %>% 
    as.data.frame() %>% 
    vegan::vegdist(method = "bray", binary = FALSE)

# https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html

dend <- lapply(dis, function(a){ 
    dend <- a %>% 
        hclust(method = 'ward.D2') %>% ### METHOD
        as.dendrogram()
    L <- labels(dend) %>% 
        tibble(id = .) %>% 
        left_join(labs, by = "id") %>% 
        mutate(l = factor(coast), 
               l = as.numeric(l)) 
    dend %>% 
        set("labels_col", L$l) %>% 
        set("labels_cex", 0.5)
    })

par(mfrow = c(2,2))
plot(dend[[1]],
    horiz = TRUE, 
    main = "Order = Oribatida\n Data = binary (Jaccard)\n Method = Ward")
plot(dend[[2]],
     horiz = TRUE, 
     main = "Order = Oribatida\n Data = numeric (Bray-Curtis)\n Method = Ward")
plot(dend[[3]],
     horiz = TRUE, 
     main = "Order = Mesostigmata\n Data = binary (Jaccard)\n Method = Ward")
plot(dend[[4]],
     horiz = TRUE, 
     main = "Order = Mesostigmata\n Data = numeric (Bray-Curtis)\n Method = Ward")
par(mfrow = c(1,1))

# + Рис. 7. Видовое богатство по сериям  ------------------------------------
# в среднем в пробе этой серии
div <- orwi %>% 
    select(-sp) %>% 
    as.list() %>% 
    lapply(function(a){data.frame(Oribatida = length(a[a>0]))}) %>% 
    map_df(rbind, .id = "id")
div <- mswi %>% 
    select(-sp) %>% 
    as.list() %>% 
    lapply(function(a){data.frame(Mesostigmata = length(a[a>0]))}) %>% 
    map_df(rbind, .id = "id") %>% 
    left_join(div, by = "id") %>% 
    inner_join(select(labs, -RH), by = "id") %>% 
    as_tibble()

p7 <- div %>% 
    select(-id) %>% 
    pivot_longer(names_to = "taxa", values_to = "nsp", -c("seria", "coast", "plants.d")) %>% 
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
               fill = taxa, 
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
ggsave("plot_7.pdf", p7, width = 297/25, height = 210/25)

# + Рис. 8. Ареалогический состав по сериям: качественный и количест --------
p8 <- rbind(mutate(mswi, Order = "Mesostigmata", .before = 1), 
            mutate(orwi, Order = "Oribatida", .before = 1), 
            mutate(rbind(mswi, orwi), .before = 1,
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

# + Табл. 5. Виды-специалисты (индикаторы) ----------------------------------

set.seed(2)
iv.s <- orws %>% 
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
# summary(iv.s)
iv.s <- iv.s$str %>% 
    as.data.frame() %>% 
    rownames_to_column("sp") %>% 
    as_tibble() %>% 
    pivot_longer(names_to = "biotop", values_to = "iv", -sp) %>% 
    group_by(sp) %>% 
    filter(iv == max(iv)) %>% 
    ungroup() %>% 
    left_join(select(rownames_to_column(iv.s$sign, "sp"), sp, p.value), 
              by = "sp") %>% 
    mutate(`iv, %_s` = round(iv*100), 
           p.value_s = round(p.value, 4),
           sign_s = case_when(p.value <= 0.001 ~ "***", 
                              p.value <= 0.01 ~ "**", 
                              p.value <= 0.05 ~ "*", 
                              TRUE ~ ""), 
           .keep = "unused", 
           .after = 2) %>% 
    left_join(select(taxa, sp, order), by = "sp") %>% 
    select(order, sp:sign_s)

set.seed(3)
iv.i <- orwi %>% 
    rbind(mswi) %>%
    mutate(total = apply(.[,-1], 1, sum), .before = 2) %>% 
    filter(total >= 5) %>% 
    select(-total) %>% 
    column_to_rownames("sp") %>% 
    t %>%
    indicspecies::multipatt(
        pull(arrange(distinct(select(labs, id, coast)), id) , coast), 
        control = how(nperm=999),
        max.order = 4,
        func = "indval",
        duleg = FALSE
    )
iv.i <- iv.i$str %>% 
    as.data.frame() %>% 
    rownames_to_column("sp") %>% 
    as_tibble() %>% 
    pivot_longer(names_to = "biotop", values_to = "iv", -sp) %>% 
    group_by(sp) %>% 
    filter(iv == max(iv)) %>% 
    ungroup() %>% 
    left_join(select(rownames_to_column(iv.i$sign, "sp"), sp, p.value), 
              by = "sp") %>% 
    mutate(`iv, %_i` = round(iv*100), 
           p.value_i = round(p.value, 4),
           sign_i = case_when(p.value <= 0.001 ~ "***", 
                              p.value <= 0.01 ~ "**", 
                              p.value <= 0.05 ~ "*", 
                              TRUE ~ ""), 
           .keep = "unused", 
           .after = 2) %>% 
    left_join(select(taxa, sp, order), by = "sp") %>% 
    select(order, sp:sign_i)

full_join(iv.s, iv.i, by = c("order", "sp", "biotop")) %>% 
    DT::datatable(., 
                  filter = 'top', 
                  extensions = c('FixedColumns',"FixedHeader"),
                  options = list(scrollX = TRUE, 
                                 paging=FALSE,
                                 fixedHeader=TRUE)) %>% 
    DT::formatStyle(
        'sp', fontStyle = list(fontStyle = 'italic'))
    
    
# Рис. 9. ССА распределения видов -----------------------------------------
# Факторы: тип берега (? 1-4), 
# тип растительности (рогозы, злаки, ситники, хвощ, лох, тростники), 
# расстояние до моря (м), RH%, C%  

"https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/ca-dca-and-cca/"
"https://gist.github.com/perrygeo/7572735"
"https://pmassicotte.github.io/stats-denmark-2019/07_rda.html#/redundancy-analysis"

library(vegan)
library(labdsv)

# + Рис. 10. Ординация мероценозов Mesostigmata и Oribatida  ----------------

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


p10a <- ggplot(M2, aes(x = axis1, y = axis2, color = plants.d)) + 
    geom_point() + 
    stat_ellipse() +
    geom_text(aes(label = axis1, x = 0, y = -0.77), color = "black", data = eig, alpha = 0.68) +
    geom_text(aes(label = axis2, x = -0.99, y = 0), color = "black", data = eig, alpha = 0.68, angle = 90) +
    facet_grid(cols = vars(type), rows = vars(taxa)) + 
    labs(x = NULL, y = NULL, color = NULL, #subtitle = "Б. Доминантные виды растений") + 
         subtitle = "A. Dominant plant species") + 
    theme(legend.text = element_text(face = "italic"))
p10b <- ggplot(M2, aes(x = axis1, y = axis2, color = coast)) + 
    geom_point() + 
    stat_ellipse() +
    geom_text(aes(label = axis1, x = 0, y = -0.63), color = "black", data = eig, alpha = 0.68) +
    geom_text(aes(label = axis2, x = -0.9, y = 0), color = "black", data = eig, alpha = 0.68, angle = 90) +
    scale_color_manual(values = c("#F8766D", "#00A9FF", "#0CB720", "#CD9600"))+
    facet_grid(cols = vars(type), rows = vars(taxa)) + 
    labs(x = NULL, y = NULL, color  = NULL, #subtitle = "A. Тип берега")
         subtitle = "B. Coast type") 
p10 <- gridExtra::grid.arrange(p10a, p10b, ncol = 1) #, left = "Суммарное обилие в серии, экз.")
ggsave("plot_10.pdf", p10, width = 210/25, height = 297/25)




