ms_prs_pval_001_scores <- scan("~/rezultatai/ms_prs_pval_001_scores.txt")
print("Filtravimo su p reikšme 0.01 įverčiai:")
ms_prs_pval_001_scores

print("Vidurkis:")
mean(ms_prs_pval_001_scores)
print("Mediana:")
median(ms_prs_pval_001_scores)
print("Minimumas:")
min(ms_prs_pval_001_scores)
print("Maksimumas:")
max(ms_prs_pval_001_scores)
print("Standartinis nuokrypis:")
sd(ms_prs_pval_001_scores)

ms_prs_pval_001_se <- sd(ms_prs_pval_001_scores) / sqrt(length(ms_prs_pval_001_scores))
print("Standartinė paklaida:")
ms_prs_pval_001_se

z_ms_prs_pval_001_scores <- as.numeric(scale(ms_prs_pval_001_scores))
print("Normalizuoti filtravimo su p reikšme 0.01 įverčiai:")
z_ms_prs_pval_001_scores

print("Normalizuotų įverčių vidurkis:")
median(z_ms_prs_pval_001_scores)
print("Normalizuotų įverčių minimumas:")
min(z_ms_prs_pval_001_scores)
print("Normalizuotų įverčių maksimumas:")
max(z_ms_prs_pval_001_scores)

z_ms_prs_pval_001_se <- sd(z_ms_prs_pval_001_scores) / sqrt(length(z_ms_prs_pval_001_scores))
print("Normalizuotų įverčių standartinė paklaida:")
z_ms_prs_pval_001_se

ms_prs_pval_005_scores <- scan("~/rezultatai/ms_prs_pval_005_scores.txt")
print("Filtravimo su p reikšme 0.05 įverčiai:")
ms_prs_pval_005_scores

print("Vidurkis:")
mean(ms_prs_pval_005_scores)
print("Mediana:")
median(ms_prs_pval_005_scores)
print("Minimumas:")
min(ms_prs_pval_005_scores)
print("Maksimumas:")
max(ms_prs_pval_005_scores)
print("Standartinis nuokrypis:")
sd(ms_prs_pval_005_scores)

ms_prs_pval_005_se <- sd(ms_prs_pval_005_scores) / sqrt(length(ms_prs_pval_005_scores))
print("Standartinė paklaida:")
ms_prs_pval_005_se

z_ms_prs_pval_005_scores <- as.numeric(scale(ms_prs_pval_005_scores))
print("Normalizuoti filtravimo su p reikšme 0.05 įverčiai:")
z_ms_prs_pval_005_scores

print("Normalizuotų įverčių vidurkis:")
median(z_ms_prs_pval_005_scores)
print("Normalizuotų įverčių minimumas:")
min(z_ms_prs_pval_005_scores)
print("Normalizuotų įverčių maksimumas:")
max(z_ms_prs_pval_005_scores)

z_ms_prs_pval_005_se <- sd(z_ms_prs_pval_005_scores) / sqrt(length(z_ms_prs_pval_005_scores))
print("Normalizuotų įverčių standartinė paklaida:")
z_ms_prs_pval_005_se

ms_prs_pval_010_scores <- scan("~/rezultatai/ms_prs_pval_010_scores.txt")
print("Filtravimo su p reikšme 0.1 įverčiai:")
ms_prs_pval_010_scores

print("Vidurkis:")
mean(ms_prs_pval_010_scores)
print("Mediana:")
median(ms_prs_pval_010_scores)
print("Minimumas:")
min(ms_prs_pval_010_scores)
print("Maksimumas:")
max(ms_prs_pval_010_scores)
print("Standartinis nuokrypis:")
sd(ms_prs_pval_010_scores)

ms_prs_pval_010_se <- sd(ms_prs_pval_010_scores) / sqrt(length(ms_prs_pval_010_scores))
print("Standartinė paklaida:")
ms_prs_pval_010_se

z_ms_prs_pval_010_scores <- as.numeric(scale(ms_prs_pval_010_scores))
print("Normalizuoti filtravimo su p reikšme 0.1 įverčiai:")
z_ms_prs_pval_010_scores

print("Normalizuotų įverčių vidurkis:")
median(z_ms_prs_pval_010_scores)
print("Normalizuotų įverčių minimumas:")
min(z_ms_prs_pval_010_scores)
print("Normalizuotų įverčių maksimumas:")
max(z_ms_prs_pval_010_scores)

z_ms_prs_pval_010_se <- sd(z_ms_prs_pval_010_scores) / sqrt(length(z_ms_prs_pval_010_scores))
print("Normalizuotų įverčių standartinė paklaida:")
z_ms_prs_pval_010_se

ms_prs_pval_050_scores <- scan("~/rezultatai/ms_prs_pval_050_scores.txt")
print("Filtravimo su p reikšme 0.5 įverčiai:")
ms_prs_pval_050_scores

print("Vidurkis:")
mean(ms_prs_pval_050_scores)
print("Mediana:")
median(ms_prs_pval_050_scores)
print("Minimumas:")
min(ms_prs_pval_050_scores)
print("Maksimumas:")
max(ms_prs_pval_050_scores)
print("Standartinis nuokrypis:")
sd(ms_prs_pval_050_scores)

ms_prs_pval_050_se <- sd(ms_prs_pval_050_scores) / sqrt(length(ms_prs_pval_050_scores))
print("Standartinė paklaida:")
ms_prs_pval_050_se

z_ms_prs_pval_050_scores <- as.numeric(scale(ms_prs_pval_050_scores))
print("Normalizuoti filtravimo su p reikšme 0.5 įverčiai:")
z_ms_prs_pval_050_scores

print("Normalizuotų įverčių vidurkis:")
median(z_ms_prs_pval_050_scores)
print("Normalizuotų įverčių minimumas:")
min(z_ms_prs_pval_050_scores)
print("Normalizuotų įverčių maksimumas:")
max(z_ms_prs_pval_050_scores)

z_ms_prs_pval_050_se <- sd(z_ms_prs_pval_050_scores) / sqrt(length(z_ms_prs_pval_050_scores))
print("Normalizuotų įverčių standartinė paklaida:")
z_ms_prs_pval_050_se

library(ggplot2)

z_scores_pval <- data.frame(
  z_ms_prs_pval_001_scores = z_ms_prs_pval_001_scores,
  z_ms_prs_pval_005_scores = z_ms_prs_pval_005_scores,
  z_ms_prs_pval_010_scores = z_ms_prs_pval_010_scores,
  z_ms_prs_pval_050_scores = z_ms_prs_pval_050_scores
)

z_scores_pval$number <- 1:nrow(z_scores_pval)

z_scores_pval$number <- factor(z_scores_pval$number, levels = z_scores_pval$number)
z_scores_pval

z_scores_pval_long <- reshape(z_scores_pval,
                              direction = "long",
                              idvar = "number",
                              varying = list(names(z_scores_pval)[1:4]),
                              v.names = "scores",
                              timevar = "p_val",
                              times = c("p_val_001","p_val_005","p_val_010","p_val_050"))
z_scores_pval_long

z_scores_pval_graph <- ggplot(z_scores_pval_long, aes(x = factor(number), y = scores, fill = p_val)) +
                         geom_col(position = "stack") +
                         coord_flip() +
                         labs(x = "Tiriamieji",
                              y = "Standartinio nuokrypio vienetai",
                              title = "Normalizuoti genetinės rizikos įverčiai") +
                         theme_minimal() +
                         theme(panel.grid.major.y = element_blank(),
                               axis.text.y = element_blank(),
                               plot.title = element_text(size = 18, hjust = 0.5, margin = margin(t = 20, b = 20)),
                               axis.title = element_text(size = 14, hjust = 0.5, margin = margin(t = 20, b = 20)),
                               axis.title.x = element_text(margin = margin(t = 20, b = 20), face = "italic"),
                               axis.title.y = element_text(margin = margin(t = 20, b = 20), face = "italic"),
                               legend.position = "bottom") +
                         scale_fill_discrete(labels = c("0.01", "0.05", "0.1", "0.5")) +
                         labs(fill = "p-reikšmė")

ggsave("~/grafikai/p_val_scores.png", plot = z_scores_pval_graph, bg = "white")

ms_prs_pval_001_LD_010_scores <- scan("~/rezultatai/ms_prs_pval_001_LD_010_scores.txt")
print("Filtravimo su p reikšme 0.01 ir r2 0.1 įverčiai:")
ms_prs_pval_001_LD_010_scores

print("Vidurkis:")
mean(ms_prs_pval_001_LD_010_scores)
print("Mediana:")
median(ms_prs_pval_001_LD_010_scores)
print("Minimumas:")
min(ms_prs_pval_001_LD_010_scores)
print("Maksimumas:")
max(ms_prs_pval_001_LD_010_scores)
print("Standartinis nuokrypis:")
sd(ms_prs_pval_001_LD_010_scores)

ms_prs_pval_001_LD_010_se <- sd(ms_prs_pval_001_LD_010_scores) / sqrt(length(ms_prs_pval_001_LD_010_scores))
print("Standartinė paklaida:")
ms_prs_pval_001_LD_010_se

z_ms_prs_pval_001_LD_010_scores <- as.numeric(scale(ms_prs_pval_001_LD_010_scores))
print("Normalizuoti filtravimo su p reikšme 0.01 ir r2 0.1 įverčiai:")
z_ms_prs_pval_001_LD_010_scores

print("Normalizuotų įverčių vidurkis:")
median(z_ms_prs_pval_001_LD_010_scores)
print("Normalizuotų įverčių minimumas:")
min(z_ms_prs_pval_001_LD_010_scores)
print("Normalizuotų įverčių maksimumas:")
max(z_ms_prs_pval_001_LD_010_scores)

z_ms_prs_pval_001_LD_010_se <- sd(z_ms_prs_pval_001_LD_010_scores) / sqrt(length(z_ms_prs_pval_001_LD_010_scores))
print("Normalizuotų įverčių standartinė paklaida:")
z_ms_prs_pval_001_LD_010_se

ms_prs_pval_005_LD_010_scores <- scan("~/rezultatai/ms_prs_pval_005_LD_010_scores.txt")
print("Filtravimo su p reikšme 0.05 ir r2 0.1 įverčiai:")
ms_prs_pval_005_LD_010_scores

print("Vidurkis:")
mean(ms_prs_pval_005_LD_010_scores)
print("Mediana:")
median(ms_prs_pval_005_LD_010_scores)
print("Minimumas:")
min(ms_prs_pval_005_LD_010_scores)
print("Maksimumas:")
max(ms_prs_pval_005_LD_010_scores)
print("Standartinis nuokrypis:")
sd(ms_prs_pval_005_LD_010_scores)

ms_prs_pval_005_LD_010_se <- sd(ms_prs_pval_005_LD_010_scores) / sqrt(length(ms_prs_pval_005_LD_010_scores))
print("Standartinė paklaida:")
ms_prs_pval_005_LD_010_se

z_ms_prs_pval_005_LD_010_scores <- as.numeric(scale(ms_prs_pval_005_LD_010_scores))
print("Normalizuoti filtravimo su p reikšme 0.05 ir r2 0.1 įverčiai:")
z_ms_prs_pval_005_LD_010_scores

print("Normalizuotų įverčių vidurkis:")
median(z_ms_prs_pval_005_LD_010_scores)
print("Normalizuotų įverčių minimumas:")
min(z_ms_prs_pval_005_LD_010_scores)
print("Normalizuotų įverčių maksimumas:")
max(z_ms_prs_pval_005_LD_010_scores)

z_ms_prs_pval_005_LD_010_se <- sd(z_ms_prs_pval_005_LD_010_scores) / sqrt(length(z_ms_prs_pval_005_LD_010_scores))
print("Normalizuotų įverčių standartinė paklaida:")
z_ms_prs_pval_005_LD_010_se

ms_prs_pval_010_LD_010_scores <- scan("~/rezultatai/ms_prs_pval_010_LD_010_scores.txt")
print("Filtravimo su p reikšme 0.1 ir r2 0.1 įverčiai:")
ms_prs_pval_010_LD_010_scores

print("Vidurkis:")
mean(ms_prs_pval_010_LD_010_scores)
print("Mediana:")
median(ms_prs_pval_010_LD_010_scores)
print("Minimumas:")
min(ms_prs_pval_010_LD_010_scores)
print("Maksimumas:")
max(ms_prs_pval_010_LD_010_scores)
print("Standartinis nuokrypis:")
sd(ms_prs_pval_010_LD_010_scores)

ms_prs_pval_010_LD_010_se <- sd(ms_prs_pval_010_LD_010_scores) / sqrt(length(ms_prs_pval_010_LD_010_scores))
print("Standartinė paklaida:")
ms_prs_pval_010_LD_010_se

z_ms_prs_pval_010_LD_010_scores <- as.numeric(scale(ms_prs_pval_010_LD_010_scores))
print("Normalizuoti filtravimo su p reikšme 0.1 ir r2 0.1 įverčiai:")
z_ms_prs_pval_010_LD_010_scores

print("Normalizuotų įverčių vidurkis:")
median(z_ms_prs_pval_010_LD_010_scores)
print("Normalizuotų įverčių minimumas:")
min(z_ms_prs_pval_010_LD_010_scores)
print("Normalizuotų įverčių maksimumas:")
max(z_ms_prs_pval_010_LD_010_scores)

z_ms_prs_pval_010_LD_010_se <- sd(z_ms_prs_pval_010_LD_010_scores) / sqrt(length(z_ms_prs_pval_010_LD_010_scores))
print("Normalizuotų įverčių standartinė paklaida:")
z_ms_prs_pval_010_LD_010_se

ms_prs_pval_050_LD_010_scores <- scan("~/rezultatai/ms_prs_pval_050_LD_010_scores.txt")
print("Filtravimo su p reikšme 0.5 ir r2 0.1 įverčiai:")
ms_prs_pval_050_LD_010_scores

print("Vidurkis:")
mean(ms_prs_pval_050_LD_010_scores)
print("Mediana:")
median(ms_prs_pval_050_LD_010_scores)
print("Minimumas:")
min(ms_prs_pval_050_LD_010_scores)
print("Maksimumas:")
max(ms_prs_pval_050_LD_010_scores)
print("Standartinis nuokrypis:")
sd(ms_prs_pval_050_LD_010_scores)

ms_prs_pval_050_LD_010_se <- sd(ms_prs_pval_050_LD_010_scores) / sqrt(length(ms_prs_pval_050_LD_010_scores))
print("Standartinė paklaida:")
ms_prs_pval_050_LD_010_se

z_ms_prs_pval_050_LD_010_scores <- as.numeric(scale(ms_prs_pval_050_LD_010_scores))
print("Normalizuoti filtravimo su p reikšme 0.5 ir r2 0.1 įverčiai:")
z_ms_prs_pval_050_LD_010_scores

print("Normalizuotų įverčių vidurkis:")
median(z_ms_prs_pval_050_LD_010_scores)
print("Normalizuotų įverčių minimumas:")
min(z_ms_prs_pval_050_LD_010_scores)
print("Normalizuotų įverčių maksimumas:")
max(z_ms_prs_pval_050_LD_010_scores)

z_ms_prs_pval_050_LD_010_se <- sd(z_ms_prs_pval_050_LD_010_scores) / sqrt(length(z_ms_prs_pval_050_LD_010_scores))
print("Normalizuotų įverčių standartinė paklaida:")
z_ms_prs_pval_050_LD_010_se

ms_prs_LD_010_scores <- scan("~/rezultatai/ms_prs_LD_010_scores.txt")
print("Filtravimo su r2 0.1 įverčiai:")
ms_prs_LD_010_scores

print("Vidurkis:")
mean(ms_prs_LD_010_scores)
print("Mediana:")
median(ms_prs_LD_010_scores)
print("Minimumas:")
min(ms_prs_LD_010_scores)
print("Maksimumas:")
max(ms_prs_LD_010_scores)
print("Standartinis nuokrypis:")
sd(ms_prs_LD_010_scores)

ms_prs_LD_010_se <- sd(ms_prs_LD_010_scores) / sqrt(length(ms_prs_LD_010_scores))
print("Standartinė paklaida:")
ms_prs_LD_010_se

z_ms_prs_LD_010_scores <- as.numeric(scale(ms_prs_LD_010_scores))
print("Normalizuoti filtravimo su r2 0.1 įverčiai:")
z_ms_prs_LD_010_scores

print("Normalizuotų įverčių vidurkis:")
median(z_ms_prs_LD_010_scores)
print("Normalizuotų įverčių minimumas:")
min(z_ms_prs_LD_010_scores)
print("Normalizuotų įverčių maksimumas:")
max(z_ms_prs_LD_010_scores)

z_ms_prs_LD_010_se <- sd(z_ms_prs_LD_010_scores) / sqrt(length(z_ms_prs_LD_010_scores))
print("Normalizuotų įverčių standartinė paklaida:")
z_ms_prs_LD_010_se

z_scores_LD <- data.frame(
  z_ms_prs_pval_001_LD_010_scores = z_ms_prs_pval_001_LD_010_scores,
  z_ms_prs_pval_005_LD_010_scores = z_ms_prs_pval_005_LD_010_scores,
  z_ms_prs_pval_010_LD_010_scores = z_ms_prs_pval_010_LD_010_scores,
  z_ms_prs_pval_050_LD_010_scores = z_ms_prs_pval_050_LD_010_scores,
  z_ms_prs_LD_010_scores = z_ms_prs_LD_010_scores
)

z_scores_LD$number <- 1:nrow(z_scores_LD)

z_scores_LD$number <- factor(z_scores_LD$number, levels = z_scores_LD$number)
z_scores_LD

z_scores_LD_long <- reshape(z_scores_LD,
                              direction = "long",
                              idvar = "number",
                              varying = list(names(z_scores_LD)[1:5]),
                              v.names = "scores",
                              timevar = "p_val",
                              times = c("p_val_001_LD_010","p_val_005_LD_010","p_val_010_LD_010","p_val_050_LD_010", "LD_010"))
z_scores_LD_long

z_scores_LD_graph <- ggplot(z_scores_LD_long, aes(x = factor(number), y = scores, fill = p_val)) +
                         geom_col(position = "stack") +
                         coord_flip() +
                         labs(x = "Tiriamieji",
                              y = "Standartinio nuokrypio vienetai",
                              title = "Normalizuoti genetinės rizikos įverčiai") +
                         theme_minimal() +
                         theme(panel.grid.major.y = element_blank(),
                               axis.text.y = element_blank(),
                               plot.title = element_text(size = 18, hjust = 0.5, margin = margin(t = 20, b = 20)),
                               axis.title = element_text(size = 14, hjust = 0.5, margin = margin(t = 20, b = 20)),
                               axis.title.x = element_text(margin = margin(t = 20, b = 20), face = "italic"),
                               axis.title.y = element_text(margin = margin(t = 20, b = 20), face = "italic"),
                               legend.position = "bottom") +
                         scale_fill_discrete(labels = c("0.01", "0.05", "0.1", "0.5", "visi VNP")) +
                         labs(fill = "p-reikšmė, r2 = 0.1")

ggsave("~/grafikai/LD_scores.png", plot = z_scores_LD_graph, bg = "white")

ms_prs_pval_010_shrink_050_scores <- scan("~/rezultatai/ms_prs_010_shrink_050_scores.txt")
print("Filtravimo su p-reikšme 0.1 ir shrink parametru 0.5 įverčiai:")
ms_prs_pval_010_shrink_050_scores

print("Vidurkis:")
mean(ms_prs_pval_010_shrink_050_scores)
print("Mediana:")
median(ms_prs_pval_010_shrink_050_scores)
print("Minimumas:")
min(ms_prs_pval_010_shrink_050_scores)
print("Maksimumas:")
max(ms_prs_pval_010_shrink_050_scores)
print("Standartinis nuokrypis:")
sd(ms_prs_pval_010_shrink_050_scores)

ms_prs_pval_010_shrink_050_se <- sd(ms_prs_pval_010_shrink_050_scores) / sqrt(length(ms_prs_pval_010_shrink_050_scores))
print("Standartinė paklaida:")
ms_prs_pval_010_shrink_050_se

z_ms_prs_pval_010_shrink_050_scores <- as.numeric(scale(ms_prs_pval_010_shrink_050_scores))
print("Normalizuoti filtravimo su p-reikšme 0.1 ir shrink parametru 0.5 įverčiai:")
z_ms_prs_pval_010_shrink_050_scores

print("Normalizuotų įverčių vidurkis:")
median(z_ms_prs_pval_010_shrink_050_scores)
print("Normalizuotų įverčių minimumas:")
min(z_ms_prs_pval_010_shrink_050_scores)
print("Normalizuotų įverčių maksimumas:")
max(z_ms_prs_pval_010_shrink_050_scores)

z_ms_prs_pval_010_shrink_050_se <- sd(z_ms_prs_pval_010_shrink_050_scores) / sqrt(length(z_ms_prs_pval_010_shrink_050_scores))
print("Normalizuotų įverčių standartinė paklaida:")
z_ms_prs_pval_010_shrink_050_se

ms_prs_pval_010_shrink_200_scores <- scan("~/rezultatai/ms_prs_010_shrink_200_scores.txt")
print("Filtravimo su p-reikšme 0.1 ir shrink parametru 2.0 įverčiai:")
ms_prs_pval_010_shrink_200_scores

print("Vidurkis:")
mean(ms_prs_pval_010_shrink_200_scores)
print("Mediana:")
median(ms_prs_pval_010_shrink_200_scores)
print("Minimumas:")
min(ms_prs_pval_010_shrink_200_scores)
print("Maksimumas:")
max(ms_prs_pval_010_shrink_200_scores)
print("Standartinis nuokrypis:")
sd(ms_prs_pval_010_shrink_200_scores)

ms_prs_pval_010_shrink_200_se <- sd(ms_prs_pval_010_shrink_200_scores) / sqrt(length(ms_prs_pval_010_shrink_200_scores))
print("Standartinė paklaida:")
ms_prs_pval_010_shrink_200_se

z_ms_prs_pval_010_shrink_200_scores <- as.numeric(scale(ms_prs_pval_010_shrink_200_scores))
print("Normalizuoti filtravimo su p-reikšme 0.1 ir shrink parametru 2.0 įverčiai:")
z_ms_prs_pval_010_shrink_200_scores

print("Normalizuotų įverčių vidurkis:")
median(z_ms_prs_pval_010_shrink_200_scores)
print("Normalizuotų įverčių minimumas:")
min(z_ms_prs_pval_010_shrink_200_scores)
print("Normalizuotų įverčių maksimumas:")
max(z_ms_prs_pval_010_shrink_200_scores)

z_ms_prs_pval_010_shrink_200_se <- sd(z_ms_prs_pval_010_shrink_200_scores) / sqrt(length(z_ms_prs_pval_010_shrink_200_scores))
print("Normalizuotų įverčių standartinė paklaida:")
z_ms_prs_pval_010_shrink_200_se

ms_prs_pval_005_shrink_050_scores <- scan("~/rezultatai/ms_prs_005_shrink_050_scores.txt")
print("Filtravimo su p-reikšme 0.05 ir shrink parametru 0.5 įverčiai:")
ms_prs_pval_005_shrink_050_scores

print("Vidurkis:")
mean(ms_prs_pval_005_shrink_050_scores)
print("Mediana:")
median(ms_prs_pval_005_shrink_050_scores)
print("Minimumas:")
min(ms_prs_pval_005_shrink_050_scores)
print("Maksimumas:")
max(ms_prs_pval_005_shrink_050_scores)
print("Standartinis nuokrypis:")
sd(ms_prs_pval_005_shrink_050_scores)

ms_prs_pval_005_shrink_050_se <- sd(ms_prs_pval_005_shrink_050_scores) / sqrt(length(ms_prs_pval_005_shrink_050_scores))
print("Standartinė paklaida:")
ms_prs_pval_005_shrink_050_se

z_ms_prs_pval_005_shrink_050_scores <- as.numeric(scale(ms_prs_pval_005_shrink_050_scores))
print("Normalizuoti filtravimo su p-reikšme 0.05 ir shrink parametru 0.5 įverčiai:")
z_ms_prs_pval_005_shrink_050_scores

print("Normalizuotų įverčių vidurkis:")
median(z_ms_prs_pval_005_shrink_050_scores)
print("Normalizuotų įverčių minimumas:")
min(z_ms_prs_pval_005_shrink_050_scores)
print("Normalizuotų įverčių maksimumas:")
max(z_ms_prs_pval_005_shrink_050_scores)

z_ms_prs_pval_005_shrink_050_se <- sd(z_ms_prs_pval_005_shrink_050_scores) / sqrt(length(z_ms_prs_pval_005_shrink_050_scores))
print("Normalizuotų įverčių standartinė paklaida:")
z_ms_prs_pval_005_shrink_050_se

ms_prs_pval_005_shrink_200_scores <- scan("~/rezultatai/ms_prs_005_shrink_200_scores.txt")
print("Filtravimo su p-reikšme 0.05 ir shrink parametru 2.0 įverčiai:")
ms_prs_pval_005_shrink_200_scores

print("Vidurkis:")
mean(ms_prs_pval_005_shrink_200_scores)
print("Mediana:")
median(ms_prs_pval_005_shrink_200_scores)
print("Minimumas:")
min(ms_prs_pval_005_shrink_200_scores)
print("Maksimumas:")
max(ms_prs_pval_005_shrink_200_scores)
print("Standartinis nuokrypis:")
sd(ms_prs_pval_005_shrink_200_scores)

ms_prs_pval_005_shrink_200_se <- sd(ms_prs_pval_005_shrink_200_scores) / sqrt(length(ms_prs_pval_005_shrink_200_scores))
print("Standartinė paklaida:")
ms_prs_pval_005_shrink_200_se

z_ms_prs_pval_005_shrink_200_scores <- as.numeric(scale(ms_prs_pval_005_shrink_200_scores))
print("Normalizuoti filtravimo su p-reikšme 0.05 ir shrink parametru 2.0 įverčiai:")
z_ms_prs_pval_005_shrink_200_scores

print("Normalizuotų įverčių vidurkis:")
median(z_ms_prs_pval_005_shrink_200_scores)
print("Normalizuotų įverčių minimumas:")
min(z_ms_prs_pval_005_shrink_200_scores)
print("Normalizuotų įverčių maksimumas:")
max(z_ms_prs_pval_005_shrink_200_scores)

z_ms_prs_pval_005_shrink_200_se <- sd(z_ms_prs_pval_005_shrink_200_scores) / sqrt(length(z_ms_prs_pval_005_shrink_200_scores))
print("Normalizuotų įverčių standartinė paklaida:")
z_ms_prs_pval_005_shrink_200_se


