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
