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



