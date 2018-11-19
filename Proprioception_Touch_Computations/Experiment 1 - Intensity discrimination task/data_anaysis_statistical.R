# Intensity - Finger vs. Finger

library(car) # Required library

# PSE data
DataPSE = read.csv("rmANOVA_data_pse.txt", sep = "")
iDataPSE = read.csv("rmANOVA_iData.txt", sep = "")
dataBindPSE = cbind(DataPSE$NearAd2, DataPSE$NearAd4, DataPSE$NearAd7, DataPSE$MiddleAd2, DataPSE$MiddleAd4, DataPSE$MiddleAd7, DataPSE$FarAd2, DataPSE$FarAd4, DataPSE$FarAd7)

# rmANOVA
lm.pse = lm(dataBindPSE ~ 1)
anova.pse = Anova(lm.pse, idata = iDataPSE, idesign = ~proprioception*distractor, type = "III")
summary(anova.pse, multivariate = FALSE)

# PSE - post-hoc analysis
DataPSE$subj_ID <- NULL
z = stack(DataPSE)
names(z) <- c("level","group")
AV1 = aov(level ~ group, data = z)
tk = TukeyHSD(AV1)
tk
plot(tk)



# JND data
DataJND = read.csv("rmANOVA_data_jnd.txt", sep = "")
iDataJND = read.csv("rmANOVA_iData.txt", sep = "")
dataBindJND = cbind(DataJND$NearAd2, DataJND$NearAd4, DataJND$NearAd7, DataJND$MiddleAd2, DataJND$MiddleAd4, DataJND$MiddleAd7, DataJND$FarAd2, DataJND$FarAd4, DataJND$FarAd7)

# rmANOVA
lm.jnd = lm(dataBindJND ~ 1)
anova.jnd = Anova(lm.jnd, idata = iDataJND, idesign = ~proprioception*distractor, type = "III")
summary(anova.jnd, multivariate = FALSE)

# JND - post-hoc analysis
DataJND$subj_ID <- NULL
z = stack(DataJND)
names(z) <- c("level","group")
AV1 = aov(level ~ group, data = z)
tk = TukeyHSD(AV1)
tk
plot(tk)
