# Finger vs. Forearm

library(car) # Required library

# PSE data
DataPSE = read.csv("rmANOVA_data_pse.txt", sep = "")
iData = read.csv("rmANOVA_iData.txt", sep = "")
dataBind = cbind(DataPSE$NearFd100, DataPSE$NearFd300, DataPSE$FarFd100, DataPSE$FarFd300)

# rmANOVA
lm_fit = lm(dataBind ~ 1)
anova_fit = Anova(lm_fit, idata = iData, idesign = ~proprioception*distractor, type = "III")
summary(anova_fit, multivariate = FALSE)

# PSE - post-hoc analysis
DataPSE$subj_ID <- NULL
z = stack(Data)
names(z) <- c("level","group")
AV1 = aov(level ~ group, data = z)
tk = TukeyHSD(AV1)
tk
plot(tk)


# JND data
DataJND = read.csv("rmANOVA_data_jnd.txt", sep = "")
iData = read.csv("rmANOVA_iData.txt", sep = "")
dataBind = cbind(DataJND$NearFd100, DataJND$NearFd300, DataJND$FarFd100, DataJND$FarFd300)

# rmANOVA
lm_fit = lm(dataBind ~ 1)
anova_fit = Anova(lm_fit, idata = iData, idesign = ~proprioception*distractor, type = "III")
summary(anova_fit, multivariate = FALSE)

# JND - post-hoc analysis
DataJND$subj_ID <- NULL
z = stack(Data)
names(z) <- c("level","group")
AV1 = aov(level ~ group, data = z)
tk = TukeyHSD(AV1)
tk
plot(tk)
