# Finger vs. Finger, Main experiment

library(car) # Required library

# PSE data
DataPSE = read.csv("rmANOVA_data_pse.txt", sep = "")
iDataPSE = read.csv("rmANOVA_iData.txt", sep = "")
dataBindPSE = cbind(DataPSE$NearFd100, DataPSE$NearFd200, DataPSE$NearFd300, DataPSE$MiddleFd100, DataPSE$MiddleFd200, DataPSE$MiddleFd300, DataPSE$FarFd100, DataPSE$FarFd200, DataPSE$FarFd300)

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
dataBindJND = cbind(DataJND$NearFd100, DataJND$NearFd200, DataJND$NearFd300, DataJND$MiddleFd100, DataJND$MiddleFd200, DataJND$MiddleFd300, DataJND$FarFd100, DataJND$FarFd200, DataJND$FarFd300)

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

