# propcut no distractor

library(car) # Required library

# PSE data
DataPSE = read.csv("rmANOVA_data_pse.txt", sep = "")
iDataPSE = read.csv("rmANOVA_iData.txt", sep = "")
dataBindPSE = cbind(DataPSE$near, DataPSE$middle, DataPSE$far)

# rmANOVA
lm.pse = lm(dataBindPSE ~ 1)
anova.pse = Anova(lm.pse, idata = iDataPSE, idesign = ~position, type = "III")
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
dataBindJND = cbind(DataJND$near, DataJND$middle, DataJND$far)

# rmANOVA
lm.jnd = lm(dataBindJND ~ 1)
anova.jnd = Anova(lm.jnd, idata = iDataJND, idesign = ~position, type = "III")
summary(anova.jnd, multivariate = FALSE)

# JND - post-hoc analysis
DataJND$subj_ID <- NULL
z = stack(DataJND)
names(z) <- c("level","group")
AV1 = aov(level ~ group, data = z)
tk = TukeyHSD(AV1)
tk
plot(tk)
