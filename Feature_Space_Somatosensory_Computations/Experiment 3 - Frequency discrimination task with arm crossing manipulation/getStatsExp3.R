# Xed vs. unXed

library(car) # Required library
# PSE data
DataPSE = read.csv("rmANOVA_data_pse.txt", sep = "")
iData = read.csv("rmANOVA_iData.txt", sep = "")
dataBind = cbind(DataPSE$XedFd100, DataPSE$XedFd300, DataPSE$unXedFd100, DataPSE$unXedFd300)
# dataBind = cbind(Data$XedFd100, Data$XedFd300, Data$XedFd0, Data$unXedFd100, Data$unXedFd300, Data$unXedFd0)

# rmANOVA
lm_fit = lm(dataBind ~ 1)
anova_fit = Anova(lm_fit, idata = iData, idesign = ~proprioception*distractor, type = "III")
summary(anova_fit, multivariate = FALSE)

# PSE - post-hoc analysis
Data$subj_ID <- NULL
z = stack(DataPSE)
names(z) <- c("level","group")
AV1 = aov(level ~ group, data = z)
tk = TukeyHSD(AV1)
tk
plot(tk)



# JND data
DataJND = read.csv("rmANOVA_data_jnd.txt", sep = "")
iData = read.csv("rmANOVA_iData.txt", sep = "")
dataBind = cbind(DataJND$XedFd100, DataJND$XedFd300, DataJND$unXedFd100, DataJND$unXedFd300)
# dataBind = cbind(Data$XedFd100, Data$XedFd300, Data$XedFd0, Data$unXedFd100, Data$unXedFd300, Data$unXedFd0)

# rmANOVA
lm_fit = lm(dataBind ~ 1)
anova_fit = Anova(lm_fit, idata = iData, idesign = ~proprioception*distractor, type = "III")
summary(anova_fit, multivariate = FALSE)

# JND - post-hoc analysis
Data$subj_ID <- NULL
z = stack(DataJND)
names(z) <- c("level","group")
AV1 = aov(level ~ group, data = z)
tk = TukeyHSD(AV1)
tk
plot(tk)

