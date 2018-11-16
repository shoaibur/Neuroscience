# Bayes Factor Analysis

# Load required library
library(BayesFactor)

# Load Log-Likelihood data
ll = read.csv("bayesFactorAnalysisData.csv", sep = ",")
# Description of the data:
# Column 1 contains LL of the Best model in frequency domain.
# Columns 2-20 contains LL of frequency combination models.
# Column 22 contains LL of the Best model in intensity domain.
# Columns 23-24 contains LL of intensity combination models.
# We will compare LL of other models with the LL of the Best model in both domains.

# Bayes Factor (BF) of frequency combination models (column indices 2-20)
for(i in 2:20){
  # Compute BF based on 2-sample t-test
  out = ttestBF(ll[,i], ll[,1])
  # Extract the BF from stats and print
  BF = exp(1)^out@bayesFactor$bf
  print(BF)
}

# Bayes Factor (BF) of intensity combination models (column indices 23-24)
for(i in 23:24){
  # Compute BF based on 2-sample t-test
  out = ttestBF(-ll[,i], -ll[,22])
  # Extract the BF from stats and print
  BF = exp(1)^out@bayesFactor$bf
  print(BF)
}


