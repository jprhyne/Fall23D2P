library(spdep)
library(sf)
library(smerc)
library(RColorBrewer)


# Load in the file that we saved the jobloss data to
load("COMap.rda")
load("pop.rda")

# Values we use throughout this file
coordsPoly = st_centroid(COMap$geometry)
coordsXY   = st_coordinates(st_centroid(COMap$geometry))
cases = COMap$Total

# Significance value. We are choosing alpha = .01.
# Originally chose alpha = .05 but this gave us a lot of clusters which makes 
# talking about them more difficult, and a smaller value allows us to be more picky
alpha = .01

# Construct neighbor data
comNb=poly2nb(COMap)
plot(COMap$geometry)
plot(comNb, coords=coordsPoly, add=TRUE)
# CEPP method
# fivenum(pop) gives: 779.0   5566.5  15023.0  46684.5 740559.0
# So we pick n* values that will cover many small counties as well 
# as values  that will cover larger ones. Max value is close to half the total
# population
nStarVec = c(800, 1500, 4000, 1e+04, 5e+04, 1e+05, 5e+05, 8e+05, 1e+06, 3e+06)
# If the file for storing our results already exists, we remove it first
if (file.exists("CEPPResults.txt"))
  file.remove("CEPPResults.txt")
for (nStar in nStarVec) {
  temp = cepp.test(coords = coordsXY,
                   cases = cases,
                   pop = pop,
                   nstar = nStar,
                   alpha = alpha)
  plot(COMap$geometry, border = "grey", axes= TRUE,
       col=color.clusters(temp))
  title(main = paste("CEPP with n* =",nStar))
  # Write these results to a file for later reference
  
  sink(file="CEPPResults.txt", append = TRUE)
  cat("--------------------------------------------------------------------------------------------\n")
  cat(paste("Summary for CEPP with n* =", nStar))
  cat("\n")
  cat("--------------------------------------------------------------------------------------------\n")
  print(summary(temp))
  cat("--------------------------------------------------------------------------------------------\n\n")
  sink()
}

# Now we do the Besag-Newell method with varying c* values. 
# fivenum(COMap$Total) produces: 11.0    58.5   210.5   627.0 10404.0
# So we will pick c* values that will cover ranges of these values

cStarVec = c(40, 80, 150, 500, 1000, 5000, 7000, 10000, 15000, 20000, 30000, 35000)
# If the file for storing our results already exists, we remove it first
if (file.exists("BNResults.txt"))
  file.remove("BNResults.txt")

for (cStar in cStarVec) {
  # Perform Besag-Newell test with c* = cStar
  temp = bn.test(coords = coordsXY,
                   cases = cases,
                   pop = pop,
                   cstar = cStar,
                   alpha = alpha)
  if (temp$clusters[[1]]$pvalue > alpha) {
    # skip over this one as we found no significant clusters
    # but first, we print that out to file and console
    # for user information
    print(paste("No significant cluster found for c* =",cStar))
    sink(file="BNResults.txt", append = TRUE)
    cat("--------------------------------------------------------------------------------------------\n")
    cat(paste("Summary for Besag-Newell with c* =", cStar))
    cat("\n")
    cat("--------------------------------------------------------------------------------------------\n")
    cat(paste("No significant cluster found for c* =",cStar))
    cat("\n")
    cat("--------------------------------------------------------------------------------------------\n\n")
    sink()
    next 
  }
  plot(COMap$geometry, border = "grey", axes= TRUE,
       col=color.clusters(temp))
  title(main = paste("BN with c* =",cStar))
  # Write these results to a file for later reference
  
  sink(file="BNResults.txt", append = TRUE)
  cat("--------------------------------------------------------------------------------------------\n")
  cat(paste("Summary for Besag-Newell with c* =", cStar))
  cat("\n")
  cat("--------------------------------------------------------------------------------------------\n")
  print(summary(temp))
  cat("--------------------------------------------------------------------------------------------\n\n")
  sink()
}

# Spatial Scan statistic
# Expected rate.
e = sum(cases)/sum(pop) * pop
# apply circular scan method
scanRes = scan.test(coords = coordsXY,
                 cases = cases,
                 pop = pop,
                 ex = e,
                 nsim = 999,
                 alpha  = alpha)

# Grab the number of detected clusters
numCluster = length(clusters(scanRes))

# We want to color each cluster separately
mycol = grDevices::hcl.colors(numCluster)

# Pretty plot
plot(COMap$geometry, border = "grey60", axes=TRUE,
     col=color.clusters(scanRes, col=mycol))
title(main="Clusters found using circular spatial scan method")
# Write out the cluster information to a filef
sink(file="ScanResults.txt")
cat("--------------------------------------------------------------------------------------------\n")
cat("Summary for spatial scan statistic")
cat("\n")
cat("--------------------------------------------------------------------------------------------\n")
print(summary(scanRes))
cat("--------------------------------------------------------------------------------------------\n\n")
sink()

# Now we will use Moran's I
# Choose to be binary. IE w_{ij} = 1 iff region i shares a border with j
# Other methods could be used, but this is mostly exploratory
w = nb2mat(comNb, style = "B")
# Convert to a list for use with some of the functions
lw = nb2listw(comNb, style = "B")

# base test w/ randomization p-value
(ir = moran.mc(cases, listw = lw, nsim = 999))
# base test w/ Monto Carlo p-value, simulating data under constant risk hypothesis
# some preliminaries
N = length(cases) # number of regions
y = cases # number of cases
n = pop #population sizes
r <- sum(y)/sum(n) # estimated risk
rni <- r * n # expected per region

# observed moran's statistic
nsim = 999
t0 = moran(y, listw = lw, n = N, S0 = Szero(lw))$I
# simulate data under CRH
tsim = numeric(nsim)
# calculate moran's i for poisson data simulated under crh
for (i in 1:nsim) {
  tsim[i] = moran(rpois(N, rni), listw = lw, n = N, S0 = Szero(lw))$I
}

# p-value for moran's i constant risk monte carlo test
pVal = (sum(tsim >= t0) + 1)/(nsim + 1)
print(pVal)
