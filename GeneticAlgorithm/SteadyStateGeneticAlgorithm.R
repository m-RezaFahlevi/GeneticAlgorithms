# Author        :   Muhammad Reza Fahlevi
# Affiliation   :   Departemen Ilmu Komputer,
#                   Fakultas Ilmu Komputer - Teknologi Informasi,
#                   Universitas Sumatera Utara, Indonesia
# References    :   García-Martínez C., Rodriguez F.J., Lozano M. (2018) Genetic Algorithms. In: Martí R., Pardalos P., Resende M. (eds) Handbook of Heuristics. Springer, Cham. https://doi.org/10.1007/978-3-319-07124-4_28


library(ggplot2) # import library for data visualization

phi <- function(x, y) exp(-0.5 * (x ** 2) + (y ** 2))

# some value of phi(x,y) for 0 <= x = y <= 0.5
phi(seq(0, 0.5, 0.05), seq(0, 0.5, 0.05))

xgen <- runif(7, min = -1, max = 1)
ygen <- runif(7, min = -1, max = 1)
popdf <- data.frame("xgenotype" = xgen, "ygenotype" = ygen, "fitness" = phi(xgen, ygen))
popdf

dd <- transform(popdf, fitness = factor(popdf$fitness))
popdfsorted <- dd[ do.call(order, dd["fitness"]), ]
popdfsorted

SIGMA <- sum(popdf$fitness)
popdfsorted["pvalue"] <- phi(popdfsorted$xgenotype, popdfsorted$ygenotype) / SIGMA
popdfsorted

# Create Data
data <- data.frame(
    group=paste(popdfsorted$pvalue * 100),
    value= popdfsorted$pvalue
)

# Basic piechart
ggplot(data, aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + ggtitle("Roulette-Wheel Selection")

# Create roulette-wheel function
roulette <- function(popdframe) {
    rouletteWheel <- runif(1, min = min(popdframe$pvalue), 
                           max = max(popdframe$pvalue))
    if (rouletteWheel > 0 && rouletteWheel <= popdframe$pvalue[1]) {
        parent <- c(popdframe$xgenotype[1], popdframe$ygenotype[1])
    } else if (rouletteWheel > popdframe$pvalue[1] && rouletteWheel <= popdframe$pvalue[2]) {
        parent <- c(popdframe$xgenotype[2], popdframe$ygenotype[2])
    } else if (rouletteWheel > popdframe$pvalue[2] && rouletteWheel <= popdframe$pvalue[3]) {
        parent <- c(popdframe$xgenotype[3], popdframe$ygenotype[3])
    } else if (rouletteWheel > popdframe$pvalue[3] && rouletteWheel <= popdframe$pvalue[4]) {
        parent <- c(popdframe$xgenotype[4], popdframe$ygenotype[4])
    } else if (rouletteWheel > popdframe$pvalue[4] && rouletteWheel <= popdframe$pvalue[5]) {
        parent <- c(popdframe$xgenotype[5], popdframe$ygenotype[5])
    } else if (rouletteWheel > popdframe$pvalue[5] && rouletteWheel <= popdframe$pvalue[6]) {
        parent <- c(popdframe$xgenotype[6], popdframe$ygenotype[6])
    } else if (rouletteWheel > popdframe$pvalue[6] && rouletteWheel <= popdframe$pvalue[7]) {
        parent <- c(popdframe$xgenotype[7], popdframe$ygenotype[7])
    }
    # print(rouletteWheel)
    return(parent)
}

# Select two chromosome
selection <- function(popdataframe) {
    stparents <- roulette(popdataframe)
    ndparents <- roulette(popdataframe)
    chromosomes <- matrix(c(stparents, ndparents),
                          nrow = 2, ncol = 2, byrow = TRUE)
    return(chromosomes)
}

# Crossover: PBX-\alpha method
pbxalpha <- function(selMat, constAlpha, stDomain, ndDomain) {
    stI <- abs(selMat[1, 1] - selMat[2, 1])
    ndI <- abs(selMat[1, 2] - selMat[2, 2])
    p11 <- max(stDomain[1], min(selMat[1, 1], selMat[2, 1]) - (constAlpha * stI))
    p12 <- min(stDomain[2], max(selMat[1, 1], selMat[2, 1]) + (constAlpha * stI))
    p21 <- max(ndDomain[1], min(selMat[1, 2], selMat[2, 2]) - (constAlpha * ndI))
    p22 <- min(ndDomain[2], max(selMat[1, 2], selMat[2, 2]) + (constAlpha * ndI))
    crossover <- matrix(c(p11, p12, p21, p22), nrow = 2, ncol = 2, byrow = TRUE)
    return(crossover)
}
offspring <- pbxalpha(selectedMatrix, runif(1, min = 0.01, max = 0.99), 
                      stDomain = c(0.01, 0.99), ndDomain = c(0.01, 0.99))
offspring

# Mutation function
computedGamma <- function(tCurrent, tMax, constBeta) (1 - (tCurrent / tMax)) ** constBeta
computedDelta <- function(theDomain, theVar, constGamma) {
    deltaVect <- c()
    for (i in seq(1, 2)) {
        drawProbability <- runif(1)
        stOutcome <- (theDomain[1] - theVar[i]) * (1 - runif(1)) ** constGamma
        ndOutcome <- (theDomain[2] - theVar[i]) * (1 - runif(1)) ** constGamma
        ifelse(drawProbability == (1 / 2),
               deltaVect <- append(deltaVect, stOutcome),
               deltaVect <- append(deltaVect, ndOutcome))
    }
    return(deltaVect)
}
mutation <- function(varDelta, varX) varDelta + varX

# is mutation occur? function
ismutate <- function(toOffspring, mutationProbability, maxgeneration) {
    mutationChance <- runif(1)
    if (mutationChance < MUTATIONRATE) {
        compGamma <- computedGamma(tCurrent = 1, tMax = maxgeneration, constBeta = 3)
        stCompDelta <- computedDelta(theDomain = c(0.01, 0.99), 
                                     theVar = offspring[1,], constGamma = compGamma)
        ndCompDelta <- computedDelta(theDomain = c(0.01, 0.99), 
                                     theVar = offspring[2,], constGamma = compGamma)
        stCompMutation <- mutation(stCompDelta, offspring[1,])
        ndCompMutation <- mutation(ndCompDelta, offspring[2,])
        compMutation <- matrix(c(stCompMutation, ndCompMutation),
                               nrow = 2, ncol = 2, byrow = TRUE)
        return(compMutation)
    } else {
        return(offspring)
    }
}


# Define the constant
NP <- 7
MAXGENERATION <- 50
MUTATIONRATE <- 0.3
STDOMAIN <- c(-0.99, 0.99)
NDDOMAIN <- c(-0.99, 0.99) 

# Repeat the process until stop criterea is found
generations[[1]] <- popdf
for (numberGeneration in seq(1, MAXGENERATION)) {
    xgenvect <- c()
    ygenvect <- c()
    while (length(xgenvect) < NP) {
        # Config the population
        workingPopulation <- generations[[numberGeneration]]
        dd <- transform(workingPopulation, 
                        fitness = factor(workingPopulation$fitness))
        workingPopulationSorted <- dd[ do.call(order, dd["fitness"]), ]
        
        SIGMA <- sum(workingPopulation$fitness)
        workingPopulationSorted["pvalue"] <- phi(workingPopulationSorted$xgenotype, workingPopulationSorted$ygenotype) / SIGMA
        
        # Selection
        selectIndividu <- selection(workingPopulationSorted)
        
        # Crossover
        getOffspring <- pbxalpha(selectIndividu,
                                 constAlpha = runif(1, min = 0.01, max = 0.99),
                                 stDomain = STDOMAIN, ndDomain = NDDOMAIN)
        # Mutation
        getOff <- ismutate(toOffspring = getOffspring,
                           mutationProbability = MUTATIONRATE,
                           maxgeneration = MAXGENERATION)
        
        # calculation the fitness by using the objective function
        stCalc <- phi(getOff[1, 1], getOff[1, 2])
        ndCalc <- phi(getOff[2, 1], getOff[2, 2])
        if (stCalc > ndCalc) {
            xgenvect <- append(xgenvect, getOff[1, 1])
            ygenvect <- append(ygenvect, getOff[1, 2])
        } else {
            xgenvect <- append(xgenvect, getOff[2, 1])
            ygenvect <- append(ygenvect, getOff[2, 2])
        }
    }
    getNewPopulation <- data.frame("xgenotype" = xgenvect,
                                   "ygenotype" = ygenvect,
                                   "fitness" = phi(xgenvect, ygenvect))
    generations[[numberGeneration + 1]] <- getNewPopulation
}

# print the last 10 generation
for (n in seq(MAXGENERATION - 10, MAXGENERATION)) {
    print(n)
    print(generations[[n]])
}
