# Author        :   Muhammad Reza Fahlevi
# Affiliation   :   Departemen Ilmu Komputer,
#                   Fakultas Ilmu Komputer - Teknologi Informasi,
#                   Universitas Sumatera Utara, Indonesia
# References    :   García-Martínez C., Rodriguez F.J., Lozano M. (2018) Genetic Algorithms. In: Martí R., Pardalos P., Resende M. (eds) Handbook of Heuristics. Springer, Cham. https://doi.org/10.1007/978-3-319-07124-4_28


library(ggplot2) # import library for data visualization
library(dplyr)

phi <- function(x, y) {
    exp(-0.5 * ((x ** 2) + (y ** 2))) %>% return()
}

config_df <- function(x_vect, y_vect) {
    phi_vect <- phi(x_vect, y_vect)
    dframe <- tibble(x_vect, y_vect, phi_vect)
    dframe %>%
        arrange(phi_vect) %>%
        slice(2:10) %>%
        return()
}

tournament <- function(df_pop) {
    df_pop %>%
        slice_sample(n = 4) %>%
        slice_max(order_by = phi_vect, n = 2) %>%
        return()
}

pbxalpha <- function(domain_interval, p1, p2, alpha) {
    constant <- (p1["x_vect"] - p2["x_vect"]) %>% abs() %>% max()
    p1_lower_bound <- (p1["x_vect"] - (alpha * constant)) %>% max(domain_interval[1])
    p1_upper_bound <- (p1["x_vect"] + (alpha * constant)) %>% min(domain_interval[2])
    
    constant <- (p1["y_vect"] - p2["y_vect"]) %>% abs() %>% max()
    p2_lower_bound <- (p2["y_vect"] - (alpha * constant)) %>% max(domain_interval[1])
    p2_upper_bound <- (p2["y_vect"] + (alpha * constant)) %>% min(domain_interval[2])
    
    p1_offspring <- c(p1_lower_bound, p1_upper_bound)
    p2_offspring <- c(p2_lower_bound, p2_upper_bound)
    
    p1_offspring %>%
        c(p2_offspring) %>%
        matrix(nrow = 2, ncol = 2, byrow = TRUE) %>%
        return()
}

gamma_ <- function(t, tmax, beta) (1 - (t / tmax)) ** beta

delta <- function(domain_interval, genotype, computed_gamma) {
    get_probability <- runif(1)
    probability_criteria <- 0.5
    normvar <- runif(1)
    if (get_probability == probability_criteria) {
        delta_value <- (domain_interval[2] - genotype) * ((1 - normvar) ** computed_gamma)
    } else {
        delta_value <- (domain_interval[1] - genotype) * ((1 - normvar) ** computed_gamma)
    }
    delta_value %>% return()
}

steady_state_genetic_algorithm <- function(population_size, mutation_probability, max_generation) {
    fitness <- c()
    x_vector <- runif(n = population_size) %>% c()
    y_vector <- runif(n = population_size) %>% c()
    
    for (generation in seq(1, max_generation)) {
        df_population <- config_df(x_vector, y_vector)
        df_population %>% print()
        fitness[generation] <- df_population[population_size - 1, 3] %>% max() # record the best
        x_vector <- c() # start a new population
        y_vector <- c()
        x_vector <- x_vector %>% append(df_population[population_size - 1, 1] %>% max()) # add the fittest to the next generation
        y_vector <- y_vector %>% append(df_population[population_size - 1, 2] %>% max())
        
        vect_size <- length(x_vector)
        while (vect_size != population_size) {
            winner <- tournament(df_population)
            p1 <- winner[1,]
            p2 <- winner[2,]
            offsprings <- pbxalpha(problem_interval, p1, p2, alpha = 0.5)
            current_mutation_rate <- runif(1)
            if (current_mutation_rate < mutation_probability) {
                p11 <- offsprings[1, 1] + delta(problem_interval, 
                                               offsprings[1, 1], gamma_(generation, max_generation, beta = 9))
                p12 <- offsprings[1, 2] + delta(problem_interval,
                                                offsprings[1, 2], gamma_(generation, max_generation, beta = 9))
                p21 <- offsprings[2, 1] + delta(problem_interval,
                                                offsprings[2, 1], gamma_(generation, max_generation, beta = 9))
                p22 <- offsprings[2, 2] + delta(problem_interval,
                                                offsprings[2, 2], gamma_(generation, max_generation, beta = 9))
                offsprings <- c(p11, p12, p21, p22) %>%
                    matrix(nrow = 2, ncol = 2, byrow = TRUE)
            }
            first_value <- offsprings[1, 1] %>% phi(offsprings[1, 2])
            second_value <- offsprings[2, 1] %>% phi(offsprings[2, 2])
            survivor <- c()
            if (first_value > second_value) {
                survivor <- offsprings[1, 1] %>% c(offsprings[1, 2])
            } else {
                survivor <- offsprings[2, 1] %>% c(offsprings[2, 2])
            }
            x_vector <- x_vector %>% append(survivor[1])
            y_vector <- y_vector %>% append(survivor[2])
            vect_size <- vect_size + 1
        }
    }
    fitness %>% return()
}

# define the parameter
problem_interval <- c(-5, 5)
NP <- 10
Pm <- 0.3
MAX_GENERATION <- 10

steady_state_genetic_algorithm(10, Pm, MAX_GENERATION)

result <- c(0.8503840, 0.9855303, 0.9983927, 0.9983927, 0.9996629,
            0.9997525, 0.9999375, 0.9999664, 0.9999841, 0.9999884)
idx <- seq(1, MAX_GENERATION)

df_result <- data.frame("generation" = idx, "fitness_value" = result)

plt <- ggplot(df_result, mapping = aes(generation, fitness_value)) + geom_line(color = "blue") + geom_point(color = "darkblue")
plt
