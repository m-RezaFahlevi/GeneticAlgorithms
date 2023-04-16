library(dplyr)
library(ggplot2)
library(gganimate)
library(readr)

# rcga_df <- read_csv("data/obt_data1681577279.csv",
#                                col_types = cols(is_trans = col_logical()))
# rcga_df <- obt_data1681577275 <- read_csv("data/obt_data1681577275.csv", 
#                                           col_types = cols(is_trans = col_logical()))
# rcga_df <- read_csv("data/obt_data1681608301.csv", 
#                                col_types = cols(is_trans = col_logical()))
# rcga_df <- read_csv("data/obt_data1681609166.csv", 
#                                col_types = cols(is_trans = col_logical()))
rcga_df <- read_csv("data/obt_data1681611939.csv", 
                               col_types = cols(is_trans = col_logical()))

wrangling_data <- function(obt_data) {
    n_iter <- seq_len(dim(obt_data)[1])
    obt_data[,3] <- n_iter
    names(obt_data)[3] <- "n_iter"
    return(obt_data)
}
rcga_df <- wrangling_data(rcga_df)
summary(rcga_df)
View(rcga_df)

binom_data <- function(obt_data) {
    n_trans <- obt_data %>%
        count(is_trans)
    no_trans <- n_trans[1,2] / sum(n_trans[,2])
    trans <- n_trans[2,2] / sum(n_trans[,2])
    binom_df <- data.frame(
        "state" = c("trans", "no_trans"),
        "probability" = c(trans, no_trans)
    )
    return(binom_df)
}
rcga_binom <- binom_data(rcga_df)
summary(rcga_binom)
View(rcga_binom)

acc_data <- function(obt_data) {
    acc_val <- obt_data %>%
        filter(is_trans) %>%
        select(fitval) %>%
        as.data.frame()
    acc_val <- acc_val[,1]
    iter_num <- length(acc_val) %>% seq_len()
    acc_fitval <- data.frame(
        "iter_num" = iter_num,
        "fitval" = acc_val
    )
    names(acc_fitval)
    head(acc_fitval)
    tail(acc_fitval)
    summary(acc_fitval)
    return(acc_fitval)
}
rcga_acc <- acc_data(rcga_df)
summary(rcga_acc)
View(rcga_acc)

datname = "obt_data1681611939.csv"

ggplot(rcga_acc[1:200,], aes(iter_num, fitval)) +
    geom_point(color = "darkblue", size = 0.7) +
    geom_line(color = "blue") +
    ggtitle(
        label = "Real-Coded Genetic Algorithm",
        subtitle = datname
    ) + xlab("nth-evaluation") +
    ylab("fitness") +
    theme_bw()

rcga_reg <- rcga_acc[1:300,]


ggplot(rcga_reg, aes(iter_num, fitval)) +
        geom_point(color = "darkblue", size = 0.8) +
        ggtitle(
            label = "Real-Coded Genetic Algorithm",
            subtitle = datname
        ) +
        xlab("nth-evaluation") + ylab("fitness") +
        theme_bw()

ggplot(rcga_reg, aes(iter_num, fitval)) +
    geom_line(color = "blue") +
    ggtitle(
        label = "Real-Coded Genetic Algorithm",
        subtitle = datname
    ) +
    xlab("nth-evaluation") + ylab("fitness") +
    theme_bw()
    
ggplot(rcga_reg, aes(iter_num, fitval)) +
        geom_point(color = "darkblue", size = 0.8) +
        stat_smooth(formula = y ~ s(x, bs="cs"), method = "gam") +
        ggtitle(
            label = "Real-Coded Genetic Algorithm",
            subtitle = datname
        ) +
        xlab("nth-evaluation") + ylab("fitness") +
        theme_bw()
    
rcga_gif <- ggplot(rcga_reg, aes(iter_num, fitval)) +
        geom_point(color = "darkblue", size = 0.8) +
        geom_line(color = "blue") +
        ggtitle(
            label = "Real-Coded Genetic Algorithm",
            subtitle = datname
        ) +
        xlab("nth-evaluation") + ylab("fitness") +
        theme_bw() +
    transition_manual(iter_num, cumulative = TRUE)
rcga_gif
anim_save("visualization/obt_data1681611939.gif")
