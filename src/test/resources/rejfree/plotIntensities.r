list.of.packages <- c("ggplot2", "reshape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="http://cran.us.r-project.org")

require(ggplot2)
require(reshape)

data <- read.csv( comment.char = '#',
                  file = 
                    "@{samples.getAbsolutePath()}",
                  header = TRUE)


reshaped <- melt(data)

plot <- ggplot(
  subset(reshaped, (variable != "lp__" & variable != "accept_stat__" & variable != "stepsize__" & variable != "treedepth__" & variable != "n_leapfrog__" & variable != "n_divergent__")), 
  aes(x = reorder(variable, value, FUN=median), y = exp(value))) + 
  coord_flip() +
  geom_boxplot()

ggsave(filename = "@{getIntensitiesOutput().getAbsolutePath()}", width = 20, height = 20)
