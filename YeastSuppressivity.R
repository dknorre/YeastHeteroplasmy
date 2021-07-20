############################################
## Load libraries
############################################

library(tidyverse)
library(rbenchmark)

############################################
## Global params
############################################
grande_duplication_time <- 90
petite_duplication_time <- c(120)
copy_number <- 20 
pathogeneity_threshold <- 10
simulation_time <- 2000
starting_CellNumber <- 100
population_limit <- 100
simulation_n <- 10
suppressivity_select <- seq(1,1.8,0.05)
starting_rhom_proportion_select <- seq(0,1,0.05)


############################################
#Initialize report data.frame
############################################

set.seed(1)
report <-tibble(sim_repeat = integer(),
                    copy_number = integer(),
                    Time=integer(), 
                    n = numeric(),
                    rho_p_ratio = numeric(),
                    Average_age=numeric())


############################################
###FUNCTIONS
############################################

## test phenotype
test_phenotype <- function(rhop, rhom) {
  if (rhop > pathogeneity_threshold) {phenotype = TRUE} # TRUE means Grande
  else {phenotype = FALSE}
  return(phenotype)
}


## initiate data.frame
restart_simulation <- function() {
CellNumber <- starting_CellNumber
CellAge <- as.integer(runif(CellNumber, 0, petite_duplication_time)) 
# upgrade this! check maximum age as a function of pathogeneity threshold
rho_p <- copy_number * (1-starting_rhom_proportion)
rho_m <- copy_number * starting_rhom_proportion
phenotype = test_phenotype(rho_p, rho_m)
population <- tibble(CellAge, 
                         rho_p,
                         rho_m, 
                        phenotype)

return(population)

}


## resolve heteroplasmy
resolve_heteroplasmy <- function(rhop, rhom) {
  mtDNA_mother <- c(rep(TRUE, 100*rhop),rep(FALSE, suppressivity * 100 *rhom))
  mtDNA_mother <- sample(mtDNA_mother, 2*rhop + 2*rhom)
  g_1 <- mtDNA_mother[1:(length(mtDNA_mother)/2)]
  g_2 <- tail(mtDNA_mother, length(mtDNA_mother)/2)
  genotype_1 <- c(sum(g_1), sum(!g_1))
  genotype_2 <- c(sum(g_2), sum(!g_2))
  return(list(genotype_1, genotype_2))
}

## ITERATION STEP
#calculation
iteration_step <- function(population) {
 
population_out <- population
cell = 1

# increment age 
population_out$CellAge <- population_out$CellAge + 1

for (cell in seq_along(1:nrow(population)))
  {

    
# budding    
  if  (population_out$phenotype[cell] == FALSE ) {
    dupl_time = petite_duplication_time
  }
  
  if  (population_out$phenotype[cell] == TRUE ) {
    dupl_time = grande_duplication_time
  }

  if (population_out$CellAge[cell] >= dupl_time) {
      segregation <- resolve_heteroplasmy(population_out$rho_p[cell],
                                         population_out$rho_m[cell])
      
      population_out[cell,] <- tibble_row(CellAge = 0, 
                                rho_p = segregation[[1]][1], 
                                rho_m = segregation[[1]][2],
                                phenotype = test_phenotype(segregation[[1]][1],
                                                           segregation[[1]][2]))
      new_cell <- tibble_row(CellAge =0,
                    rho_p = segregation[[2]][1], 
                    rho_m = segregation[[2]][2],
                    phenotype = test_phenotype(segregation[[2]][1],
                                               segregation[[2]][2])
                    )
      population_out <- bind_rows(population_out, new_cell)    }
  }

# test population size
  if (nrow(population_out) > population_limit) {
    population_out <- sample_n(population_out, population_limit)
  }

  return(population_out)
  
}

############################################
##MAIN LOOP
############################################

for (starting_rhom_proportion in starting_rhom_proportion_select) {
  pathogeneity_threshold <- copy_number
for (suppressivity in suppressivity_select) { 
for (sim_repeat in seq_along(1:simulation_n)) {
  
population <- restart_simulation()
print(c(sim_repeat, "--", copy_number, "-- ", suppressivity)) 

for (minute in seq(1:simulation_time))
{
population <- iteration_step(population)

#reporting per time
if (minute %% 60 == 0) {
iteration_at_start <- tibble(sim_repeat = sim_repeat,
                             iterator = starting_rhom_proportion,
                             iterator2 = suppressivity,
                             Time=minute, 
                             n = nrow(population),
                             rho_p_ratio = sum(population$rho_p)/(sum(population$rho_m) + 
                                                                    sum(population$rho_p)),
                             Average_age=mean(population$CellAge),
                             fixed_rhom = sum(population$rho_p == 0),
                             fixed_rhop = sum(population$rho_m == 0))
report <- bind_rows(report, iteration_at_start)  
}
    
}
}
}
}  


############

write.csv2(report, file = "simulation_report.csv")

#########IMAGING#



ggplot(report, aes(x = Time/60, y = rho_p_ratio, group = sim_repeat)) +
  geom_ribbon(aes(ymin = 1-fixed_rhop/100), ymax = 1, alpha = 0.05, fill = "blue") +
  geom_ribbon(aes(ymin = 0, ymax = fixed_rhom/100), alpha = 0.05, fill = "red") +
  geom_line(alpha = 0.5) +
  ylim(0,1) +
  facet_grid(iterator ~ iterator2) +
  xlab("Time, h") + 
  ggtitle(paste(
              c(" grande_duplication_time ", 
                grande_duplication_time,
                "\n pathogeneity_threshold",
                pathogeneity_threshold, 
              "\n copy_number_select",
                copy_number,
              "\n suppressivity_select",
                suppressivity_select),
                collapse = " - "
                )) +
  theme_bw()

############### draw heatmap for specific time
slice_fhp <- report %>% filter(Time == 1800) %>% 
  group_by(iterator, iterator2) %>% summarise(fixed_rhom = mean(fixed_rhom))

slice_fhp %>% ggplot(aes(x = iterator,y = iterator2, fill = fixed_rhom))+
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "red") + 
  scale_x_continuous() +
  scale_y_continuous(breaks = suppressivity_select) +
  theme_minimal()





