### Version 0.0.4 

############################################
## Load libraries
############################################

library(tidyverse)
library(rbenchmark)
library(dqrng)

############################################
## Global params
############################################
simulation_params <- expand.grid(
  grande_duplication_time = 90,
  petite_duplication_time = 120,
  copy_number = 20, 
  pathogeneity_threshold = 0.5,
  simulation_time = 500,
  starting_CellNumber = 100,
  population_limit = 500,
  simulation_n = seq(1),
  suppressivity_select = seq(1,1.8,0.05),
  starting_rhom_proportion_select = 0.5)


############################################
#Initialize report data.frame
############################################
report <-tibble(sim_repeat = integer(),
                # simulation params to report
                    grande_duplication_time = integer(),
                    petite_duplication_time = integer(),
                    pathogeneity_threshold = numeric(),
                    copy_number = integer(),
                    suppressivity = numeric(),
                    starting_rhom_proportion = numeric(),
                # simulation results to report
                    time=integer(), 
                    n = numeric(),
                    rho_p_ratio = numeric(),
                    Average_age=numeric())


############################################
###FUNCTIONS
############################################

## test phenotype
test_phenotype <- function(rhop, rhom) {
  if (rhop > simrun$copy_number * simrun$pathogeneity_threshold) {
    phenotype = TRUE
    } # TRUE means Grande
  else {phenotype = FALSE}
  
  if (phenotype == TRUE) {tau12 = simrun$grande_duplication_time}
  else {tau12 = simrun$petite_duplication_time}
  
  return(list(phenotype,tau12))
}


## initiate data.frame
restart_simulation <- function() {
  CellNumber <- simrun$starting_CellNumber
  rho_p <- as.integer(simrun$copy_number * (1-simrun$starting_rhom_proportion))
  rho_m <- as.integer(simrun$copy_number * simrun$starting_rhom_proportion)
  phenotype = test_phenotype(rho_p, rho_m)[[1]]
  tau12 = as.integer(test_phenotype(rho_p, rho_m)[[2]])
  CellAge <- as.integer(runif(CellNumber, 0, tau12)) 

  population <- tibble(CellAge, 
                       rho_p,
                       rho_m, 
                       phenotype, 
                       tau12)

return(population)

}


## resolve heteroplasmy
resolve_heteroplasmy <- function(rhop, rhom) {
      # mind the multiplication factor 20 is set equal to copy number
  mtDNA_mother <- c(rep(TRUE, rhop*20),
                    rep(FALSE, simrun$suppressivity * 20 *rhom)) %>%
                dqsample(2*(rhop + rhom))
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

# increment age 
population_out$CellAge <- population_out$CellAge + 1L

for (cell in seq_along(1:nrow(population)))
  {

# resolving mtDNA segregation   
  if (population_out$CellAge[cell] >= population_out$tau12[cell]) {
      segregation <- resolve_heteroplasmy(population_out$rho_p[cell],
                                         population_out$rho_m[cell])
      
      population_out[cell,] <- tibble_row(CellAge = 0, 
                                rho_p = segregation[[1]][1], 
                                rho_m = segregation[[1]][2],
                                phenotype = test_phenotype(segregation[[1]][1],
                                                           segregation[[1]][2])[[1]],
                                tau12 = as.integer(test_phenotype(segregation[[1]][1],
                                                       segregation[[1]][2])[[2]]))
      
      
      new_cell <- tibble_row(CellAge =0,
                    rho_p = segregation[[2]][1], 
                    rho_m = segregation[[2]][2],
                    phenotype = test_phenotype(segregation[[2]][1],
                                               segregation[[2]][2])[[1]],
                    tau12 = as.integer(test_phenotype(segregation[[2]][1],
                                               segregation[[2]][2])[[2]])
                    ) 
      
      
      population_out <- bind_rows(population_out, new_cell)    }
  }

# test population size
  if (nrow(population_out) > simrun$population_limit) {
    population_out <- slice_sample(population_out, n = simrun$population_limit)
  }

  return(population_out)
  
}


############################################
##MAIN LOOP
############################################
set.seed(1)

for (rowname in seq(1:nrow(simulation_params[1:3,]))) {
  simrun <- simulation_params[rowname,]
  start_time <- Sys.time()
  
   
# population <- restart_simulation()
print(c(rowname, "of", 
        nrow(simulation_params)), quote = FALSE) 
population <- restart_simulation()

for (minute in seq(1:simrun$simulation_time))
{
  
population <- iteration_step(population)
head(population)

####### reporting per time
if (minute %% 30 == 0) {

  
iteration_result <- tibble(sim_repeat = simrun$simulation_n,
                            grande_duplication_time = simrun$grande_duplication_time,
                            petite_duplication_time = simrun$petite_duplication_time,
                            pathogeneity_threshold = simrun$pathogeneity_threshold,
                            copy_number = simrun$copy_number,
                            suppressivity = simrun$suppressivity,
                            starting_rhom_proportion = simrun$starting_rhom_proportion,
                            time=minute,
                            n = nrow(population),
                            rho_p_ratio = sum(population$rho_p)/(sum(population$rho_m) +
                                                                    sum(population$rho_p)),
                            Average_age=mean(population$CellAge),
                            fixed_rhom = sum(population$rho_p == 0),
                            fixed_rhop = sum(population$rho_m == 0))

report <- bind_rows(report, iteration_result)

}


}

time1 <- Sys.time()
print(time1 - start_time)

}

############################################
## GENERATE REPORT FILE
############################################
write.csv2(report, file = "simulation_report.csv")



