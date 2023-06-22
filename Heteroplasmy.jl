using StatsBase
using Dates
using Distributions
using Pipe
using DataFrames
using CSV
cd("/home/dima/Desktop/R/")

##############
# Simulation inputs
##############


function expand_grid(;iters...) 
   ## Courtesy to sjmgarnier
   ## https://discourse.julialang.org/t/function-like-expand-grid-in-r/4350/19
    var_names = collect(keys(iters))
    var_itr = [1:length(x) for x in iters.data]
    var_ix = vcat([collect(x)' for x in Iterators.product(var_itr...)]...)
    out = DataFrame()
    for i = 1:length(var_names)
        out[:,var_names[i]] = collect(iters[i])[var_ix[:,i]]
    end
    return out
end

##### Set parameters for simulations here
  simulation_params = expand_grid(
    grande_duplication_time = 90, 
    petite_duplication_time = 360,
    copy_number = 20,
    pathogeneity_threshold = 0.5,
    simulation_time = 5000,
    population_limit = 200,
    simulation_n = 1:10,
    suppressivity = 1:0.2:4,
    starting_rhom_proportion = 0.05:0.05:0.95,
    mutation_rate = 0
  )
  
#Initialize report DataFrame
report_df = DataFrame(
  sim_repeat = Int16[],
  # simulation params to report
  grande_duplication_time = Int16[],
  petite_duplication_time = Int16[],
  population_limit = Int16[],
  pathogeneity_threshold = Float16[],
  copy_number = Int16[],
  suppressivity = Float16[],
  starting_rhom_proportion = Float16[],
  mutation_rate = Float16[],
  # simulation results to report
  time = Int16[],
  rho_p_ratio = Float16[],
  fixed_rhom = Int16[],
  fixed_rhop = Int16[],
  Average_age = Float16[],  
)

p = Binomial(20)
binomial_noise  = [i-10 for i in rand(p,100)]

##############
#FUNCTIONS
##############



#Resolve heteroplasmy & mutagenesis
function resolve_heteroplasmy(rhop,rhom, simrun)
  if rhop == 0 
    return (0, rhom, 0, rhom) 
  end
  
  if rhom == 0 && simrun.mutation_rate == 0
    return (rhop,0, rhop, 0 )
  end
    
  if (1- simrun.mutation_rate)^rhop < rand()
    rhop = rhop - 1
    rhom = rhom + 1
  end
  
  mtDNA_mother = 
   @pipe vcat(fill(true, rhop*20),fill(false, Int(floor(rhom*simrun.suppressivity*20)))) |>
   @pipe  sample(_, 2*(rhop+rhom))
    
   g_1 = sum(mtDNA_mother[1:rhop+rhom])
   g_2 = sum(mtDNA_mother[rhop+rhom+1:end])

  return (g_1, rhop+rhom - g_1, g_2, rhop+rhom - g_2)
end
  
# Test phenotype
function test_phenotype(rhop, simrun)
  if rhop > simrun.copy_number * simrun.pathogeneity_threshold
    phenotype = true
  # TRUE means Grande
  else 
    phenotype = false
  end
  
  if phenotype == true
    tau12 = simrun.grande_duplication_time + sample(binomial_noise)
    # sample provide randomization that prevents culture synchronyzation (founder's effects)
  else 
    tau12 = simrun.petite_duplication_time + sample(binomial_noise)
    # sample provide randomization that prevents culture synchronyzation (founder's effects)
  end

  return (phenotype, tau12)
end

#Restart simulation
mutable struct yeast_cell
  rho_p::Int16
  rho_m::Int16
  phenotype::Bool
  tau12::Int16
  cellAge::Int16
end

function restart_simulation(simrun)
  population = []
  for i in 1:simrun.population_limit
    rhom = round(simrun.copy_number*simrun.starting_rhom_proportion)
    rhop = simrun.copy_number-rhom
    phenotype = test_phenotype(rhop, simrun)
    if phenotype == true
      max_age = simrun.grande_duplication_time
    else max_age = simrun.petite_duplication_time
    end
  
    new_cell = yeast_cell(
      rhop,
      rhom,     
      phenotype[1],
      phenotype[2],
      rand(1:max_age)
    )
    population = push!(population,new_cell)

  end
  println("new round...")
  return(population)
end


##############
# Itearation step!
##############
function iteration_step!(population, simrun)
  for cell in population
    cell.cellAge += 1
    if cell.cellAge >= cell.tau12
      new_proportions = resolve_heteroplasmy(cell.rho_p, cell.rho_m, simrun)
     
      # rewrite mother cell
      mother_phenotype = test_phenotype(new_proportions[1], simrun)
     
      cell.cellAge = 0
      cell.rho_p = new_proportions[1]
      cell.rho_m = new_proportions[2]
      cell.phenotype = mother_phenotype[1]
      cell.tau12 = mother_phenotype[2]


      #add daughter cell
      bud_phenotype = test_phenotype(new_proportions[3], simrun)

      new_cell = yeast_cell(
        new_proportions[3],
        new_proportions[4],
        bud_phenotype[1],
        bud_phenotype[2],
        0
      )
 
      push!(population, new_cell)    
    end
   
  end
      
  if length(population) > simrun.population_limit 
    population = population[sample(axes(population,1), simrun.population_limit, replace = false, ordered= true)]
  end
  return(population)  

end

function append_report!(population,report_df,simrun,simtime)
  new_row = DataFrame(
  sim_repeat = simrun.simulation_n,
  # simulation params to report
  population_limit = simrun.population_limit,
  grande_duplication_time = simrun.grande_duplication_time,
  petite_duplication_time = simrun.petite_duplication_time,
  pathogeneity_threshold = simrun.pathogeneity_threshold,
  copy_number = simrun.copy_number,
  suppressivity = simrun.suppressivity,
  starting_rhom_proportion = simrun.starting_rhom_proportion,
  mutation_rate = simrun.mutation_rate,
  # simulation results to report
  time = simtime,
  rho_p_ratio = mean([i.rho_p for i in population])/simrun.copy_number,
  fixed_rhom = count(i-> i == 0, [i.rho_p for i in population]),
  fixed_rhop = count(i-> i == 0, [i.rho_m for i in population]),
  Average_age = mean([i.cellAge for i in population]),
)

append!(report_df, new_row)

end
  
 
##############
# MAIN LOOP
##############
simrun = simulation_params[1, :]

counter = 0

for simrun in eachrow(simulation_params)
  counter +=1 
  population = restart_simulation(simrun)
  print(counter)
  print(" of ")
  println(size(simulation_params)[1])
  for simtime in 0:simrun.simulation_time
    population = iteration_step!(population, simrun)
    if simtime % 30 == 0 
      append_report!(population, report_df, simrun, simtime)
    end
  end
end

report_time = string(now())
mkdir(report_time)
cd(report_time)
CSV.write("simulation_params.csv", describe(simulation_params, :min, :max))
CSV.write("results.csv", report_df)




