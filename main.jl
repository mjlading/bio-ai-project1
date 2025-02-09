include("simple_ga.jl")

using Plots

const GENOME_LENGTH = 5
POPULATION_SIZE = 100
MUTATION_RATE = 1 / (GENOME_LENGTH * 3)
CROSSOVER_RATE = 0.8 #TODO
N_GENERATIONS = 100

# println("Generating random population")
population = generate_random_population(POPULATION_SIZE)
# println.(population.individuals)

fitness_history = []
entropy_history = []

# Main evolution cycle
for generation in 1:N_GENERATIONS

    # println("\nGeneration $generation")

    parents = select_parents_roulette(population)
    offspring = crossover_single_point(parents)
    post_mutation_offspring = mutate(offspring, MUTATION_RATE)

    average_offspring_fitness = mean(i -> i.fitness, post_mutation_offspring)
    # println("average offspring fitness: ", average_offspring_fitness)
    push!(fitness_history, average_offspring_fitness)
    push!(entropy_history, shannon_entropy(post_mutation_offspring))

    replace_generation!(population, post_mutation_offspring)
end

println.(population.individuals)

best_individual = get_best_individual(population)
@show best_individual

@show length(fitness_history)
plot(1:N_GENERATIONS, fitness_history, xlabel="Generation", ylabel="Avg. fitness", title="Fitness over generations")

entropy_scaled = entropy_history .* (maximum(fitness_history) / maximum(entropy_history))
plot!(1:N_GENERATIONS, entropy_scaled, label="Entropy (scaled)", linestyle=:dash, color=:red)

hline!([961], label="Max Fitness")
