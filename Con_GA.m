clear
clc

gene_range = [0 100];
number_of_genes = 3;
population_size = 40; % must be multiple of 4
mutation_rate = 0.02;
generation_limit = 1000;

population = create_initial_population(population_size, gene_range, number_of_genes);
generation_count = 0;
while true
    fitness_values = calc_fitness(population);
    generation = [population fitness_values];
    generation = sortrows(generation, number_of_genes + 1, 'descend'); % sort by fitness values
    
    % check termination criteria
    if generation_count > generation_limit
        break;
    end
    
    % create next generation
    sorted_population = generation(:, 1:end-1);
    next_population = perform_selection_pairing_mating(sorted_population);
    muted_population = perform_mutation(next_population, gene_range, mutation_rate);
    population = muted_population;
    generation_count = generation_count + 1;
end

generation(1, 1:end-1)
generation(1, end)
%%
% www.wolframalpha.com
% Maximize[{2xz*exp(-x)-2y^3+y^2-3z^3,0<=x<=100 && 0<=y<=100 && 0<=z<=100}, {x,y,z}]
clear
clc
e = exp(1)
[1 1/3 sqrt(2/e)/3]
fitness_function([1 1/3 sqrt(2/e)/3])

%%
function population = create_initial_population(population_size, gene_range, number_of_genes)
    range = gene_range(2) - gene_range(1);
    lower_limit = gene_range(1);
    population = rand(population_size, number_of_genes) * range + lower_limit;
end

function fitness_values = calc_fitness(population)
    population_size = size(population, 1);
    fitness_values = zeros(population_size, 1);
    for i = 1:population_size
        fitness_values(i) = fitness_function(population(i, :));
    end
end

function fitness = fitness_function(individual)
    % maximize problem => we need fitness function, not cost
    x = individual(1); y = individual(2); z = individual(3);
    fitness = 2*x*z*exp(-x)-2*y^3+y^2-3*z^3;
end

function next_population = perform_selection_pairing_mating(sorted_population)
    population_size = size(sorted_population, 1);
    number_of_genes = size(sorted_population, 2);
    half_population_size = 0.5 * population_size;
    next_population = zeros(population_size, number_of_genes);
    
    % fittest half selection
    next_population(1:half_population_size, :) = sorted_population(1:half_population_size, :);
    
    % pairing and mating
    for i = 1:2:half_population_size
        parents = sorted_population(i:i+1,:);
        offsprings = perform_mating(parents);
        next_population(half_population_size+i+0,:) = offsprings(1,:);
        next_population(half_population_size+i+1,:) = offsprings(2,:);
    end
end

function offsprings = perform_mating(parents)
    offsprings = parents;
    alpha = rand(1);
    offsprings(1, 1) = alpha*parents(1, 1) + (1-alpha)*parents(2, 1);
    offsprings(2, 1) = alpha*parents(2, 1) + (1-alpha)*parents(1, 1);
end

function muted_population = perform_mutation(next_population, gene_range, mutation_rate)
    % find which genes should be muted
    muted_population = next_population;
    temp = rand(size(next_population));
    mutant_indices = find(temp < mutation_rate);
    % perform reset mutation
    range = gene_range(2) - gene_range(1);
    lower_limit = gene_range(1);
    number_of_mutants = length(mutant_indices);
    mutants = rand(number_of_mutants, 1) * range + lower_limit;
    muted_population(mutant_indices) = mutants;
end