function fitness = beamforming_fitness(weights, theta, desired_pattern)
    array_pattern = array_pattern(weights, theta);
    fitness = -sum((abs(array_pattern) - abs(desired_pattern)).^2);
end