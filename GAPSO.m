% 参数设置
num_elements = 10;      % 阵列元素数量
theta = linspace(-pi/2, pi/2, 181); % 角度范围
desired_pattern = cos(2*pi*0.5*sin(theta)); % 期望波束图

% 遗传算法参数
ga_population_size = 50;
ga_generation = 100;

% 粒子群算法参数
pso_swarm_size = 50;
pso_iterations = 100;

% 生成随机阵列权重
initial_weights = rand(1, num_elements);

% 遗传算法优化
ga_options = optimoptions('ga', 'PopulationSize', ga_population_size, 'MaxGenerations', ga_generation);
ga_fitness_function = @(weights) -beamforming_fitness(weights, theta, desired_pattern);
ga_best_weights = ga(ga_fitness_function, num_elements, [], [], [], [], zeros(1, num_elements), ones(1, num_elements), [], ga_options);

% 粒子群算法优化
pso_options = optimoptions('particleswarm', 'SwarmSize', pso_swarm_size, 'MaxIterations', pso_iterations);
pso_fitness_function = @(weights) -beamforming_fitness(weights, theta, desired_pattern);
pso_best_weights = particleswarm(pso_fitness_function, num_elements, zeros(1, num_elements), ones(1, num_elements), pso_options);

% 计算波束图
ga_beam_pattern = array_pattern(ga_best_weights, theta);
pso_beam_pattern = array_pattern(pso_best_weights, theta);

% 绘制结果
figure;
subplot(2, 1, 1);
plot(theta, 20*log10(abs(ga_beam_pattern)), 'r', 'LineWidth', 2);
hold on;
plot(theta, 20*log10(abs(desired_pattern)), 'k--', 'LineWidth', 2);
title('GA Beamforming');
legend('Actual Pattern', 'Desired Pattern');
xlabel('Angle (rad)');
ylabel('Magnitude (dB)');

subplot(2, 1, 2);
plot(theta, 20*log10(abs(pso_beam_pattern)), 'b', 'LineWidth', 2);
hold on;
plot(theta, 20*log10(abs(desired_pattern)), 'k--', 'LineWidth', 2);
title('PSO Beamforming');
legend('Actual Pattern', 'Desired Pattern');
xlabel('Angle (rad)');
ylabel('Magnitude (dB)');

function fitness = beamforming_fitness(weights, theta, desired_pattern)
    array_pattern = array_pattern_func(weights, theta);
    fitness = -sum((abs(array_pattern) - abs(desired_pattern)).^2);
end

function pattern = array_pattern_func(weights, theta)
    array_factor = zeros(size(theta));
    for i = 1:length(weights)
        array_factor = array_factor + weights(i) * exp(1i * (i-1) * theta);
    end
    pattern = abs(array_factor);
end
