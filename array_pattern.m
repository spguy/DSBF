function pattern = array_pattern(weights, theta)
    array_factor = zeros(size(theta));
    for i = 1:length(weights)
        array_factor = array_factor + weights(i) * exp(1i * (i-1) * theta);
    end
    pattern = abs(array_factor);
end