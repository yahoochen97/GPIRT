function [ps] = sigmoid_fn(fs)
% Sigmoid link function
    ps = 1./(1+exp(-fs));
end

