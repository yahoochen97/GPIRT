function [ll] = log_like_fn(fs, ys)
    % summed log likelihood
   ps = sigmoid_fn(fs);
   ll = sum(ys.*log(ps)+log(1-ps).*(1-ys));
end

