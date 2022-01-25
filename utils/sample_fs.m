function [next_fs, ll] = sample_fs(cur_fs, prior_fs, ys)
    prior_sample = mvnrnd(prior_fs.mu, prior_fs.cov);
    % define log likelihood functions
    ll_fs = @(fs) log_like_fn(fs, ys);
    [next_fs, ll] = elliptical_slice(cur_fs, prior_sample, ll_fs);
end

function [ll] = log_like_fn(fs, ys)
% summed log likelihood
   ps = sigmoid_fn(fs);
   ll = sum(ys.*log(ps)+log(1-ps).*(1-ys));
end
