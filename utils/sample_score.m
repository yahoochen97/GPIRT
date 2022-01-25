function [next_score, ll] = sample_score(cur_score, prior_score, cur_fs, ys)
    prior_sample = normrnd(prior_score.mu, prior_score.cov);
    % define log likelihood functions
    ll_score = @(score) log_like_fn(score, cur_score, cur_fs, ys);
    [next_score, ll] = elliptical_slice(cur_score, prior_sample, ll_score);
end

function [ll] = log_like_fn(score, cur_score, cur_fs, ys)
% summed log likelihood
    meanfunc = {@meanZero};
    covfunc = {@covSEiso};
    likfunc = {@likGauss};
    hyp = struct('mean', [], 'cov', [0, log(1)], 'lik', -10);

   [~,~,mu,s2] = gp(hyp, @infExact, meanfunc,...
                covfunc, likfunc, cur_score*ones(size(cur_fs)), cur_fs, score);
   fs = normrnd(mu, s2);
   ps = sigmoid_fn(fs);
   ll = sum(ys.*log(ps)+log(1-ps).*(1-ys));
end

