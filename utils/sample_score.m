function [next_score] = sample_score(xs, fstars, train_y, prior_score)
    [N,M]= size(train_y);
    next_score = zeros(M,1);
    for j=1:M
        xs_post = zeros(size(xs));
        for k=1:numel(xs)
            xs_post(k) = log(normpdf(xs(k),prior_score.mu, sqrt(prior_score.cov)))...
                +log_like_fn(fstars(:,k),train_y(:,j));
        end
        xs_post = cumsum(exp(xs_post));
        xs_post = (xs_post-xs_post(1))./(xs_post(end)-xs_post(1));
        next_score(j) = xs(end);
        r = unifrnd(0,1);
        for k=1:numel(xs)
           if xs_post(k)>r
               next_score(j) = xs(k); break;
           end
        end
    end
end

% function [ll] = log_like_fn(new_scores, cur_scores, cur_fs, ys)
% % summed log likelihood
%    model;
%    [N,~]= size(cur_fs);
%    ll = 0;
%    for i=1:N
%        [~,~,mu,~] = gp(hyp, @infExact, meanfunc,...
%                     covfunc, likfunc, cur_scores, cur_fs(i,:)', new_scores);
%        ps = sigmoid_fn(mu)';
%        ll = ll + sum(ys(i,:).*log(ps)+log(1-ps).*(1-ys(i,:)));
%    end
%    ll = ll / N;
%    
% end

