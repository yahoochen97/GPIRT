% illustrate changes in measurement w.r.t prior variance 
prior_stds = linspace(0.1,2,50)';

% p(y=1|x)=Phi(x), p(x)~N(0,sigma^2)
ys = [ones(9,1);zeros(1,1)];

% log p(x|ys) = log p(x) + sum log p(y_i|x) + C
% MAP estimation:

MAP = zeros(numel(prior_stds),1);
for i=1:numel(prior_stds)
   prior_std = prior_stds(i);
   xs = linspace(-3*prior_std,3*prior_std,100)';
   tmp = 1./(1+exp(-xs));
   lp = log(normpdf(xs,0,prior_std)) + mean(ys*log(tmp)'+(1-ys)*log(1-tmp)',1)';
   [~,idx] = max(lp);
   MAP(i) = xs(idx);
end

plot(prior_stds,MAP);