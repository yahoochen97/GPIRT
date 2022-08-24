% add gpml path
gpml_path = "/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/gpml-matlab-v3.6-2015-07-07";
addpath(gpml_path);
addpath("utils");
startup;
rng("default");

x = [(1:9)';(11:30)'];
y = x + sin(2*x) + normrnd(0,1,size(x));
xs = (10)';

hyp =  struct('mean', [], 'cov', [log(2) 0], 'lik', log(0.01));

mu = covSEiso(hyp.cov,xs,x)*inv(covSEiso(hyp.cov,x,x)+exp(2*hyp.lik)*eye(numel(x)))*y;
s2 = covSEiso(hyp.cov,xs,xs)-covSEiso(hyp.cov,xs,x)*inv(covSEiso(hyp.cov,x,x)+exp(2*hyp.lik)*eye(numel(x)))*covSEiso(hyp.cov,x,xs);

disp(sqrt(diag(s2)));

% [~,~, mu, s2] = gp(hyp, @infExact, [],...
%             @covSEiso, @likGauss, x, y, xs);

sample = mvnrnd(mu,s2);
        
f = [mu+2*sqrt(diag(s2)); flipdim(mu-2*sqrt(diag(s2)),1)];
  fill([xs; flipdim(xs,1)], f, [7 7 7]/8);
  hold on; plot(xs, mu); plot(x, y, '+'); plot(xs,sample, '-g');

cur_mu = 0; 
cur_s2 = exp(2*hyp.lik);  
for i=1:numel(xs)
    cur_mu = 0;
end
 