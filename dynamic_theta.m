% load gpml
% add gpml path
gpml_path = "/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/gpml-matlab-v3.6-2015-07-07";
addpath(gpml_path);
startup;
rng("default");

t = 10;
likfunc = {@likGauss}; 
sn = 0.01;
hypSEiso.lik = log(sn);
x = (1:t)';
covfunc = {@covSEiso};
ell = 5;
sf = 1;
hypSEiso.cov = log([ell; sf]);
meanfunc = {@meanZero};
hypSEiso.mean = [];
K = feval(covfunc{:}, hypSEiso.cov, x) + 0.001*eye(t);
mu = feval(meanfunc{:}, hypSEiso.mean, x);
y = chol(K)'*normrnd(0,1,t,1) + mu + exp(hypSEiso.lik)*normrnd(0,1,t,1);

xs = linspace(1,t,100)';
[~,~,fs,s2] = gp(hypSEiso, @infExact, meanfunc, covfunc, likfunc, x, y, xs);

figure(1);
set(gca, 'FontSize', 24);
f = [fs+2*sqrt(s2); flipdim(fs-2*sqrt(s2),1)];
fill([xs; flipdim(xs,1)], f, [7 7 7]/8);
hold on; plot(xs, fs, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12);
grid on;
xlabel('Input (x)');
ylabel('Output (y)');
title("MAP");

% item responses
m = 20;
coefs = normrnd(0,1,[m,1]);
ps = 1./(1+exp(-coefs*x'));
responses = (rand([m,t])<ps);

theta = zeros(t,1);
iters = 100;
for iter=1:iters
    for i=1:t
        
    end
end

