% Gaussian process regression with hyperparameter sampling demo

% load gpml
% add gpml path
gpml_path = "/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/gpml-matlab-v3.6-2015-07-07";
addpath(gpml_path);
startup;
rng("default");

% prepare the plotting environment
clear; close all;
write_fig = 0;

% Set the sample size.
n = 20;

% Likelihood function for the sample and priors for hyperparameters
likfunc = {@likGauss}; 
sn = 0.1;
hypSEiso.lik = log(sn);

% Generate n inputs for a data generating process
x = unifrnd(-2,2,n,1);

% Covariance function and priors for hyperparameters
covfunc = {@covSEiso};
ell = 1/2;
sf = 1;
% The isotropic squared exponential function has two hyperparamters.
hypSEiso.cov = log([ell; sf]);

% Mean function and priors for hyperparameters
meanfunc = {@meanSum, {@meanLinear, @meanConst}};
% The preceding mean function is the sum of a linear combination of the
% covariates plus a constant. If there are D dimensions in the data, this
% linear+constant mean function has D+1 hyperparameters.
hypSEiso.mean = [ones(size(x,2), 1), 0];
% The preceding hyperparameters set ones(size(x,2), 1) as the prior for
% each covariate in the linear portion of the mean function, and the zero
% sets the prior for the constant.

% Generate n outputs of a data generating process
K = feval(covfunc{:}, hypSEiso.cov, x) + 0.001*eye(n);
mu = feval(meanfunc{:}, hypSEiso.mean, x);
y = chol(K)'*normrnd(0,1,n,1) + mu + exp(hypSEiso.lik)*normrnd(0,1,n,1);

% Make a grid in the x dimension on which to predict values
z = linspace(-2, 2, 101)';

% Define hyperpriors
prior.cov = {{@priorTransform,@exp,@exp,@log,{@priorGamma,1,1}},... % Gamma prior on ell
               {@priorTransform,@exp,@exp,@log,{@priorGamma,1,2}}}; % Gamma prior on sf
prior.mean = { {@priorGauss ,1,1}, ... % Gaussian prior on mean slope
                {@priorGauss ,0,1}}; % Gaussian prior on mean constant
prior.lik = {{@priorTransform,@exp,@exp,@log,{@priorGamma,0.01,10}}}; % Gamma prior on noise

inffunc = {@infPrior, @infExact, prior};
p.method = 'LBFGS';
p.length = 100;

% Find the maximum a posteriori estimates for the model hyperparameters.
hypSEiso = minimize_v2(hypSEiso, @gp, p, inffunc, meanfunc, covfunc, likfunc, x, y);
fprintf("MAP ell: %.3f\n", exp(hypSEiso.cov(1)));
fprintf("MAP sf: %.3f\n", exp(hypSEiso.cov(2)));
fprintf("MAP slope: %.3f\n", (hypSEiso.mean(1)));
fprintf("MAP constant: %.3f\n", (hypSEiso.mean(2)));
fprintf("MAP noise: %.3f\n", exp(hypSEiso.lik));

% Report negative log marginal likelihood of the data with the optimized
% model hyperparameters.
nlml2 = gp(hypSEiso, inffunc, meanfunc, covfunc, likfunc, x, y);
% Make predictions
[m, s2] = gp(hypSEiso, inffunc, meanfunc, covfunc, likfunc, x, y, z);

% Plot predicted values and credible intervals
figure(1);
set(gca, 'FontSize', 24);
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
fill([z; flipdim(z,1)], f, [7 7 7]/8);
hold on; plot(z, m, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12);
grid on;
xlabel('Input (x)');
ylabel('Output (y)');
title("MAP");
if write_fig, print -depsc f3.eps; end;

% full hyperparameter sampling
% sampler parameters
num_chains  = 3;
num_samples = 1000;
burn_in     = 500;
jitter      = 1;

% setup sampler
% turn hypSEiso from struct to vector
% define objective function and its gradient
hyp = unwrap(hypSEiso);
f = @(unwrapped_hyp) nll(unwrapped_hyp, hypSEiso, inffunc, meanfunc, covfunc, x, y, likfunc);  
  
% create and tune sampler with jitter
hmc = hmcSampler(f, hyp + randn(size(hyp))*jitter);

tic;
[hmc, tune_info] = ...
    tuneSampler(hmc, ...
                'verbositylevel', 2, ...
                'numprint', 10, ...
                'numstepsizetuningiterations', 100, ...
                'numstepslimit', 500);
toc;

% use default seed for hmc sampler
rng('default');
tic;
for i=1:num_chains
[chains{i}, endpoint, acceptance_ratio] = ...
  drawSamples(hmc, ...
              'start', hyp + jitter * randn(size(hyp)), ...
              'burnin', burn_in, ...
              'numsamples', num_samples, ...
              'verbositylevel', 1, ...
              'numprint', 10);
end
toc;

% show diagnostics among chains
disp(diagnostics(hmc, chains));

% iterate all posterior samples
clear mus;
clear s2s;
for i=1:size(chains{1},1)
    hyp = rewrap(hypSEiso, chains{1}(i,:));

     [~,~,mu,s2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x, y, z);
    
    mus{i} = mu;
    s2s{i} = s2;
end

% plot posterior mean under different hyperparameter samples
figure(4);
for i=1:20:1000, plot(z, mus{i}); hold on; end

% compute gaussian mixture
gmm_mean = mean(cell2mat(mus),2);
gmm_s2 = mean(cell2mat(s2s),2);
gmm_var = gmm_s2 + mean(cell2mat(mus).^2,2) - gmm_mean.^2;

figure(2);
set(gca, 'FontSize', 24);
f = [gmm_mean+1.96*sqrt(gmm_var); flip(gmm_mean-1.96*sqrt(gmm_var),1)];
fill([z; flipdim(z,1)], f, [7 7 7]/8);
hold on; plot(z, gmm_mean, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12);
grid on;
xlabel('Input (x)');
ylabel('Output (y)');
title("Marginalized hyperparameter");

% plot hyperparameter posterior samples
figure(3);
c = exp(chains{1});
c(:, 4) = log(c(:,5));
c(:, 5) = log(c(:,5));
plotmatrix(c);

function [f,g] = nll(unwrapped_hyp, hypSEiso, inference_method, mean_function, covariance_function, x, y, likfunction)
% compute negative log likelihood and its gradient    
    hyp = rewrap(hypSEiso,unwrapped_hyp);
    [f, g] = gp(hyp, inference_method, mean_function,...
      covariance_function, likfunction, x, y);
  g = unwrap(g);
  f = -f;
  g = -g;
end