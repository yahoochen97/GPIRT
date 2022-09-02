addpath("~/Documents/Washu/CSE515T/Code/Gaussian Process/gpml-matlab-v3.6-2015-07-07");
startup;
data = readmatrix("./data/CIRI_theta.csv");
% data = data ./ std(data,0,"all","omitnan");
n = size(data, 1);
horizon = size(data, 2);
x = (1:horizon)';

% joint model for all countries
meanfunc = [];                    % empty: don't use a mean function
covfunc = @covSEiso;              % Squared Exponental covariance function
likfunc = @likGauss;              % Gaussian likelihood

hyp = struct('mean', [], 'cov', [0 0], 'lik', log(0.1));
prior.cov = {[], ...
             @priorDelta};
prior.mean = {};
prior.lik = {@priorDelta};

p.method = 'LBFGS';
p.length = 100;

inference_method = {@infPrior, @infExact, prior};

hyp = minimize_v2(hyp, @gp_sum, p, inference_method, meanfunc, ...
                    covfunc, likfunc, x, data);

disp(exp(hyp.cov));
disp(exp(hyp.lik));

function [nlZ, dnlZ] = gp_sum(hyp, inf, mean, cov, lik, x,data) 
    n = size(x, 1);
    nlZ = 0;
    dnlZ = unwrap(hyp);
    for i=1:n
        y = data(i,:)';
        mask = isnan(y);
       [this_nlZ, this_dnlZ ] = gp(hyp, inf, mean, cov, lik, x(~mask),y(~mask));
       nlZ = nlZ + this_nlZ;
       dnlZ = dnlZ + unwrap(this_dnlZ);
    end
    dnlZ = rewrap(hyp,dnlZ);
    dnlZ.lik = 0;
end