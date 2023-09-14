waddpath("~/Documents/Washu/CSE515T/Code/Gaussian Process/gpml-matlab-v3.6-2015-07-07");
startup;
data = readmatrix("./data/CIRI_theta.csv");
file_name = './results/gpirt_Supreme_Court_dynamic.csv';
opts = detectImportOptions(file_name);
data = readtable(file_name);
data = data(strcmp(data.type,"Martin-Quinn"),:);
data_mq = zeros(numel(unique(data.justice_id)),numel(unique(data.year)));
YEARS = unique(data.year);
JUSTICES = unique(data.justice_id);
for i=1:numel(unique(data.justice_id))
   for j=1:numel(unique(data.year))
      if numel(data(strcmp(data.year, YEARS(j)) & strcmp(data.justice_id, JUSTICES(i)), 'score').score)
        data_mq(i,j) = data(strcmp(data.year, YEARS(j)) & strcmp(data.justice_id, JUSTICES(i)), 'score').score;
      end
   end
end
data_mq(data_mq==0) = NaN;
data = data_mq;
% data = data ./ std(data,0,"all","omitnan");
n = size(data, 1);
horizon = size(data, 2);
x = (1:horizon)';

% joint model for all countries
meanfunc = [];                    % empty: don't use a mean function
covfunc = {@covMaterniso,5};              % Squared Exponental covariance function
likfunc = @likGauss;              % Gaussian likelihood

hyp = struct('mean', [], 'cov', [0 0], 'lik', log(0.01));
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
    n = size(data,1);
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