% two way fixed effects using Gaussian process regression demo

% load gpml
% add gpml path
gpml_path = "/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/gpml-matlab-v3.6-2015-07-07";
addpath(gpml_path);
startup;
rng("default");

% prepare the plotting environment
clear; close all;
write_fig = 0;

% Set the panel size.
n = 10;
T = 30;
noise = 0.1;

% generate unit fixed effects
unit_mean = 0;
unit_std = 0.5;
unit_effects = normrnd(unit_mean, unit_std, [n,1]);

% generate time random effects
time_ls = 10;
time_os = 1;
K = covSEiso([log(time_ls), log(time_os)], (1:T)') + 0.00001*eye(T);
mu = zeros(T,1);
time_effects = chol(K)'*normrnd(0,1,T,1) + mu;

% plot time random effects
figure(1);
set(gca, 'FontSize', 24);
plot(1:T, time_effects, 'LineWidth', 2);
grid on;
xlabel('time');
ylabel('effect');
title("time random effect");
if write_fig 
    filename = "../results/timeRE.pdf";
    print(fig, filename, '-dpdf','-r300');
    close;
end

x = zeros(n*T, 2); % (unit, time)
y = zeros(n*T, 1); % unit effect + time effect + noise 
for i=1:n
    for j=1:T
        x((i-1)*T+j, 1) = i;
        x((i-1)*T+j, 2) = j;
        y((i-1)*T+j) = time_effects(j) + unit_effects(i) + normrnd(0, noise);
    end
end

% plot data
figure(2);
set(gca, 'FontSize', 24);
for i=1:n
    plot(1:T, y(((i-1)*T+1):(i*T)), 'LineWidth', 2);
    hold on;
end
grid on;
xlabel('time');
ylabel('observation');
title("data");
if write_fig 
    filename = "../results/2fedata.pdf";
    print(fig, filename, '-dpdf','-r300');
    close;
end

% Mean function
meanfunc = {@meanZero};
hyp.mean = [];

% Covariance function and priors for hyperparameters
% cov([i,t],[j,t']) = [i==j]*unit_std**2 + cov(t,t')
cov_unit = {@covMask, {[1,0], {@covSEiso}}};
cov_time = {@covMask, {[0,1], {@covSEiso}}};
covfunc = {@covSum, {cov_unit,cov_time}};
unit_ls = 0.01; % ls of 0.01 basically means any i!=j will have 0 covariance
% The isotropic squared exponential function has two hyperparamters.
hyp.cov = log([unit_ls; unit_std; time_ls; time_os]);

% Likelihood function for the sample and priors for hyperparameters
likfunc = {@likGauss}; 
hyp.lik = [log(noise)];

% Define hyperpriors
prior.cov = {@priorDelta,... % fixed unit ls
                {@priorTransform,@exp,@exp,@log,{@priorGamma,1,2}}, ... % Gamma prior on unit os
                {@priorTransform,@exp,@exp,@log,{@priorGamma,1,1}}, ... % Gamma prior on time ls
               {@priorTransform,@exp,@exp,@log,{@priorGamma,1,2}}}; % Gamma prior on time os
prior.mean = {}; % no mean hyp
prior.lik = {{@priorTransform,@exp,@exp,@log,{@priorGamma,0.01,10}}}; % Gamma prior on noise

inffunc = {@infPrior, @infExact, prior};
p.method = 'LBFGS';
p.length = 100;

% Find the maximum a posteriori estimates for the model hyperparameters.
hyp = minimize_v2(hyp, @gp, p, inffunc, meanfunc, covfunc, likfunc, x, y);
fprintf("MAP unit ls: %.3f\n", exp(hyp.cov(1)));
fprintf("MAP unit os: %.3f\n", exp(hyp.cov(2)));
fprintf("MAP time ls: %.3f\n", exp(hyp.cov(3)));
fprintf("MAP time os: %.3f\n", exp(hyp.cov(4)));
fprintf("MAP noise: %.3f\n", exp(hyp.lik));

% Report negative log marginal likelihood of the data with the optimized
% model hyperparameters.
nlml2 = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x, y);
% Make predictions
[m, s2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x, y, x);

% Plot predicted values and credible intervals for observation
figure(3);
set(gca, 'FontSize', 24);
for i=1:n
    idx = ((i-1)*T+1):(i*T);
    f = [m(idx)+2*sqrt(s2(idx)); flipdim(m(idx)-2*sqrt(s2(idx)),1)];
    fill([(1:T)'; (T:-1:1)'], f, [7 7 7]/8);
    hold on; plot((1:T)', m(idx), 'LineWidth', 2); plot((1:T)', y(idx), '+', 'MarkerSize', 12);
    hold on;
end
grid on;
xlabel('time');
ylabel('y');
title("MAP");
if write_fig 
    filename = "../results/2feMAP.pdf";
    print(fig, filename, '-dpdf','-r300');
    close;
end

% Plot predicted values and credible intervals for time effect
% basically this is posterior of f conditioning on observing f+g+noise
prior_time_cov = feval(cov_time{:}, hyp.cov(3:4), x(1:T,:));
time_y_cov = feval(cov_time{:}, hyp.cov(3:4), x(1:T,:), x);
y_cov = feval(covfunc{:}, hyp.cov, x) + exp(2*hyp.lik)*eye(T*n);

post_time_mu = time_y_cov*(y_cov\y);
post_time_cov = prior_time_cov - time_y_cov*(y_cov\time_y_cov');

figure(4);
set(gca, 'FontSize', 24);
f = [post_time_mu+2*sqrt(diag(post_time_cov));...
    flipdim(post_time_mu-2*sqrt(diag(post_time_cov)),1)];
fill([(1:T)'; (T:-1:1)'], f, [7 7 7]/8);
hold on; plot((1:T)', post_time_mu, 'LineWidth', 2); plot((1:T)', time_effects, '+', 'MarkerSize', 12);
grid on;
xlabel('time');
ylabel('effect');
title("time random effect MAP");
if write_fig 
    filename = "../results/timeREMAP.pdf";
    print(fig, filename, '-dpdf','-r300');
    close;
end

% Plot predicted values and credible intervals for unit effect
% basically this is posterior of f conditioning on observing f+g+noise
prior_unit_cov = feval(cov_unit{:}, hyp.cov(1:2), x(1:T:n*T,1));
unit_y_cov = feval(cov_unit{:}, hyp.cov(1:2), x(1:T:n*T,1), x(:,1));
y_cov = feval(covfunc{:}, hyp.cov, x) + exp(2*hyp.lik)*eye(T*n);

post_unit_mu = unit_y_cov*(y_cov\y);
post_unit_cov = prior_unit_cov - unit_y_cov*(y_cov\unit_y_cov');

figure(5);
set(gca, 'FontSize', 24);
f = [post_unit_mu+2*sqrt(diag(post_unit_cov));...
    flipdim(post_unit_mu-2*sqrt(diag(post_unit_cov)),1)];
fill([(1:n)'; (n:-1:1)'], f, [7 7 7]/8);
hold on; plot((1:n)', post_unit_mu, 'LineWidth', 2); plot((1:n)', unit_effects, '+', 'MarkerSize', 12);
grid on;
xlabel('time');
ylabel('effect');
title("unit fixed effect MAP");
if write_fig 
    filename = "../results/unitFEMAP.pdf";
    print(fig, filename, '-dpdf','-r300');
    close;
end
