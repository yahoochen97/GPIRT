addpath("../results");

fig = figure(1);
tiledlayout(2,2,'Padding', 'none', 'TileSpacing', 'compact');

FONTSIZE = 16;

nexttile;
x = linspace(-3,3,100);
y = 1./(1+exp(-3*(x+1)));
plot(x,y);
ylim([0,1]);
set(gca,'TickLength',[0 0], 'ytick', 0:0.2:1);
xlabel("{\theta}", 'FontSize', FONTSIZE);
ylabel("Pr(y = 1)", 'FontSize', FONTSIZE);
title("Saturate and monotonic IRF", 'FontSize', FONTSIZE);

nexttile;
x = linspace(-3,3,100);
y0 = exp(-2*0*(x)-2);
y1 = exp(-2*1*(x)+5);
y2 = exp(-2*2*(x)+5);
y3 = exp(-2*3*(x)+2);
y = (y1+y2)./(y0+y1+y2+y3);
plot(x,y);
ylim([0,1]);
set(gca,'TickLength',[0 0], 'ytick', 0:0.2:1);
xlabel("{\theta}", 'FontSize', FONTSIZE);
ylabel("Pr(y = 1)", 'FontSize', FONTSIZE);
title("Non-monotonic but saturate", 'FontSize', FONTSIZE);


nexttile;
x = linspace(-3,3,100);
y = 0.2+0.6*1./(1+exp(-3*(x-1)));
plot(x,y);
ylim([0,1]);
set(gca,'TickLength',[0 0], 'ytick', 0:0.2:1);
xlabel("{\theta}", 'FontSize', FONTSIZE);
ylabel("Pr(y = 1)", 'FontSize', FONTSIZE);
title("Monotonic but non-saturate", 'FontSize', FONTSIZE);

nexttile;
x = linspace(-3,3,100);
y = 0.1 + 0.1./(1+exp(-3*(x)-3*sin(x*3)));
plot(x,y);
ylim([0,1]);
set(gca,'TickLength',[0 0], 'ytick', 0:0.2:1);
xlabel("{\theta}", 'FontSize', FONTSIZE);
ylabel("Pr(y = 1)", 'FontSize', FONTSIZE);
title("Non-saturate and non-monotonic", 'FontSize', FONTSIZE);


% nexttile;
% x = linspace(-3,3,100);
% y = 1./(1+exp(-3*(x)));
% plot(x,y);
% ylim([0,1]);
% set(gca,'TickLength',[0 0], 'ytick', 0:0.2:1);
% xlabel("{\theta}", 'FontSize', FONTSIZE);
% ylabel("Pr(y = 1)", 'FontSize', FONTSIZE);
% title("Symmetric IRF", 'FontSize', FONTSIZE);
% 
% nexttile;
% x = linspace(-3,3,100);
% y = 0.8./(1+exp(-2*(x+1)));
% plot(x,y);
% ylim([0,1]);
% set(gca,'TickLength',[0 0], 'ytick', 0:0.2:1);
% xlabel("{\theta}", 'FontSize', FONTSIZE);
% ylabel("Pr(y = 1)", 'FontSize', FONTSIZE);
% title("Asymmetric IRF", 'FontSize', FONTSIZE);

set(fig, 'PaperPosition', [0 0 10 4]); 
set(fig, 'PaperSize', [10 4]); 

filename = "./figures/irf.pdf";
print(fig, filename, '-dpdf','-r300', '-fillpage');
close;


gpml_path = "/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/gpml-matlab-v3.6-2015-07-07";
addpath(gpml_path);
startup;

meanfunc = {@meanConst};            % constant mean
covfunc = {@covMaterniso,5};              % Squared Exponental covariance function
likfunc = {@likGauss};              % Gaussian likelihood
hyp = struct('mean', [0], 'cov', [0 0], 'lik', log(0.1));

rng(12345);
x = [0;1;2;3;4;5;6;8;9;10];
y = [2;1.5;1;1.2;0.7;0.5;0.7;1.2;1.2;1.5]-1;
n = 50;
noise = 0.1;
xs = linspace(0,10,n)';
% mu = feval(meanfunc{:}, hyp.mean, xs);
% K = feval(covfunc{:}, hyp.cov, xs);
% target_f = mvnrnd(mu,K,1);
[target_f, ~] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y, xs);

x = [7*rand(n*0.8,1);8+2*rand(n*0.2,1)];
[y, ~] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, xs, target_f, x);
y = y + noise*normrnd(0,1,n,1);

fig = figure(1);
% scatter(x,y,24,'ro', 'filled'); hold on;

[~,~,mu, s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y, xs);
f = [mu+1.92*sqrt(s2); flipdim(mu-1.92*sqrt(s2),1)];
fill([xs; flipdim(xs,1)], f, [64, 214, 247] / 255, ...
     'facealpha', 0.5, ...
     'edgecolor', 'none'); hold on;
plot(xs,mu,'m', 'LineWidth',3);
% plot(xs,target_f,'k--','LineWidth',3); 
 

title("GP dynamic latent traits", 'FontSize', FONTSIZE);
 xlabel("t", 'FontSize', FONTSIZE);
 ylabel("{\theta}", 'FontSize', FONTSIZE);
%  legend('training data', '2{\sigma} credible region',...
%      'prediction' ,'target function','Location', 'Best', 'FontSize', FONTSIZE, 'NumColumns', 2);

legend('2{\sigma} credible region',...
     'posterior mean','Location', 'Best', 'FontSize', FONTSIZE, 'NumColumns', 2);


set(fig, 'PaperPosition', [-0.5 0.1 9 4]); 
set(fig, 'PaperSize', [7.8 4.1]); 

filename = "./figures/gpexample.pdf";
print(fig, filename, '-dpdf','-r300');
close;

fig = figure(1);

% 1PL
x = linspace(-3,3,100);
y = 1./(1+exp(-2*x));
plot(x,y, 'b-.'); hold on;


% 3PL
x = linspace(-3,3,100);
y = 0.2+0.8./(1+exp(-3*(x-1)));
plot(x,y, 'r--');

% GGUM
y0 = exp(-2*0*(x)-2);
y1 = exp(-2*1*(x)+2);
y2 = exp(-2*2*(x)+2);
y3 = exp(-2*3*(x)-2);
y = (y1+y2)./(y0+y1+y2+y3);
plot(x,y, 'k:');


ylim([0,1]);
set(gca,'TickLength',[0 0], 'ytick', 0:0.2:1);
xlabel("{\theta}", 'FontSize', FONTSIZE);
ylabel("Pr(y = 1)", 'FontSize', FONTSIZE);
title("Standard IRT models", 'FontSize', FONTSIZE);

legend('Rasch','3PL','GGUM','Location', 'Best', 'FontSize', FONTSIZE, 'NumColumns', 1);

set(fig, 'PaperPosition', [-0.5 -0.2 9 6]); 
set(fig, 'PaperSize', [7.9 5.8]); 

filename = "./figures/standardirfs.pdf";
print(fig, filename, '-dpdf','-r300');
close;

fig = figure(1);

addpath("~/Documents/Washu/CSE515T/Code/Gaussian Process/gpml-matlab-v3.6-2015-07-07");
startup;

x = linspace(-4,4,21)';
xs = linspace(-3,3,101)';
f = 4*x;
p = normcdf(f);
y = 2*binornd(1,p)-1;

meanfunc = {@meanZero};   
covfunc = {@covSEiso};             
likfunc = {@likErf};
inffunc = {@infEP};
% hyp.mean = [-2;1;1];
hyp.mean = [];
hyp.cov = [log(1);log(2)];
hyp.lik = [];

prior.cov  = {@priorDelta, @priorDelta};  
prior.lik  = {};
prior.mean = {};

[~,~,fmu,fs2] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, x, y, xs);

plot(xs,normcdf(fmu./sqrt(1+fs2)), 'r--'); hold on;


% 3PL
f = -x.^2+4;
p = normcdf(f);
y = 2*binornd(1,p)-1;
[~,~,fmu,fs2] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, x, y, xs);

plot(xs,normcdf(fmu./sqrt(1+fs2)), 'b-.');


% GGUM
f = x.^2/8+x/8-1;
p = normcdf(f);
y = 2*binornd(1,p)-1;
[~,~,fmu,fs2] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, x, y, xs);

plot(xs,normcdf(fmu./sqrt(1+fs2)), 'k:');


ylim([0,1]);
set(gca,'TickLength',[0 0], 'ytick', 0:0.2:1);
xlabel("{\theta}", 'FontSize', FONTSIZE);
ylabel("Pr(y = 1)", 'FontSize', FONTSIZE);
title("GP response functions", 'FontSize', FONTSIZE);

legend('Logistic','Nonmonotonic','Nonsaturate','Location', 'Best', 'FontSize', FONTSIZE, 'NumColumns', 1);


set(fig, 'PaperPosition', [-0.5 -0.2 9 6]); 
set(fig, 'PaperSize', [7.9 5.8]); 

filename = "./figures/gpirfs.pdf";
print(fig, filename, '-dpdf','-r300');
close;
