addpath("../results");

fig = figure(1);
tiledlayout(2,2,'Padding', 'none', 'TileSpacing', 'compact');

FONTSIZE = 16;

nexttile;
x = linspace(-3,3,100);
y = 1./(1+exp(-3*(x-1)));
plot(x,y);
ylim([0,1]);
set(gca,'TickLength',[0 0], 'ytick', 0:0.2:1);
xlabel("{\theta}", 'FontSize', FONTSIZE);
ylabel("Pr(y = 1)", 'FontSize', FONTSIZE);
title("Monotonic IRF", 'FontSize', FONTSIZE);

nexttile;
x = linspace(-3,3,100);
y0 = exp(-2*0*(x));
y1 = exp(-2*1*(x)+2);
y2 = exp(-2*2*(x));
y3 = exp(-2*3*(x));
y = (y1+y2)./(y0+y1+y2+y3);
plot(x,y);
ylim([0,1]);
set(gca,'TickLength',[0 0], 'ytick', 0:0.2:1);
xlabel("{\theta}", 'FontSize', FONTSIZE);
ylabel("Pr(y = 1)", 'FontSize', FONTSIZE);
title("Non-monotonic IRF", 'FontSize', FONTSIZE);

nexttile;
x = linspace(-3,3,100);
y = 1./(1+exp(-3*(x+1)));
plot(x,y);
ylim([0,1]);
set(gca,'TickLength',[0 0], 'ytick', 0:0.2:1);
xlabel("{\theta}", 'FontSize', FONTSIZE);
ylabel("Pr(y = 1)", 'FontSize', FONTSIZE);
title("Saturate IRF", 'FontSize', FONTSIZE);

nexttile;
x = linspace(-3,3,100);
y = 0.1 + 0.1./(1+exp(-3*(x)));
plot(x,y);
ylim([0,1]);
set(gca,'TickLength',[0 0], 'ytick', 0:0.2:1);
xlabel("{\theta}", 'FontSize', FONTSIZE);
ylabel("Pr(y = 1)", 'FontSize', FONTSIZE);
title("Non-saturate IRF", 'FontSize', FONTSIZE);


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
covfunc = {@covSEiso};              % Squared Exponental covariance function
likfunc = {@likGauss};              % Gaussian likelihood
hyp = struct('mean', [0], 'cov', [0 0], 'lik', log(0.1));

rng(12345);
x = [0;1;2;3;4;5;6;8;9;10];
y = [2;1.5;1;1.2;0.7;0.5;0.7;1.2;1.2;1.5];
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
scatter(x,y,24,'ro', 'filled'); hold on;

[~,~,mu, s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y, xs);
f = [mu+1.92*sqrt(s2); flipdim(mu-1.92*sqrt(s2),1)];
fill([xs; flipdim(xs,1)], f, [64, 214, 247] / 255, ...
     'facealpha', 0.5, ...
     'edgecolor', 'none');
plot(xs,mu,'m', 'LineWidth',3);
plot(xs,target_f,'k--','LineWidth',3); 
 
 xlabel("{\theta}", 'FontSize', FONTSIZE);
 ylabel("f({\theta})", 'FontSize', FONTSIZE);
 legend('training data', '2{\sigma} credible region',...
     'prediction' ,'target function','Location', 'Best', 'FontSize', FONTSIZE, 'NumColumns', 2);

set(fig, 'PaperPosition', [0 0 12 4]); 
set(fig, 'PaperSize', [12 4]); 

filename = "./figures/gp.pdf";
print(fig, filename, '-dpdf','-r300', '-fillpage');
close;

fig = figure(1);

file_name = '../results/mq_dynamic.csv';
opts = detectImportOptions(file_name);
mq = readtable(file_name);
mq = sortrows(mq, ["year"]);

liberal = mq.avgscore(strcmp(mq.ideology,"liberal"), :);
conservative = mq.avgscore(strcmp(mq.ideology,"conservative"), :);
years = unique(mq.year);

plot(liberal, years); hold on;
plot(conservative, years);

