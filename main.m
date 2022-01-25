% add gpml path
gpml_path = "/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/gpml-matlab-v3.6-2015-07-07";
addpath(gpml_path);
addpath("utils");
startup;
rng("default");

% synthetic data
N=1; % num of items
M=200; % num of respondents
scores = normrnd(0,1,[M,1]);
% scores = linspace(-5,5,M)';

% define item response functions
meanfunc = {@meanZero};
covfunc = {@covSEiso};
likfunc = {@likGauss};
hyp = struct('mean', [], 'cov', [0, log(1)], 'lik', -10);

% define grid of scores
x = linspace(-2,2,10);
xs = -5:0.01:5;
% score_pdfs = normpdf(xs, 0, 1);
% scores = score_grid(gendist(score_pdfs, M, 1))';

mu = feval(meanfunc{:},hyp.mean,x');
sigma = feval(covfunc{:},hyp.cov, x');
y = mvnrnd(mu, sigma, N);

[~,~, ys, ~] = gp(hyp, @infExact, meanfunc,...
   covfunc, @likGauss, x', y(1,:)', xs');

plot(xs,sigmoid_fn(ys)); hold on;

% define training data (item, response)
train_x = [reshape(repmat([1:N], M,1), [],1),...
    repmat((1:M)',N,1)];

train_f = zeros(N,M);
for i=1:N
    [~,~, train_f(i,:), ~] = gp(hyp, @infExact, meanfunc,...
        covfunc, @likGauss, x', y(i,:)', scores);
%     train_f((i-1)*M+1:i*M) = ys;
end

% bernoulli likelihood after sigmoid transform
% p(y=1|f) = sigmoid(f) 
train_y = sigmoid_fn(train_f);
tmp = rand(size(train_y));
train_y = (tmp<train_y);

% bernoulli post
hyp.lik = [];
[~,~,bmu, ~] = gp(hyp, @infEP, meanfunc,...
                covfunc, @likLogistic, scores, 2*train_y'-1, xs');

plot(xs,sigmoid_fn(bmu)); hold on;
hyp.lik = -10;

% define prior for latent scores
prior_score = struct;
prior_score.mu = 0;
prior_score.cov = 1;

% priors on fs are the same
prior_fs = struct;
prior_fs.mu = zeros(M,1);
prior_fs.cov = feval(covfunc{:},hyp.cov, scores);

% initialize sampler
% cur_scores = score_grid(gendist(score_pdfs, M, 1))';
cur_scores = zeros(M,1);
for j=1:M
   cur_scores(j) = normrnd(prior_score.mu, prior_score.cov);
end
cur_scores = scores;
cur_fs = zeros(N,M);
cur_fstars = zeros(N,length(xs));

for i=1:N
   cur_fs(i,:) = mvnrnd(prior_fs.mu, prior_fs.cov);
end

% Gibbs sampling
ITERS=1000;
lls = zeros(ITERS,1);
fs_samples = zeros(ITERS, M);
score_samples = zeros(ITERS, M);

for iter=1:ITERS
    % sample new fs given scores, fs and ys
    for i=1:N
        [cur_fs(i,:), ll] = sample_fs(cur_fs(i,:), prior_fs, train_y(i,:));  
        fs_samples(iter, :) = cur_fs(i,:);
    end
    
    % plot the evolution
%     figure(1);
%     h1 = plot(xs,sigmoid_fn(ys)); hold on;
%     set(gcf, 'Color', 'w');
%     h2 = scatter(cur_scores, sigmoid_fn(cur_fs(i,:)), '*');
%     title(sprintf('iter %d', iter));
%     drawnow;
%     
%     
%     delete(h2);
    
    % sample new scores given fs, scores and ys
%     for j=1:M
%         [cur_scores(j),ll] = sample_score(cur_scores(j), prior_score, cur_fs(:,j), train_y(:,j));    
%     end
%     score_samples(iter, :) = cur_scores;
    
    lls(iter) = ll;
end

scatter(cur_scores, mean(sigmoid_fn(fs_samples(900:end,:))));

% plot(1:ITERS, lls);

% extend fs to fstars
% for i=1:N
%     [~,~, cur_fstars(i,:), ~] = gp(hyp, @infExact, meanfunc,...
%             covfunc, @likGauss, cur_scores, cur_fs(i,:)', xs');
%     % plot irf
%     plot(xs,sigmoid_fn(cur_fstars(i,:)));
% end

legend({'true IRF', 'post', 'sampled IRF'},'Location','southwest');
