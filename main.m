clear; close all;

% add gpml path
gpml_path = "/Users/yahoo/Documents/WashU/CSE515T/Code/Gaussian Process/gpml-matlab-v3.6-2015-07-07";
addpath(gpml_path);
addpath("utils");
startup;
rng("default");

% synthetic data
synthetic_data;

% plot the first true IRF
for i=1:2
    figure(i);
    plot(xs,sigmoid_fn(true_fs(i,:))); hold on;
    ylim([0,1]);
end

% plot first true fs
for i=1:2
    figure(i);
    scatter(scores, sigmoid_fn(train_f(i,:)));
    scatter(scores, train_y(i,:), 'x');
end

% initialize sampler
cur_scores = sort(normrnd(0,1,M,1));
prior_fs.cov = feval(covfunc{:},hyp.cov, cur_scores);
cur_fs = zeros(N,M);

for i=1:N
   
   cur_fs(i,:) = mvnrnd(prior_fs.mu, prior_fs.cov);
end

% Gibbs sampling
WARMUP = 100;
SAMPLE = 100;
ITERS = WARMUP + SAMPLE;
fs_samples = zeros(ITERS, N, M);
fstar_samples = zeros(ITERS, N, numel(xs));
score_samples = zeros(ITERS, M);
lls = zeros(ITERS,N);
h = waitbar(0,'Please wait...');
for iter=1:ITERS 
    waitbar(iter / ITERS);
    % sample new fs given scores, fs and ys
    
    for i=1:N
        [cur_fs(i,:), lls(iter,i)] = sample_fs(cur_fs(i,:), prior_fs, train_y(i,:)); 
    end

    % get fstar
    fstars = zeros(N, numel(xs));
    for i=1:N
        [~,~, fstars(i,:), ~] = gp(hyp, @infExact, meanfunc,...
            covfunc, @likGauss, cur_scores, cur_fs(i,:)', xs);
        fstar_samples(iter, i,:) = fstars(i,:);
    end
    
    % sample new scores given fs, scores and ys
    cur_scores = sample_score(xs, fstars, train_y, prior_score);
    score_samples(iter, :) = cur_scores;
    
    % update fs cov based on current score
    prior_fs.cov = feval(covfunc{:},hyp.cov, cur_scores);
    
    idx = int64((cur_scores-min(xs))/0.01)+1;
    for i=1:N
        cur_fs(i,:) = fstars(i,idx);
        fs_samples(iter, i, :) = cur_fs(i,:);
    end
   
     % plot the evolution
     for i=1:2
        if mod(iter,5)==1 && 1
            figure(i);
            if exist('h2', 'var'), delete(h1);delete(h2); end
            set(gcf, 'Color', 'w');
            h1 = scatter(cur_scores, sigmoid_fn(cur_fs(i,:)), '*');
            h2 = plot(xs, sigmoid_fn(fstars(i,:)));
            drawnow;
            pause(0.5);
        end
     end
end
close(h);

fs_samples = fs_samples((end-SAMPLE):end,:,:);
fstar_samples = fstar_samples((end-SAMPLE):end,:,:);
score_samples = score_samples((end-SAMPLE):end,:);

for i=1:2
    figure(i);

    scatter(mean(score_samples,1), sigmoid_fn(mean(fs_samples(:,i,:),1)));
    tmp = mean(fstar_samples(:,i,:),1);
    plot(xs, sigmoid_fn(reshape(tmp,numel(xs),1)));
    legend({'true IRF', 'true fs', 'response', 'fs mean', 'fstar'},'Location','southwest');
end
