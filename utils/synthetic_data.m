% synthetic data
N=10; % num of items
M=40; % num of respondents

% define item response functions
model;
coefs = unifrnd(-1,1, N, 3);

% define grid of scores
xs = (-3:0.01:3)';
scores = linspace(-2,2,M);

true_fs = zeros(N,length(xs));
for i=1:N
    true_fs(i,:) = score_fn(coefs, i, xs);
end

train_f = zeros(N,M);
for i=1:N
      train_f(i,:) = score_fn(coefs, i, scores);
end

% bernoulli likelihood after sigmoid transform
% p(y=1|f) = sigmoid(f) 
train_y = zeros(N,M);
for i=1:N
   for j=1:M
       train_y(i,j) = (rand()<sigmoid_fn(train_f(i,j)));
   end
end


% define standard normal prior for latent scores
prior_score = struct;
prior_score.mu = 0;
prior_score.cov = 1;


% priors on fs are the same
prior_fs = struct;
prior_fs.mu = zeros(M,1);
prior_fs.cov = feval(covfunc{:},hyp.cov, scores');


function thetas = score_fn(coefs, i, xs)
%     model;
%     [~,~,thetas,~] = gp(hyp, @infExact, meanfunc,...
%                     covfunc, likfunc, x', y', xs);
    thetas = coefs(i,1)+coefs(i,2)*xs+coefs(i,3)*xs.^2/4;
end