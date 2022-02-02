% define gp model
% meanfunc = {@meanSum, {@meanConstant, @meanPoly}};
meanfunc = {@meanZero};
covfunc = {@covSEiso};
likfunc = {@likGauss};
hyp = struct('mean', [], 'cov', [log(1)/2, log(1)], 'lik', log(0.01));
