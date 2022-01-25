import math
import torch
import gpytorch
from matplotlib import pyplot as plt
from gpytorch.models import ExactGP
from gpytorch.likelihoods import GaussianLikelihood
from gpytorch.means import ConstantMean, ZeroMean
from gpytorch.kernels import ScaleKernel, RBFKernel

N=5 # num of items
M=10 # num of respondents 

true_scores=2*torch.normal(mean=0,std=1,size=(1,M)).reshape((-1,))-1 # true scores [-1,1]

# define IRFs
theta_grids = torch.linspace(-2,2,10)

class GPModel(ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super().__init__(train_x, train_y, likelihood)
        self.mean_module = ConstantMean()
        self.covar_module = ScaleKernel(RBFKernel())

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)

likelihood = gpytorch.likelihoods.GaussianLikelihood()
model = GPModel(theta_grids, torch.sin(theta_grids), likelihood)
model.covar_module.outputscale = 3.0
model.covar_module.base_kernel.lengthscale = 1.0
model.likelihood.noise=0.01

true_fs = torch.zeros(size=(N,len(theta_grids)))
output=model(theta_grids)
for i in range(N):
    true_fs[i]=output.sample()

xs = torch.linspace(-2,2,100)
def IRF(xs, i):
    irtlikelihood = gpytorch.likelihoods.GaussianLikelihood()
    irtmodel = GPModel(theta_grids, true_fs[i], irtlikelihood)
    irtmodel.covar_module.outputscale = 3.0
    irtmodel.covar_module.base_kernel.lengthscale = 1.0
    irtmodel.eval()
    irtlikelihood.eval()
    with torch.no_grad():
        f_preds = irtmodel(xs)
    return f_preds.mean.detach()

# plot true IRFs
for i in range(N):
    plt.plot(xs, torch.sigmoid(IRF(xs, i)))
plt.ylim([0,1])
plt.show()

# define data
train_x = torch.zeros(N*M, 3) # (item, respondent, score)
train_y = torch.zeros(N*M) # (response)
for i in range(N):
    train_y[(M*i):(M*(i+1))] = torch.sigmoid(IRF(true_scores, i))
    for j in range(M):
        train_x[M*i+j, 0] = i
        train_x[M*i+j, 1] = j
        train_x[M*i+j, 2] = true_scores[j]

# define GPIRT model

class GPIRTModel(ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super().__init__(train_x, train_y, likelihood)
        self.mean_module = ZeroMean()
        self.covar_module = RBFKernel()

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
