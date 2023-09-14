data {
  int<lower=1> N; // number of respondents
  int<lower=1> M; // number of itemss
  int<lower=1> H; // number of horizons
  int<lower=1> C; // number of response categories
  int<lower=1,upper=C>[N,M,H] y; // ordinal responses
}

transformed data {
  real delta = 1e-9; // jitter term on covariance diagonals
  real rho = 1; // output scale of IRFs and latent positions
  real l_x = H/2; // length scale of latent positions
  real l_f = 1; // length scale of IRFs
  real
}

parameters {
  real x[N,H]; // latent positions
  vector[N] eta;
}

model {
  vector[N] f;
  {
    matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_quad(x, alpha, rho);

    // diagonal elements
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;

    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }

  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();
  eta ~ std_normal();

  y ~ normal(f, sigma);
}
