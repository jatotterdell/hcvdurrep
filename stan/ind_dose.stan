data {
  int N;
  vector[N] d;
  int n[N];
  int y[N];
  real mu;
  real<lower=0> sigma;
}
parameters {
  vector[N] eta;
}
transformed parameters {
  vector[N] theta;
  theta = inv_logit(eta);
}
model {
  target += normal_lpdf(eta | mu, sigma)
          + binomial_logit_lpmf(y | n, eta);
}
generated quantities {
  vector[N-1] zeta = theta[1:(N-1)] - theta[N];
}
