data {
  int N;
  vector[N] d;
  int n[N];
  int y[N];
  real mu;
  real<lower=0> sigma;
  real<lower=0> lambda;
}
transformed data {
  vector[N-1] dd = d[2:N] - d[1:(N-1)];

}
parameters {
  real eta1;
  vector<lower=0>[N - 1] u;
  vector<lower=0>[N - 1] tau;
  real<lower=0> gamma;
}
transformed parameters {
  vector[N] eta;
  vector[N] theta;
  eta[1] = eta1;
  for(i in 2:N) {
    eta[i] = eta[i - 1] + u[i-1];
  }
  theta = inv_logit(eta);
}
model {
  eta1 ~ normal(mu, sigma);
  gamma ~ cauchy(0, lambda);
  tau ~ cauchy(0, gamma);
  u ~ normal(0, tau);
  y ~ binomial_logit(n, eta);
}
generated quantities {
  vector[N-1] zeta = theta[1:(N-1)] - theta[N];
}
