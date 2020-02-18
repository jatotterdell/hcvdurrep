data {
  int N;
  vector[N] d;
  int n[N];
  int y[N];
  real mu;
  real<lower=0> sigma;
  real<lower=0> tau_n;
  real<lower=0> tau_m;
}
transformed data {
  vector[N-1] dd = d[2:N] - d[1:(N-1)];

}
parameters {
  real eta1;
  vector<lower=0>[N - 1] u;
  real<lower=0> tausq;
}
transformed parameters {
  vector[N] eta;
  vector[N] theta;
  real<lower=0> tau = sqrt(tausq);
  eta[1] = mu + sigma*eta1;
  for(i in 2:N) {
    eta[i] = eta[i - 1] + u[i-1];
  }
  theta = inv_logit(eta);
}
model {
  target += normal_lpdf(eta1 | 0, 1)
          + normal_lpdf(u | 0, tau) - normal_lccdf(0 | 0, tau)
          + inv_gamma_lpdf(tau | tau_n * inv(2), tau_m^2 * tau_n * inv(2))
          + binomial_logit_lpmf(y | n, eta);
}
generated quantities {
  vector[N-1] zeta = theta[1:(N-1)] - theta[N];
}
