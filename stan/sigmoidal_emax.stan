data {
  int N;
  vector[N] d;
  int n[N];
  int y[N];
  vector[4] mu;
  vector<lower=0>[4] sigma;
}
parameters {
  real beta1;
  real beta2;
  real<lower=0> beta3;
  real<lower=0> beta4;
}
transformed parameters {
  vector[N] eta;
  vector[N] theta;
  for(i in 1:N)
    eta[i] = beta1 + beta2 * pow(d[i], beta4) / (pow(d[i], beta4) + pow(beta3, beta4));
  theta = inv_logit(eta);
}
model {
  target += normal_lpdf(beta1 | mu[1], sigma[1])
          + normal_lpdf(beta2 | mu[2], sigma[2])
          + normal_lpdf(beta3 | mu[3], sigma[3]) - normal_lccdf(0 | mu[3], sigma[3])
          + normal_lpdf(beta4 | mu[4], sigma[4]) - normal_lccdf(0 | mu[4], sigma[4])
          + binomial_logit_lpmf(y | n, eta);
}
generated quantities {
  vector[N-1] zeta = theta[1:(N-1)] - theta[N];
}
