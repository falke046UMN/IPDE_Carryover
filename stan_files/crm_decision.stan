
data {
  int<lower=0> D;
  array[D] int y;
  array[D] int n;
  array[D] real d;
}

parameters {
  real theta;
}

model {
  // prior
  theta ~ normal(0,sqrt(2));

  // model
  y ~ binomial(n, d.^exp(theta));
}
