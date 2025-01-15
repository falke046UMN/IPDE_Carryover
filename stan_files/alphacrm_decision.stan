functions {
  real linear_interp(real x, real[] x_pred, real[] y_pred, int D){
    real out;
    int i = 1;
    for(k in 2:(D-1)){i = i + (x_pred[k] <= x);}
    out = y_pred[i] + (y_pred[i + 1] - y_pred[i]) * (x - x_pred[i])/(x_pred[i + 1] - x_pred[i]);
    out = fmin(out, 0.98);
    return(out);
  }
}

data {
  int<lower=0> D;
  array[D] int y;
  array[D] int y_ipde;
  array[D] int n;
  array[D] int n_ipde;  
  array[D] real d;
  array[D] real d1_mg;
  real<lower=0> alpha_a;
  real<lower=0> alpha_b;
}

transformed data{
  array[D] real d1_mg_prev;
  d1_mg_prev[1] = 0;
  for(i in 2:D){
    d1_mg_prev[i] = d1_mg[(i-1)];
  }
}

parameters {
  real theta;
  real<lower=0,upper=1> alpha;
}

transformed parameters {
  array[D] real<lower=0,upper=1> d_ipde;
  for(i in 1:D){
   d_ipde[i] = linear_interp(d1_mg[i] + alpha*d1_mg_prev[i], d1_mg, d, D);
  }
}

model {
  // prior
  theta ~ normal(0, sqrt(2));
  alpha ~ beta(alpha_a, alpha_b);

  // model
  y ~ binomial(n, d.^exp(theta));
  y_ipde ~ binomial(n_ipde, d_ipde.^exp(theta));
}

