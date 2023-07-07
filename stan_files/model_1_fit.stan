data {
    int N;  # number of samplings
    int rating_difference[N];
    int time_control[N];
    int mistakes_white_model[N];
    int mistakes_black_model[N];

}

parameters {
    real<lower=-0.0016, upper=-0.0014> gamma_white;
    real<lower=0.0014, upper=0.0016> gamma_black;
    real<lower=1.13, upper=1.14> time_controll_coeff_blitz;
    real<lower=0.76, upper=0.78> time_controll_coeff_rapid;
    real<lower=0.20, upper=0.22> time_controll_coeff_classic;
}

transformed parameters {
    real time_controll_coeff[3]; // 1 - blitz, 2 - rapid, 3 - classic
    real gammas[2];

    gammas[1] = gamma_white;
    gammas[2] = gamma_black;
    time_controll_coeff[1] = time_controll_coeff_blitz;
    time_controll_coeff[2] = time_controll_coeff_rapid;
    time_controll_coeff[3] = time_controll_coeff_classic;
}

model {
    gammas[1] ~ normal(-0.00145, 0.0001);
    gammas[2] ~ normal(0.00137, 0.0001); 

    time_controll_coeff[1] ~ normal(1.234, 0.01);
    time_controll_coeff[2] ~ normal(0.836, 0.01);
    time_controll_coeff[3] ~ normal( 0.295, 0.01);
    
    for (ind in 1:N){
        mistakes_white_model[ind] ~ poisson_log(time_controll_coeff[time_control[ind]] + gammas[1] * rating_difference[ind]);
        mistakes_black_model[ind] ~ poisson_log(time_controll_coeff[time_control[ind]] + gammas[2] * rating_difference[ind]);
    }
}


generated quantities {
    int<lower=0, upper=50> mistakes_white_pred[N];
    int<lower=0, upper=50> mistakes_black_pred[N];
    vector [N] log_lik;

    for (ind in 1:N){
        mistakes_white_pred[ind] = poisson_log_rng(time_controll_coeff[time_control[ind]] + gammas[1] * rating_difference[ind]);
        mistakes_black_pred[ind] = poisson_log_rng(time_controll_coeff[time_control[ind]] + gammas[2] * rating_difference[ind]);
        log_lik[ind] = poisson_log_lpmf(mistakes_white_model[ind] | time_controll_coeff[time_control[ind]] + gammas[1] * rating_difference[ind]);
    }
}
