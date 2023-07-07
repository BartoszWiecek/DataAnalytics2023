data {
    int N;  # number of samplings
    int rating_difference[N];
    int time_control[N];
    int opening[N];
    int mistakes_white_model[N];
    int mistakes_black_model[N];

}

parameters {
    real<lower=-0.0015, upper=-0.0013> gamma_white;
    real<lower=0.0013, upper=0.00155> gamma_black;
    real<lower=1.26, upper=1.27> time_controll_coeff_blitz;
    real<lower=0.45, upper=0.48> time_controll_coeff_rapid;
    real<lower=0.35, upper=0.37> time_controll_coeff_classic;
    real<lower=0.571874, upper=0.591874> opening_coeff_sicilian;
    real<lower=0.349679, upper=0.369679> opening_coeff_petroff;
    real<lower=0.317872, upper=0.337872> opening_coeff_queen_pawn;
    real<lower=0.317622, upper=0.337622> opening_coeff_italian;

}

transformed parameters {
    real time_controll_coeff[3]; // 1 - blitz, 2 - rapid, 3 - classic
    real opening_coeff[4]; // 1 - Sicilian Defense, 2 - Petroff Defense, 3 - Queen Pawn Game, 4 - Italian Game
    real gammas[2]; // 1 - gamma_white, 2 - gamma_black

    gammas[1] = gamma_white;
    gammas[2] = gamma_black;
    time_controll_coeff[1] = time_controll_coeff_blitz;
    time_controll_coeff[2] = time_controll_coeff_rapid;
    time_controll_coeff[3] = time_controll_coeff_classic;

        
    opening_coeff[1] = opening_coeff_sicilian;
    opening_coeff[2] = opening_coeff_petroff;
    opening_coeff[3] = opening_coeff_queen_pawn;
    opening_coeff[4] = opening_coeff_italian;
}

model {
    gammas[1] ~ normal(-0.00141, 0.0001);
    gammas[2] ~ normal(0.00154, 0.0001); 

    time_controll_coeff[1] ~ normal(1.281, 0.01);
    time_controll_coeff[2] ~ normal(0.840, 0.01);
    time_controll_coeff[3] ~ normal(0.269, 0.01);

    opening_coeff[1] ~ normal(0.675, 0.05);
    opening_coeff[2] ~ normal(0.305, 0.05);
    opening_coeff[3] ~ normal(0.389, 0.05);
    opening_coeff[4] ~ normal(0.288, 0.05);

    for (ind in 1:N){
        mistakes_white_model[ind] ~ poisson_log(time_controll_coeff[time_control[ind]] + opening_coeff[opening[ind]] + gammas[1] * rating_difference[ind]);
        mistakes_black_model[ind] ~ poisson_log(time_controll_coeff[time_control[ind]] + opening_coeff[opening[ind]] + gammas[2] * rating_difference[ind]);
    }
}


generated quantities {
    int<lower=0, upper=50> mistakes_white_pred[N];
    int<lower=0, upper=50> mistakes_black_pred[N];
    vector [N] log_lik;

    for (ind in 1:N){
        mistakes_white_pred[ind] = poisson_log_rng(time_controll_coeff[time_control[ind]] + opening_coeff[opening[ind]] + gammas[1] * rating_difference[ind]);
        mistakes_black_pred[ind] = poisson_log_rng(time_controll_coeff[time_control[ind]] + opening_coeff[opening[ind]] + gammas[2] * rating_difference[ind]);
        log_lik[ind] = poisson_log_lpmf(mistakes_white_model[ind] | time_controll_coeff[time_control[ind]] + opening_coeff[opening[ind]] + gammas[1] * rating_difference[ind]);
    }
}
