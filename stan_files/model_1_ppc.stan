data {
    int N;  # number of samplings
    int rating_difference[N];
    int time_control[N];

}

generated quantities {
    real gamma_white = normal_rng(-0.0015, 0.0001);
    real gamma_black = normal_rng(0.0015, 0.0001);

    real time_controll_coeff[3]; // 1 - blitz, 2 - rapid, 3 - classic
    time_controll_coeff[1] = normal_rng(1.3, 0.1);
    time_controll_coeff[2] = normal_rng(0.8, 0.1);
    time_controll_coeff[3] = normal_rng(0.4, 0.1);

    int<lower=0, upper=50> mistakes_white[N];
    int<lower=0, upper=50> mistakes_black[N];

    for (ind in 1:N){
        mistakes_white[ind] = poisson_log_rng(time_controll_coeff[time_control[ind]] + gamma_white * rating_difference[ind]);
        mistakes_black[ind] = poisson_log_rng(time_controll_coeff[time_control[ind]] + gamma_black * rating_difference[ind]);
    }
}
