functions {
    //Zellner (eqn 2.17 page 22)   
    real se_lpdf(vector std, vector n, real sig ) {
        return -sum((n+1))*log(sig)-0.5*sum(n.*std.*std)/(sig*sig);
    }
}
data {
    /* none of the cmip5 models have the dissicnat variable
       so we only know the difference, ΔCant - ΔCnat */
    int<lower=0> N1a;        // number of cmip5 models with nⱼ = 1
    real<lower=0> d1a[N1a];  // (ΔCant - ΔCnat) of j-th cmip5 models with nⱼ = 1
    
    int<lower=0> N1b;        // number of cmip5 models with nⱼ > 1
    vector[N1b] d1b;         // (ΔCant - ΔCnat) of j-th cmip5 model with nⱼ > 1
    vector[N1b] sig1b;       // std(ΔCant - ΔCnat) of j-th cmip5 model
    vector[N1b] n1b;         // number of ensemble members for j-th model (i.e. nⱼ)

    /* some of the cmip6 models do not have the dissicnat variable
       for those model we only know the diffference, ΔCant - ΔCnat */
    int<lower=0> N2a;        // number of cmip6(w.o. dissicnat) models with nⱼ = 1
    real<lower=0> d2a[N2a];  // (ΔCant - ΔCnat) of j-th cmip6(w.o. dissicnat) model with nⱼ = 1
    
    int<lower=0> N2b;        // number of cmip6 models with nⱼ > 1
    vector[N2b] d2b;         // (ΔCant - ΔCnat) of j-th cmip6(w.o. dissicnat) model with nⱼ > 1
    vector[N2b] sig2b;        // std(ΔCant - ΔCnat) of j-th cmip6(w.o. dissicnat) model;
    vector[N2b] n2b;         // number of ensemble members for j-th model (i.e. nⱼ)

    /* some of the cmip6 models have the dissicnat variable
       for those models we can know ΔCant and ΔCnat separately
        all those models have more than one ensemble member, i.e. nⱼ > 1 */
    int<lower=0> N3;         // number of cmip6 models with the dissicnat variable
    vector[N3] d3;           // ΔCant of the j-th cmip6(w. dissicnat) model 
    vector[N3] sig3;         // std(ΔCant) of the j-th cmip6(w. dissicnat) model 
    vector[N3] n3;           // number of ensemble members for the j-th cmip6(w. dissicnat) model (i.e. nⱼ)

    int<lower=0> N4;         // number of cmip6 models with the dissicnat variable 
    vector[N4] d4;           // (ΔCnat) of the j-th cmip6(w. dissicnat) model 
    vector[N4] sig4;         // std(ΔCnat) of the j-th cmip6(w. dissicnat) model
    vector[N4] n4;           // number of ensemble members for the j-th cmip6(w. dissicnat) model (i.e. nⱼ)
    }
parameters {
    real<lower=0> sigXMa;  // standard deviation of the random ΔCant effects (i.e. of δaⱼ) that are specific to the j-th model (i.e. fixed for the j-th model ensemble )
    real<lower=0> sigXMn;  // standard deviation of the random ΔCnat effects (i.e. of δnⱼ) that are specific to the j-th model (i.e. fixed for the j-th model ensemble )
    real<lower=0> sigSMa;  // standard deviation of the random ΔCant effects (i.e. of ϵaⱼᵢ) that fluctuate for each ensemble member due to natural climate variability
    real<lower=0> sigSMn;  // standard deviation of the random ΔCnat effects (i.e. of ϵnⱼᵢ) that fluctuate for each ensemble member due to natural climate variability
    real DeltaCant;        // multimodel mean ΔCant (i.e. μΔCant) fixed for all models
    real DeltaCnat;        // multimodel mean ΔCnat (i.e. μΔCnat) fixed for all models
    real GRDeltaCant;       // GR19 ΔCant estimate
    vector[N1a] dX1a_a;    // δaⱼ : ΔCant-random-effect for the j-th dissicnat model with nⱼ = 1
    vector[N1a] dX1n_a;    // δnⱼ : ΔCnat-random-effect for the j-th dissicnat model with nⱼ =  1
    vector[N1b] dX1a_b;    // δaⱼ : ΔCant-random-effect fixed for all the ensemble members of the j-th cmip5 model that has nⱼ > 1
    vector[N1b] dX1n_b;    // δnⱼ : ΔCnat-random-effect fixed for all the ensemble members of the j-th cmip5 model that has nⱼ > 1
    vector[N2a] dX2a_a;    // δaⱼ : ΔCant-random-effect fixed for the j-th cmip6(w.o. dissicnat) model that has nⱼ = 1
    vector[N2a] dX2n_a;    // δnⱼ : ΔCnat-random-effect fixed for the j-th cmip6(w.o. dissicnat) model that has nⱼ = 1 
    vector[N2b] dX2a_b;    // δaⱼ : ΔCant-random-effect fixed for all the ensemble members of the j-th cmip6(w.o. dissicnat) model that has nⱼ > 1
    vector[N2b] dX2n_b;    // δnⱼ : ΔCnat-random-effect fixed for all the ensemble members of the j-th cmip6(w.o. dissicnat) model that has nⱼ > 1
    vector[N3] dX3a;       // δaⱼ : ΔCant-random-effect fixed for all the ensemble members of the j-th cmip6(w. dissicnat) model 
    vector[N4] dX4n;       // δnⱼ : ΔCnat-random-effect fixed for all the ensemble members of the j-th cmip6(w. dissicnat) model
    }
model {
    target += log(1/sigXMa);   // prior density for std(δaⱼ)
    target += log(1/sigXMn);   // prior density for std(δnⱼ)
    target += log(1/sigSMa);   // prior density for std(ϵaⱼᵢ)
    target += log(1/sigSMn);   // prior density for std(ϵnⱼᵢ)
//   target += inv_gamma_lpdf( sigXMa | 0.001, 0.001 );
//   target += inv_gamma_lpdf( sigXMn | 0.001, 0.001 );
//   target += inv_gamma_lpdf( sigSMa | 0.001, 0.001 );   // prior density for std(ϵaⱼᵢ)
//   target += inv_gamma_lpdf( sigSMn | 0.001, 0.001 );   // prior density for std(ϵnⱼᵢ)
    target += se_lpdf( sig1b | n1b,  sqrt( square( sigSMa ) + square( sigSMn ) ) );  // log likelihood for standard errors of cmip5 models with nⱼ > 1
    target += se_lpdf( sig2b | n2b,  sqrt( square( sigSMa ) + square( sigSMn ) ) );  // log likelihood for standard errors of cmip6(w.o. dissicnat) models with nⱼ > 1
    target += se_lpdf( sig3 | n3,   sigSMa ); // log likelihood for standard errors of ΔCant for cmip6(w. dissicnat) models
    target += se_lpdf( sig4 | n4,   sigSMn ); // log likelihood for standard errors of (ΔCnat) for cmip6(w. dissicnat) models
    target += normal_lpdf( dX1a_a | 0, sigXMa ); // log probability density for δaⱼ of cmip5 models with nⱼ = 1 
    target += normal_lpdf( dX1n_a | 0, sigXMn ); // log probability density for δnⱼ of cmip5 models with nⱼ = 1
    target += normal_lpdf( dX1a_b | 0, sigXMa ); // log probability density for δaⱼ of cmip5 models with nⱼ > 1
    target += normal_lpdf( dX1n_b | 0, sigXMn ); // log probability density for δnⱼ of cmip5 models with nⱼ > 1
    target += normal_lpdf( dX2a_a | 0, sigXMa ); // log probability density for δaⱼ of cmip6(w.o. dissicnat) models with nⱼ = 1
    target += normal_lpdf( dX2n_a | 0, sigXMn ); // log probability density for δnⱼ of cmip6(w.o. dissicnat) models with nⱼ = 1
    target += normal_lpdf( dX2a_b | 0, sigXMa ); // log probability denstiy for δaⱼ of cmip6(w.o. dissicnat) models with nⱼ > 1
    target += normal_lpdf( dX2n_b | 0, sigXMn ); // log probability density for δnⱼ of cmip6(w.o. dissicnat) models with nⱼ > 1
    target += normal_lpdf( dX3a| 0, sigXMa );    // log probability density for δaⱼ of cmip6(w. dissicnat) models
    target += normal_lpdf( dX4n| 0, sigXMn );    // log probability density for δnⱼ of cmip6(w. dissicnat) models
    target += normal_lpdf( GRDeltaCant | 33.0, 4.0); // log probability density for GR19ΔCant estimate
   
    // log likelihood for ΔCant-ΔCnat of cmip5 models with nⱼ = 1
    target += normal_lpdf( d1a | DeltaCant - DeltaCnat + dX1a_a - dX1n_a , sqrt( square( sigSMa ) + square( sigSMn ) ) );
    
    // log likelihood for ΔCant-ΔCnat of cmip5 models with nⱼ > 1
    target += normal_lpdf( d1b | DeltaCant - DeltaCnat + dX1a_b + dX1n_b , sqrt( square( sigSMa ) + square( sigSMn ) )./sqrt(n1b));
    
    // log likelihood for ΔCant-ΔCnat of cmip6(w.o. dissicnat) models with nⱼ = 1
    target += normal_lpdf( d2a | DeltaCant - DeltaCnat + dX2a_a + dX2n_a , sqrt( square( sigSMa ) + square( sigSMn ) ) );

    // log likelihood for ΔCant-ΔCnat of cmip6(w.o. dissicnat) models with nⱼ >1
    target += normal_lpdf( d2b | DeltaCant - DeltaCnat + dX2a_b + dX2n_b , sqrt( square( sigSMa ) + square( sigSMn ) )./sqrt(n2b) );
    
    // log likelihood for ΔCant of cmip6(w. dissicnat) models
    target += normal_lpdf( d3    | DeltaCant + dX3a , sigSMa./sqrt(n3) );
    
    // log likelihood for ΔCnat of cmip6(w. dissicnat) models
    target += normal_lpdf( d4 | DeltaCnat + dX4n , sigSMn./sqrt(n4) );

}