using StanSample, MCMCChains, KernelDensity
using Plots, StatsPlots, Printf, Dates
ProjDir = @__DIR__
cmipmod ="
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
"
cmipData = Dict( 
    "N1a"   => 5,
    "d1a"   => [ 26.47, 25.36, 26.68, 28.02, 29.86 ], # hist-pi (Pg C) cmip5
    "N1b"   => 6,
    "d1b"   => [ 24.42, 26.26, 26.85, 28.77, 29.01, 26.95 ], # hist-pi (Pg C) cmip5
    "sig1b" => [ 0.25,   0.35,  0.18,  0.41,  0.23,  0.26 ],   # cmip5 ensemble stds
    "n1b"   => [    5,      4,     3,     6,     3,     3 ],  # number of ensemble members
    "N2a"   => 2,
    "d2a"   => [32.12, 29.08 ], #hist -pi (Pg C) 6 models without dissicnat cmip6 w. no dissicnat
    "N2b"   => 4,
    "d2b"   => [27.37,  24.28, 25.61, 26.37 ], #hist -pi (Pg C) 6 models without dissicnat cmip6 w. no dissicnat
    "sig2b" => [ 0.41,   0.51,  0.45,  0.62 ], # cmip6 ensemble stds
    "n2b"   => [   11,     35,    10,    30 ], # number of ensemble members
    "N3"    => 4,
    "d3"    => [ 26.64, 26.25, 27.97, 27.88 ], #ΔCant (i.e. hist-dissicnat) 4 cmip6 models with dissicnat
    "sig3"  => [  0.15,  0.18,  0.26,  0.12 ], # ensemble stds
    "n3"    => [     9,     3,    25,     3 ], # number of ensemble members
    "N4"    => 4,
    "d4"    => [ 0.62, 0.48, 1.07, 0.72],  # ΔCnat for the 4 cmip6 models with dissicnat
    "sig4"  => [ 0.36, 0.37, 0.54, 0.26 ], # ensemble stds
    "n4"    => [    9,    3,   25,    3 ]
)

    # Keep tmpdir across multiple runs to prevent re-compilation
tmpdir = joinpath(ProjDir, "tmp")
isdir(tmpdir) &&  rm(tmpdir; recursive=true)

sm = SampleModel("cmip", cmipmod, [4]; tmpdir = tmpdir )
#sm.method = StanSample.Sample(15000, 15000, false, 1, StanSample.Adapt(true, 0.05, 0.8, 0.75, 10.0, 75, 50, 25), StanSample.Hmc(StanSample.Nuts(10), StanSample.diag_e(), 1.0, 1.0))
sm.method = StanSample.Sample(num_samples=200000)
rc = stan_sample(sm; data = cmipData)


if success(rc)
    println()
    read_summary(sm, true)
    chns = read_samples(sm)
end


m1 = [ cmipData["d1a"]; cmipData["d1b"]; cmipData["d2a"]; cmipData["d2b"]; cmipData["d3"]-cmipData["d4"] ];
m2 = [ fill(NaN,21-4,1); cmipData["d3"]];
m3 = [ fill(NaN,21-4,1); cmipData["d4"]];
m4 = [ mean(repeat(chns[:DeltaCant]',5)+chns[:dX1a_a],dims=2); mean(repeat(chns[:DeltaCant]',6)+chns[:dX1a_b],dims=2); mean(repeat(chns[:DeltaCant]',2)+chns[:dX2a_a],dims=2); mean(repeat(chns[:DeltaCant]',4)+chns[:dX2a_b],dims=2); mean(repeat(chns[:DeltaCant]',4)+chns[:dX3a],dims=2)]; 
m5 = [ mean(repeat(chns[:DeltaCnat]',5)+chns[:dX1n_a],dims=2); mean(repeat(chns[:DeltaCnat]',6)+chns[:dX1n_b],dims=2); mean(repeat(chns[:DeltaCnat]',2)+chns[:dX2n_a],dims=2); mean(repeat(chns[:DeltaCnat]',4)+chns[:dX2n_b],dims=2); mean(repeat(chns[:DeltaCnat]',4)+chns[:dX4n],dims=2)];
M = [m1 m2 m3 m4 m5];
s1 = [ fill(NaN,5,1); cmipData["sig1b"]; fill(NaN,2,1); cmipData["sig2b"]; fill(NaN,4,1)  ];
s2 = [ fill(NaN,21-4,1); cmipData["sig3"]];
s3 = [ fill(NaN,21-4,1); cmipData["sig4"]];
s4 = [ std(repeat(chns[:DeltaCant]',5)+chns[:dX1a_a],dims=2); std(repeat(chns[:DeltaCant]',6)+chns[:dX1a_b],dims=2); std(repeat(chns[:DeltaCant]',2)+chns[:dX2a_a],dims=2); std(repeat(chns[:DeltaCant]',4)+chns[:dX2a_b],dims=2); std(repeat(chns[:DeltaCant]',4)+chns[:dX3a],dims=2)]; 
s5 = [ std(repeat(chns[:DeltaCnat]',5)+chns[:dX1n_a],dims=2); std(repeat(chns[:DeltaCnat]',6)+chns[:dX1n_b],dims=2); std(repeat(chns[:DeltaCnat]',2)+chns[:dX2n_a],dims=2); std(repeat(chns[:DeltaCnat]',4)+chns[:dX2n_b],dims=2); std(repeat(chns[:DeltaCnat]',4)+chns[:dX4n],dims=2)];
S = [s1 s2 s3 s4 s5];

fid = open("results_"*string(now())*".txt","w");
@printf(fid,"         Bayesian Hierarchical Model of the Carbon Storage in  CMIP5 and CMIP6  ESMs\n");
@printf(fid,"      _________________________________________________________________________________\n\n");
@printf(fid,"                                 ΔCant  = %6.3f ± %6.3f  \n                                 ΔCnat  = %6.3f ± %6.3f \n", 
    mean(chns[:DeltaCant]), std(chns[:DeltaCant]), mean(chns[:DeltaCnat]), std(chns[:DeltaCnat]));
    @printf(fid,"                                σ(δant) = %6.3f ± %6.3f \n",          
    mean(chns[:sigXMa]), std(chns[:sigXMa]));
<<<<<<< HEAD
    @printf(fid,"                                σ(δnat) = %6.3f ± %6.3f \n",          
    mean(chns[:sigXMn]), std(chns[:sigXMn]));
    @printf(fid,"                                τ(ϵant) = %6.3f ± %6.3f \n",          
    mean(chns[:sigSMa]), std(chns[:sigSMa])); 
    @printf(fid,"                                τ(ϵnat) = %6.3f ± %6.3f \n",          
=======
    @printf(fid,"                                σ(δant) = %6.3f ± %6.3f \n",          
    mean(chns[:sigXMn]), std(chns[:sigXMn]));
    @printf(fid,"                                σ(δant) = %6.3f ± %6.3f \n",          
    mean(chns[:sigSMa]), std(chns[:sigSMa])); 
    @printf(fid,"                                σ(δant) = %6.3f ± %6.3f \n",          
>>>>>>> 68edbe80cb1cc1086358370e11cb889f05801b2c
    mean(chns[:sigSMn]), std(chns[:sigSMn]));

for i in 1:21
    if i == 1
        @printf(fid,"\n\n");
        @printf(fid, "                      Posterior ΔCant estimate for the cmip5 models with nⱼ = 1\n\n");
        @printf(fid,"      ΔCant - ΔCnat           ΔCant                ΔCnat       |       ΔCant                 ΔCnat     \n")
        @printf(fid,"     -----data-----      ------data-----      -----data-----   |  ----BHM post---      ----BHM post--- \n")
    elseif i == 6
        @printf(fid,"\n\n");
        @printf(fid, "                      Posterior ΔCant estimate for the cmip5 models with nⱼ > 1\n\n");
        @printf(fid,"      ΔCant - ΔCnat           ΔCant                ΔCnat       |       ΔCant                 ΔCnat     \n")
        @printf(fid,"     -----data-----      ------data-----      -----data-----   |   ----BHM post---      ----BHM post--- \n")
    elseif i == 12
        @printf(fid,"\n\n");
        @printf(fid, "              Posterior ΔCant estimate for the cmip6(w.o. dissicnat) models with nⱼ = 1\n\n");
        @printf(fid,"      ΔCant - ΔCnat           ΔCant                ΔCnat       |       ΔCant                 ΔCnat     \n")
        @printf(fid,"     -----data-----      ------data-----      -----data-----   |  ----BHM post---      ----BHM post--- \n")
    elseif i == 14
        @printf(fid,"\n\n");
        @printf(fid, "              Posterior ΔCant estimate for the cmip6(w.o. dissicnat) models with nⱼ > 1\n\n");
        @printf(fid,"      ΔCant - ΔCnat           ΔCant                ΔCnat       |       ΔCant                 ΔCnat     \n")
        @printf(fid,"     -----data-----      ------data-----      -----data-----   |  ----BHM post---      ----BHM post--- \n")
    elseif i == 18
        @printf(fid,"\n\n");
        @printf(fid, "              Posterior ΔCant estimate for the cmip6(w. dissicnat) models with nⱼ > 1\n\n");
        @printf(fid,"      ΔCant - ΔCnat           ΔCant                ΔCnat       |       ΔCant                 ΔCnat     \n")
        @printf(fid,"     -----data-----      ------data-----      -----data-----   |  ----BHM post---      ----BHM post--- \n")
    end 
    @printf(fid,"%2i   %6.3f ± %6.3f     %6.3f ± %6.3f     %6.3f ± %6.3f   |  %6.3f ± %6.3f      %6.3f ± %6.3f     \n",
    i,M[i,1],S[i,1],M[i,2],S[i,2],M[i,3],S[i,3],M[i,4],S[i,4],M[i,5],S[i,5]);
end
@printf(fid,"%s",now())
close(fid)

#boxplot(repeat([1,2,3,4],inner=60000),[chns[:sigXMa];chns[:sigSMa];chns[:sigXMn];chns[:sigSMn]])
U = kde(chns[:DeltaCant]-chns[:GRDeltaCant])
p1 = plot();
p1 = plot(U.x,U.density);
function cumtrapz(x,y)
    out = similar(y)
    out[1] = 0
    for i in 2:length(x)
        out[i] = out[i-1] + 0.5*(x[i]-x[i-1])*(y[i]+y[i-1])
    end
    return out
end
cdf = cumtrapz(U.x,U.density);
a,iz = findmin( abs.(U.x) ); 
p2 = plot();
p2 = plot!(U.x,cdf,xlabel = "ΔCant bias (CMIP - GR2019)   [Pg C]",ylabel = "Cumulative Probability",legend = false, lw = 2, color = :black, yticks = collect(0:0.1:1));

p2 = plot!([0,0,0],[0,0.5,cdf[iz]], lw = 1, linestyle = :dash, color = :black);
p2 = plot!([-30,-10,0], [ cdf[iz], cdf[iz], cdf[iz] ], lw = 1, linestyle = :dash, color = :black );

display(p1)
display(p2)

