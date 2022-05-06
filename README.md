# Bayesian-Hierarchical-Model-for-the-CMIP5-6-anthropogenic-C-uptake

See 
"Evaluation ofthe CMIPmodels withtheInternational Ocean Model Benchmarktool: Rates of contemporary ocean carbon uptake linked with vertical temperature gradients and transport to the ocean interior" by Fu, Moore, Primeau, Collier, Oluwaseum, Ogunro, Hoffman, and Randerson"


![image](https://user-images.githubusercontent.com/7294603/165800464-96b582ba-39ba-45a1-976f-dd8d3656982d.png)
![bias_cdf](https://user-images.githubusercontent.com/7294603/165798190-e5fe2ed0-9ab3-45b7-a4ee-6d68db033869.png)

```
         Bayesian Hierarchical Model of the Carbon Storage in  CMIP5 and CMIP6  ESMs
      _________________________________________________________________________________

                                 ΔCant  = 27.789 ±  0.447
                                 ΔCnat  =  0.752 ±  0.221
                                σ(δant) =  1.829 ±  0.357
                                σ(δnat) =  0.278 ±  0.310
                                τ(ϵant) =  0.215 ±  0.022
                                τ(ϵnat) =  0.439 ±  0.030
                                 ΔCant(CMIP5)  = 27.887 ±  1.080
                                 ΔCant(CMIP6)  = 27.751 ±  1.960
                                 ΔCnat(CMIP5)  =  0.748 ±  0.342
                                 ΔCnat(CMIP6)  =  0.774 ±  0.428
                                 ΔCtot(CMIP5)  = 28.636 ±  1.184
                                 ΔCtot(CMIP6)  = 28.524 ±  2.054


                      Posterior ΔCant estimate for the cmip5 models with nⱼ = 1

      ΔCant - ΔCnat           ΔCant                ΔCnat       |       ΔCant                 ΔCnat
     -----data-----      ------data-----      -----data-----   |  ----BHM post---      ----BHM post---
 1   26.470 ±    NaN        NaN ±    NaN        NaN ±    NaN   |  27.252 ±  0.588       0.773 ±  0.379
 2   25.360 ±    NaN        NaN ±    NaN        NaN ±    NaN   |  26.302 ±  0.596       0.832 ±  0.389
 3   26.680 ±    NaN        NaN ±    NaN        NaN ±    NaN   |  27.465 ±  0.582       0.774 ±  0.386
 4   28.020 ±    NaN        NaN ±    NaN        NaN ±    NaN   |  28.709 ±  0.648       0.716 ±  0.438
 5   29.860 ±    NaN        NaN ±    NaN        NaN ±    NaN   |  30.300 ±  0.729       0.625 ±  0.573


                      Posterior ΔCant estimate for the cmip5 models with nⱼ > 1

      ΔCant - ΔCnat           ΔCant                ΔCnat       |       ΔCant                 ΔCnat
     -----data-----      ------data-----      -----data-----   |   ----BHM post---      ----BHM post---
 6   24.420 ±  0.250        NaN ±    NaN        NaN ±    NaN   |  25.347 ±  0.483       0.620 ±  0.575
 7   26.260 ±  0.350        NaN ±    NaN        NaN ±    NaN   |  27.077 ±  0.442       0.710 ±  0.442
 8   26.850 ±  0.180        NaN ±    NaN        NaN ±    NaN   |  27.594 ±  0.481       0.761 ±  0.419
 9   28.770 ±  0.410        NaN ±    NaN        NaN ±    NaN   |  29.425 ±  0.518       0.829 ±  0.412
10   29.010 ±  0.230        NaN ±    NaN        NaN ±    NaN   |  29.627 ±  0.577       0.842 ±  0.421
11   26.950 ±  0.260        NaN ±    NaN        NaN ±    NaN   |  27.720 ±  0.478       0.744 ±  0.406


              Posterior ΔCant estimate for the cmip6(w.o. dissicnat) models with nⱼ = 1

      ΔCant - ΔCnat           ΔCant                ΔCnat       |       ΔCant                 ΔCnat
     -----data-----      ------data-----      -----data-----   |  ----BHM post---      ----BHM post---
12   32.120 ±    NaN        NaN ±    NaN        NaN ±    NaN   |  32.302 ±  0.922       0.982 ±  0.647
13   29.080 ±    NaN        NaN ±    NaN        NaN ±    NaN   |  29.576 ±  0.661       0.844 ±  0.424


              Posterior ΔCant estimate for the cmip6(w.o. dissicnat) models with nⱼ > 1

      ΔCant - ΔCnat           ΔCant                ΔCnat       |       ΔCant                 ΔCnat
     -----data-----      ------data-----      -----data-----   |  ----BHM post---      ----BHM post---
14   27.370 ±  0.410        NaN ±    NaN        NaN ±    NaN   |  28.105 ±  0.418       0.763 ±  0.392
15   24.280 ±  0.510        NaN ±    NaN        NaN ±    NaN   |  25.175 ±  0.468       0.612 ±  0.610
16   25.610 ±  0.450        NaN ±    NaN        NaN ±    NaN   |  26.441 ±  0.407       0.681 ±  0.488
17   26.370 ±  0.620        NaN ±    NaN        NaN ±    NaN   |  27.157 ±  0.378       0.714 ±  0.448


              Posterior ΔCant estimate for the cmip6(w. dissicnat) models with nⱼ > 1

      ΔCant - ΔCnat           ΔCant                ΔCnat       |       ΔCant                 ΔCnat
     -----data-----      ------data-----      -----data-----   |  ----BHM post---      ----BHM post---
18   26.020 ±  0.390     26.640 ±  0.150      0.620 ±  0.360   |  26.643 ±  0.072       0.702 ±  0.152
19   25.770 ±  0.411     26.250 ±  0.180      0.480 ±  0.370   |  26.261 ±  0.124       0.671 ±  0.216
20   26.900 ±  0.599     27.970 ±  0.260      1.070 ±  0.540   |  27.969 ±  0.044       1.001 ±  0.099
21   27.160 ±  0.286     27.880 ±  0.120      0.720 ±  0.260   |  27.877 ±  0.123       0.768 ±  0.185
2022-05-05T17:04:57.069
```


