eDNA study design simulations
================
Roy Martin
March 07, 2025

This document describes and provides code for a probabilistic (“prior
predictive”) simulation of eDNA detection data for a hypothetical target
organism across a user-defined distribution of sites, samples within
sites, and replicate qPCR observations of those samples. In addition,
sites were distributed among two land use types, which was parameterized
to have an effect on occupancy rates.

``` r
# Lets now set up our R environment to run the simulation.
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(RColorBrewer)
library(truncnorm)
library(stringr)
library(tidyverse)
library(rstan)
library(loo)
library(bayesplot)
library(tidybayes)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)
options(loo.cores = 4)

options(max.print = 9999)
```

# 1 Model description

The model is a three-level occupancy model based on Nichols et al.
([2008](#ref-Nichols_etal_2008)) and Schmidt et al.
([2013](#ref-Schmidt_etal_2013)). The basic model below consists of a
sequence of three coupled Bernoulli trials for describing the array of
nested data,
![y\_{i,j,k}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Cj%2Ck%7D "y_{i,j,k}"),
where we have an observed detection
(![y\_{ijk} = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bijk%7D%20%3D%201 "y_{ijk} = 1"))
or not
(![y\_{ijk} = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bijk%7D%20%3D%200 "y_{ijk} = 0"))
for each PCR replicate
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k"),
nested in a water sample
![j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j "j"),
nested in a sample visit day
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")
at a site
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i").
If we plan to completely cross site and day, we would have
![S \times T = 30 \times 17 = 510](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S%20%5Ctimes%20T%20%3D%2030%20%5Ctimes%2017%20%3D%20510 "S \times T = 30 \times 17 = 510")
total sampling occasions.

![\textbf{Level-3: latent state of occupancy}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BLevel-3%3A%20latent%20state%20of%20occupancy%7D "\textbf{Level-3: latent state of occupancy}")

![z](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z "z")
= latent state (0 or 1) of occupancy
![\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi "\psi")
= probability of occupancy

![z\_{i} \sim Bernoulli(\psi_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bi%7D%20%5Csim%20Bernoulli%28%5Cpsi_i%29 "z_{i} \sim Bernoulli(\psi_i)")

![\textbf{Level-2: availability parameter}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BLevel-2%3A%20availability%20parameter%7D "\textbf{Level-2: availability parameter}")

![a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a "a")
= latent state (0 or 1) of availability (in water sample) for detection
via qPCR
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")
= probability of capturing DNA in water sample and it being available
for observation, conditional on occupancy.

![a\_{ij} \| z_i \sim Bernoulli( z_i\theta\_{ij} )](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_%7Bij%7D%20%7C%20z_i%20%5Csim%20Bernoulli%28%20z_i%5Ctheta_%7Bij%7D%20%29 "a_{ij} | z_i \sim Bernoulli( z_i\theta_{ij} )")

![\textbf{Level-1: observations}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BLevel-1%3A%20observations%7D "\textbf{Level-1: observations}")

![y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y "y")
= observed detection (0 or 1) via qPCR
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
= probability of detecting via qPCR, conditional on it being available
in water sample.

![y\_{ijk} \| a\_{ij} \sim Bernoulli( a\_{ij} p\_{ijk} )](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bijk%7D%20%7C%20a_%7Bij%7D%20%5Csim%20Bernoulli%28%20a_%7Bij%7D%20p_%7Bijk%7D%20%29 "y_{ijk} | a_{ij} \sim Bernoulli( a_{ij} p_{ijk} )")

\# Prior predictive simulation

## 1.1 Dimensions

Set the dimensions of the hypothetical data based on the number of
sites, water samples, and qPCR replicates.

``` r
nag <- 15 # number of ag sites in study
nurban <- 15 # number of urban sites in study
nsite <- nag + nurban # number of sites in study
nsamp <- 5 # number of replicate water samples per site
nrep <- 4 # number of replicate qPCR samples per water sample (per site and date)
nsim <- 1e3 # number of simulated draws from joint prior predictive distribution
```

The number of sites is 30. The total number of water samples is 150. The
total number of qPCR runs is 600.

Next create containers at each level of the model to hold the simulated
values.

## 1.2 Priors

Set priors for the parameters to be estimated from data.

### 1.2.1 ![\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi "\psi")

The linear predictor for the occupancy parameter includes an intercept
term, a land use effect, and a random site effect:

![z\_{i} \sim Bernoulli(\psi_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bi%7D%20%5Csim%20Bernoulli%28%5Cpsi_i%29 "z_{i} \sim Bernoulli(\psi_i)")

![logit( \psi_i ) = \alpha\_\psi + \beta\*X_i + \gamma\_{\psi_i}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;logit%28%20%5Cpsi_i%20%29%20%3D%20%5Calpha_%5Cpsi%20%2B%20%5Cbeta%2AX_i%20%2B%20%5Cgamma_%7B%5Cpsi_i%7D "logit( \psi_i ) = \alpha_\psi + \beta*X_i + \gamma_{\psi_i}")

#### 1.2.1.1 Intercept

``` r
# prior for log-odds-scale intercept parameter for psi
loc_a_psi = 0
scale_a_psi = 0.5
a_psi <- rnorm(nsim, loc_a_psi, scale_a_psi)
```

The intercept in this model can be interpreted as the baseline
probability of occupancy for a hypothetical “ag” site. It is centered
over 0.5 and varies around that baseline by according to the scale
parameter representing *a priori* uncertainty in that parameter.
<img src="multilevel_occupancy_design_scenario_files/figure-gfm/prior_psi_prob-1.png" style="display: block; margin: auto;" />

#### 1.2.1.2 Land use effect

``` r
# prior for log-odds-scale intercept parameter for psi
loc_beta_psi = 1
scale_beta_psi = 0.2
beta_psi <- rnorm(nsim, loc_beta_psi, scale_beta_psi)
```

The land use effect is centered over 1 on the logit scale, but in the
context of the baseline (intercept above) and on the probability scale,
it shifts the baseline probability centered over 0.5 for “ag” sites to a
probability closer to 0.75 for “urban” sites.
<img src="multilevel_occupancy_design_scenario_files/figure-gfm/prior_psi_effect_prob-1.png" style="display: block; margin: auto;" />

#### 1.2.1.3 ![\sigma\_{\psi}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cpsi%7D "\sigma_{\psi}")

``` r
# prior for scale parameter of varying site effects
loc_sigma_gamma_psi = 0
scale_sigma_gamma_psi = 0.3
sigma_gamma_psi <- rtruncnorm(nsim, a = 0, mean = loc_sigma_gamma_psi, sd = scale_sigma_gamma_psi)
```

The sigma parameter here defines scale of (“random”) site-to-site
effects on
![\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi "\psi")
not captured by the “fixed” intercept and land use effects. Effectively,
![\gamma\_{\psi}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cgamma_%7B%5Cpsi%7D "\gamma_{\psi}")
defines unmeasured process variation between sites and
![\sigma\_{\psi}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cpsi%7D "\sigma_{\psi}")
defines the scale of that variation.

<img src="multilevel_occupancy_design_scenario_files/figure-gfm/prior_sigma_gamma_psi_plot-1.png" style="display: block; margin: auto;" />

### 1.2.2 ![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")

![a\_{ij} \| z_i \sim Bernoulli( z_i\theta\_{ij} )](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_%7Bij%7D%20%7C%20z_i%20%5Csim%20Bernoulli%28%20z_i%5Ctheta_%7Bij%7D%20%29 "a_{ij} | z_i \sim Bernoulli( z_i\theta_{ij} )")

![logit( \theta\_{ij} ) = \alpha\_\theta + \gamma\_{\theta_i}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;logit%28%20%5Ctheta_%7Bij%7D%20%29%20%3D%20%5Calpha_%5Ctheta%20%2B%20%5Cgamma_%7B%5Ctheta_i%7D "logit( \theta_{ij} ) = \alpha_\theta + \gamma_{\theta_i}")

#### 1.2.2.1 Intercept

``` r
# prior for log-odds-scale intercept parameter for psi
loc_a_theta = 1.5
scale_a_theta = 0.25
a_theta <- rnorm(nsim, loc_a_theta, scale_a_theta)
```

The probability of obtaining the target eDNA, conditional on it being in
the site, was centered over about 0.8.
<img src="multilevel_occupancy_design_scenario_files/figure-gfm/prior_theta_prob-1.png" style="display: block; margin: auto;" />

#### 1.2.2.2 ![\sigma\_{\theta}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Ctheta%7D "\sigma_{\theta}")

``` r
# prior for scale parameter of varying site effects
loc_sigma_gamma_theta = 0
scale_sigma_gamma_theta = 0.1
sigma_gamma_theta <- rtruncnorm(nsim, a = 0, mean = loc_sigma_gamma_theta, sd = scale_sigma_gamma_theta)
```

Again,
![\sigma\_{\theta}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Ctheta%7D "\sigma_{\theta}")
defines the scale of unmeasured process variation around the fixed
intercept.
<img src="multilevel_occupancy_design_scenario_files/figure-gfm/prior_sigma_gamma_theta_plot-1.png" style="display: block; margin: auto;" />

### 1.2.3 p

![logit( p\_{ijk} ) = \alpha_p + \gamma\_{p_i}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;logit%28%20p_%7Bijk%7D%20%29%20%3D%20%5Calpha_p%20%2B%20%5Cgamma_%7Bp_i%7D "logit( p_{ijk} ) = \alpha_p + \gamma_{p_i}")

#### 1.2.3.1 Intercept

``` r
# prior for log-odds-scale intercept parameter for p
loc_a_p = 1.25
scale_a_p = 0.25
a_p <- rnorm(nsim, loc_a_p, scale_a_p)
```

The probability of detecting eDNA conditional on it being available in
the water sample is also centered over about 0.8.
<img src="multilevel_occupancy_design_scenario_files/figure-gfm/prior_p_prob-1.png" style="display: block; margin: auto;" />

#### 1.2.3.2 Sigma site

``` r
# prior for scale parameter of varying site effects
loc_sigma_gamma_p = 0
scale_sigma_gamma_p = 0.05
sigma_gamma_p <- rtruncnorm(nsim, a = 0, mean = loc_sigma_gamma_p, sd = scale_sigma_gamma_p)
```

The scale of unmeasured process variation in
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
is set pretty low.
<img src="multilevel_occupancy_design_scenario_files/figure-gfm/prior_sigma_gamma_p_plot-1.png" style="display: block; margin: auto;" />

## 1.3 Simulate data

Simulate draws from the prior predictive distribution based on our model
and priors defined above.

``` r
# Dummy variable for land use
X_land <- c(rep(0, nag), rep(1, nurban))
# Simulate varying effects
gamma_psi <- array(NA, dim = c(nsite, nsim))
gamma_theta <- array(NA, dim = c(nsite, nsim))
gamma_p <- array(NA, dim = c(nsite, nsim))

for (i in 1:nsite){
  gamma_psi[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_psi)
  gamma_theta[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_theta)
  gamma_p[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_p)
}

# Run the prior predictive simulation through likelihood
for(i in 1:nsite){
    psi_i[i, ] <- plogis(a_psi + beta_psi * X_land[i] + gamma_psi[i, ]) # back-transform
    z_i[i, ] <- rbinom(nsim, 1, psi_i[i, ])
  
  for (j in 1:nsamp) {
    theta_ij[i, j, ] <- plogis( a_theta + gamma_theta[i, ]) 
    a_ij[i, j, ] <- rbinom(nsim, 1, z_i[i, ] * theta_ij[i, j, ])
   
    for (k in 1:nrep) {
      p_ijk[i, j, k, ] <- plogis(a_p + gamma_p[i, ])
      y_ijk[i, j, k, ] <- rbinom(nsim, 1, a_ij[ i, j, ] * p_ijk[ i, j, k, ])
      }
    }
  }
```

### 1.3.1 Occupancy: Ag sites

using the prior predictive simulation above, we can summarize some of
the implications of our beliefs in terms of practical conditions in this
hypothetical system. For example, give our chosen model and our priors
beliefs about its parameters and our uncertainty in their true values,
the simulation indicats that we believe that anywhere between 1 and 15
of the 15 “Ag” sites are occupied, with the most likely value being
about 8.

    ## Warning: `qplot()` was deprecated in ggplot2 3.4.0.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

<img src="multilevel_occupancy_design_scenario_files/figure-gfm/distribution_occupancy_ag_sites-1.png" style="display: block; margin: auto;" />

### 1.3.2 Occupancy: Urban sites

Our simulation also indicates that we believe anywhere between 1 and 15
“urban” sites are occupied, with the most likely value being 11 and most
of probability is for a larger number of occupied sites relative to the
“ag” sites.
<img src="multilevel_occupancy_design_scenario_files/figure-gfm/distribution_occupancy_urban_sites-1.png" style="display: block; margin: auto;" />

### 1.3.3 Frequency of detections

We can look at similar implications for the probability of detections at
the qPCR scale, etc.
![](multilevel_occupancy_design_scenario_files/figure-gfm/distribution_pcr_detects-1.png)<!-- -->

### 1.3.4 Replicate datasets

We can also use the simulation to look at the hypothetical world of
potential datasets that could be collected given our model and priors.
Below we show 20 such possible datasets. And the exact parameters for
![\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi "\psi"),
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta"),
and
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
that they were generated from. Note that even for those fixed parameter
values, the dataset drawn could still look different than what was drawn
here.

<img src="multilevel_occupancy_design_scenario_files/figure-gfm/frequency_plot-1.png" style="display: block; margin: auto;" />

# 2 Design scenarios

## 2.1 30 sites

In this first scenario, we take the same priors from above, but we fix
the true land use effect to
![\beta\_{psi} = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_%7Bpsi%7D%20%3D%201 "\beta_{psi} = 1").
Then we manipulate the number of sites, samples, qPCR reps and see how
well our Stan model recovers the true parameter estimate and its
uncertainty. Below, we fix
![\beta\_{psi}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_%7Bpsi%7D "\beta_{psi}")
to be exactly 1 with no uncertainty and re-simulate a new batch of
hypothetical datasets. Then we randomly select a single simulated
dataset and use a
![\textbf{Stan}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BStan%7D "\textbf{Stan}")
model (below) to estimate the parameters of the model conditional on the
(hypothetical) data. To get a more rigorous summary of our design and
how things like sample sizes and prior uncertainty affect our ability to
detect the effect, we would want iterate this process thousands of times
and then summarize things like: how many times did our parameter
estimate for
![\beta\_{psi}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_%7Bpsi%7D "\beta_{psi}")
include the true value (1)? How many times did our estimate have the
right sign (positive vs. negative effect)? How many times or what
proportion of our estimates did not include zero? etc.

``` r
nag <- 15 # number of ag sites in study
nurban <- 15 # number of urban sites in study
nsite <- nag + nurban # number of sites in study
nsamp <- 5 # number of replicate water samples per site
nrep <- 4 # number of replicate qPCR samples per water sample (per site and date)
nsim <- 1e3 # number of simulated draws from joint prior predictive distribution
```

``` r
nsite
```

    ## [1] 30

``` r
nsite * nsamp
```

    ## [1] 150

The total number of qPCR runs is equat to the number of sites, times the
number of water samples, times the number of qPCR reps per water sample.
This number is shown below and perhaps gives a more practical sense of
the effort required to try to estimate this particular hypothetical
effect.

``` r
nsite * nsamp * nrep
```

    ## [1] 600

Fixed land use effect

``` r
# fixed slope effect of land use
beta_psi <- 1
```

### 2.1.1 Simulate draws

``` r
# Dummy variable for land use
X_land <- c(rep(0, nag), rep(1, nurban))
# Simulate varying effects
gamma_psi <- array(NA, dim = c(nsite, nsim))
gamma_theta <- array(NA, dim = c(nsite, nsim))
gamma_p <- array(NA, dim = c(nsite, nsim))

for (i in 1:nsite){
  gamma_psi[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_psi)
  gamma_theta[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_theta)
  gamma_p[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_p)
}

# Run the prior predictive simulation through likelihood
for(i in 1:nsite){
    psi_i[i, ] <- plogis(a_psi + beta_psi * X_land[i] + gamma_psi[i, ]) # back-transform
    z_i[i, ] <- rbinom(nsim, 1, psi_i[i, ])
  
  for (j in 1:nsamp) {
    theta_ij[i, j, ] <- plogis( a_theta + gamma_theta[i, ]) 
    a_ij[i, j, ] <- rbinom( nsim, 1, z_i[i, ] * theta_ij[i, j, ])
   
    for (k in 1:nrep) {
      p_ijk[i, j, k, ] <- plogis( a_p + gamma_p[i, ])
      y_ijk[i, j, k, ] <- rbinom( nsim, 1, a_ij[ i, j, ] * p_ijk[ i, j, k, ] )
      }
    }
  }
```

### 2.1.2 Pick a random dataset

Below we pick a single random dataset from the world of possible
datasets given our assumptions outlined above. The code below picks the
dataset and then we also print out the true parameter values used to
generate this specific datset. After fitting the model with this dataset
below, we can use these values to get a sense of how well the model
recovered the parameters.

``` r
sel <- sample(1:1000, 1, replace = FALSE)
fake_data <- y_ijk[, , , sel] # draw the fake dataset from 1000 simulations

data.frame(
  a_psi = a_psi[sel],
  beta_psi = beta_psi,
  sigma_gamma_psi = sigma_gamma_psi[[sel]],
  a_theta = a_theta[sel],
  sigma_gamma_theta = sigma_gamma_theta[[sel]],
  a_p = a_p[sel],
  sigma_gamma_p = sigma_gamma_p[[sel]]
) %>% t()
```

    ##                         [,1]
    ## a_psi             0.37145515
    ## beta_psi          1.00000000
    ## sigma_gamma_psi   0.21854123
    ## a_theta           1.56150788
    ## sigma_gamma_theta 0.05894662
    ## a_p               1.15557712
    ## sigma_gamma_p     0.04289127

### 2.1.3 A model in Stan

``` stan
data {
  int<lower = 1> nsite;
  int<lower = 1> nsamp;
  int<lower = 1> nrep;
  int<lower = 0, upper = 1> y[nsite, nsamp, nrep];
  int<lower = 0, upper =1> X_Land[nsite];
  int<lower=1> site[nsite];
  int<lower = 1> n_possible;
  matrix<lower = 0, upper = 1>[n_possible, nsamp] alpha_potential;
}

transformed data {
  int<lower = 0, upper = 1> known_present[nsite];
  int<lower = 0, upper = 1> known_available[nsite, nsamp];
  for (i in 1:nsite) {
    known_present[i] = 0;
    for (j in 1:nsamp) {
      known_available[i, j] = 0;
      for (k in 1:nrep) {
        if (y[i, j, k] == 1) {
         known_present[i] = 1;
         known_available[i, j] = 1;
        }
      }
    }
  }
}

parameters {
  real a_psi ;
  real a_theta;
  real a_p;
  real beta; // beta psi
  real<lower = 0> sigma_gamma_psi;
  real<lower = 0> sigma_gamma_theta;
  real<lower = 0> sigma_gamma_p;
  vector[nsite] z_psi;
  vector[nsite] z_theta;
  vector[nsite] z_p;
}

transformed parameters {
  vector[nsite] log_lik;
  real<lower = 0, upper = 1> psi[nsite];
  real<lower = 0, upper = 1> theta[nsite, nsamp];
  real<lower = 0, upper = 1> p[nsite, nsamp, nrep];
  vector[nsite] gamma_psi = z_psi * sigma_gamma_psi;
  vector[nsite] gamma_theta = z_theta * sigma_gamma_theta;
  vector[nsite] gamma_pdet = z_p * sigma_gamma_p;

  // linear predictor
  for(i in 1:nsite){
    psi[i] = inv_logit(a_psi + beta * X_Land[i] + gamma_psi[site[i]]); // linear predictor (logit scale) for psi
    for (j in 1:nsamp){
      theta[i, j] = inv_logit(a_theta + gamma_theta[site[i]]); // linear predictor (logit scale) for theta
      for(k in 1:nrep){
        p[i, j, k] = inv_logit(a_p + gamma_pdet[site[i]]); // linear predictor (logit scale) for p
        } // k
      } // j
    } // i
  
  {
    vector[nsamp] tmp_lp;
    matrix[n_possible, nsamp] tmp_poss;
    vector[n_possible + 1] sum_poss;
    
    for (i in 1:nsite) {
      if (known_present[i]) {
        for (j in 1:nsamp) {
           if (known_available[i, j]) {
             // present in site and available for water sample
             tmp_lp[j] = log(theta[i, j]) + bernoulli_lpmf(y[i, j, ] | p[i, j, ]);
          
             } else {
               // present, possibly unavailable for water sample
               tmp_lp[j] = log_sum_exp(
                 log(theta[i, j]) + bernoulli_lpmf(y[i, j, ] | p[i, j, ]), 
                 log1m(theta[i, j])
                 );
               }
        } // j( 1 )
        log_lik[i] = log(psi[i]) + sum(tmp_lp);
      } else {
        // could be present or absent (was never detected)
        // and there are 2^ntime possible combinations
        // of alpha_{i, j} that are relevant if z_i = 1
        for (jj in 1:n_possible) {
          for (j in 1:nsamp) {
            if (alpha_potential[jj, j] == 0) {
              // not available
              tmp_poss[jj, j] = log1m(theta[i, j]);
            } else {
              // available but not detected
              tmp_poss[jj, j] = log(theta[i , j]) + bernoulli_lpmf(y[i, j, ] | p[i, j, ]);
            }
          }
          sum_poss[jj] = log(psi[i]) + sum(tmp_poss[jj, ]);
        } // j( 2 )
        sum_poss[n_possible + 1] = log1m(psi[i]);
        log_lik[i] = log_sum_exp(sum_poss);
      }
    } // i
  }
}

model {
  // priors
  target += normal_lpdf(a_psi | 0, 2);
  target += normal_lpdf(a_theta | 0, 2);
  target += normal_lpdf(a_p | 0, 2);
  target += normal_lpdf(beta | 0, 2);
  target += normal_lpdf(sigma_gamma_psi | 0, 1);
  target += normal_lpdf(sigma_gamma_theta | 0, 1);
  target += normal_lpdf(sigma_gamma_p | 0, 1);
  target += normal_lpdf(z_psi | 0, 1);
  target += normal_lpdf(z_theta | 0, 1);
  target += normal_lpdf(z_p | 0, 1);
  
  // add log-likelihood
  target += sum(log_lik);
}
```

``` r
# potential combinations of alpha that can lead to all-zero capture history
alpha_potential <- expand.grid(rep(list(c(0, 1)), nsamp))

stan_d <- list(nsite = nsite, 
               nsamp = nsamp, 
               nrep = nrep, 
               X_Land = X_land,
               site = seq(1, nsite, 1),
               y = fake_data, 
               n_possible = 2 ^ nsamp, 
               alpha_potential = alpha_potential)
```

### 2.1.4 Fit model

``` r
fit_mod1 <- sampling(
  object = mocc,
  data = stan_d,
  chains = 4,
  iter = 2000,
  cores = 4,
  thin = 1
  )
```

### 2.1.5 Print posterior summary

``` r
print(fit_mod1, 
       pars = c("a_psi", "a_theta", "a_p", "beta"), 
       digits_summary = 3)
```

    ## Inference for Stan model: anon_model.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##          mean se_mean    sd   2.5%   25%   50%   75% 97.5% n_eff Rhat
    ## a_psi   0.979   0.011 0.635 -0.169 0.528 0.947 1.378 2.349  3423    1
    ## a_theta 1.803   0.006 0.330  1.207 1.583 1.778 2.002 2.516  3441    1
    ## a_p     1.246   0.002 0.136  0.985 1.154 1.243 1.334 1.527  5620    1
    ## beta    1.852   0.017 1.079 -0.068 1.127 1.790 2.507 4.224  3871    1
    ## 
    ## Samples were drawn using NUTS(diag_e) at Fri Mar  7 13:07:10 2025.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

![(0.996/0.5) ^ 2 = 4.0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%280.996%2F0.5%29%20%5E%202%20%3D%204.0 "(0.996/0.5) ^ 2 = 4.0")
![30 \times 4.0 = 120](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;30%20%5Ctimes%204.0%20%3D%20120 "30 \times 4.0 = 120")
Need about 120 sites.

``` r
#first extract posteriors
posteriors_fit <- rstan::extract(fit_mod1)

#plot
qplot(x = posteriors_fit$psi[ ,1], geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = plogis(a_psi[sel]), color = 'red', size = 2) +
  xlab(expression(psi))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_2-1.png)<!-- -->

``` r
qplot(x = posteriors_fit$theta[ ,1 ,1], geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = plogis(a_theta[sel]), color = 'red', size = 2) +
  xlab(expression(theta))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_2-2.png)<!-- -->

``` r
qplot(x = posteriors_fit$p[ , 1, 1 ,1], geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = plogis(a_p[sel]), color = 'red', size = 2) +
  xlab(expression(p))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_2-3.png)<!-- -->

``` r
qplot(x = posteriors_fit$beta, geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = beta_psi, color = 'red', size = 2) +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed", size = 2) +
  xlab(expression(Beta))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_2-4.png)<!-- -->

## 2.2 200 sites

Increase number of sites from
![15 \times 2 = 30](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;15%20%5Ctimes%202%20%3D%2030 "15 \times 2 = 30")
to
![100 \times 2 = 200](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;100%20%5Ctimes%202%20%3D%20200 "100 \times 2 = 200")

``` r
nag <- 100 # number of ag sites in study
nurban <- 100 # number of urban sites in study
nsite <- nag + nurban # number of sites in study
nsamp <- 5 # number of replicate water samples per site
nrep <- 4 # number of replicate qPCR samples per water sample (per site and date)
nsim <- 1e3 # number of simulated draws from joint prior predictive distribution
```

``` r
nsite
```

    ## [1] 200

``` r
nsite * nsamp
```

    ## [1] 1000

``` r
nsite * nsamp * nrep
```

    ## [1] 4000

### 2.2.1 Simulate draws

``` r
# Run the prior predictive simulation through likelihood
# Dummy variable for land use
X_land <- c(rep(0, nag), rep(1, nurban))
# Simulate varying effects
gamma_psi <- array(NA, dim = c(nsite, nsim))
gamma_theta <- array(NA, dim = c(nsite, nsim))
gamma_p <- array(NA, dim = c(nsite, nsim))

for (i in 1:nsite){
  gamma_psi[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_psi)
  gamma_theta[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_theta)
  gamma_p[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_p)
}

for(i in 1:nsite){
    psi_i[i, ] <- plogis(a_psi + beta_psi * X_land[i] + gamma_psi[i, ]) # back-transform
    z_i[i, ] <- rbinom(nsim, 1, psi_i[i, ])
  
  for (j in 1:nsamp) {
    theta_ij[i, j, ] <- plogis(a_theta + gamma_theta[i, ]) 
    a_ij[i, j, ] <- rbinom(nsim, 1, z_i[i, ] * theta_ij[i, j, ])
   
    for (k in 1:nrep) {
      p_ijk[i, j, k, ] <- plogis(a_p + gamma_p[i, ])
      y_ijk[i, j, k, ] <- rbinom(nsim, 1, a_ij[i, j, ] * p_ijk[i, j, k, ])
      }
    }
  }
```

### 2.2.2 Pick draws at random to analyze

``` r
sel <- sample(1:1000, 1, replace = FALSE)
fake_data <- y_ijk[, , , sel] # draw the fake dataset from 1000 simulations
data.frame(
  a_psi = a_psi[sel],
  beta_psi = beta_psi,
  sigma_gamma_psi = sigma_gamma_psi[[sel]],
  a_theta = a_theta[sel],
  sigma_gamma_theta = sigma_gamma_theta[[sel]],
  a_p = a_p[sel],
  sigma_gamma_p = sigma_gamma_p[[sel]]
) %>% t()
```

    ##                          [,1]
    ## a_psi             0.082110173
    ## beta_psi          1.000000000
    ## sigma_gamma_psi   0.061594874
    ## a_theta           1.297113667
    ## sigma_gamma_theta 0.025462763
    ## a_p               1.316779560
    ## sigma_gamma_p     0.006781541

``` r
# potential combinations of alpha that can lead to all-zero capture history
alpha_potential <- expand.grid(rep(list(c(0, 1)), nsamp))

stan_d <- list(nsite = nsite, 
               nsamp = nsamp, 
               nrep = nrep, 
               X_Land = X_land,
               site = seq(1, nsite, 1),
               y = fake_data, 
               n_possible = 2 ^ nsamp, 
               alpha_potential = alpha_potential)
```

### 2.2.3 Fit Stan model

``` r
fit_mod2 <- sampling(
  object = mocc,
  data = stan_d,
  chains = 4,
  iter = 2000,
  cores = 4,
  thin = 1
  )
```

    ## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
    ## Running the chains for more iterations may help. See
    ## https://mc-stan.org/misc/warnings.html#bulk-ess

### 2.2.4 Print posterior summary

``` r
print(fit_mod2, 
      pars = c("a_psi", "a_theta", "a_p", "beta"), 
      digits_summary = 3)
```

    ## Inference for Stan model: anon_model.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##          mean se_mean    sd   2.5%    25%   50%   75% 97.5% n_eff  Rhat
    ## a_psi   0.072   0.004 0.252 -0.416 -0.093 0.070 0.232 0.587  4302 0.999
    ## a_theta 1.478   0.002 0.131  1.233  1.387 1.473 1.563 1.750  3118 1.000
    ## a_p     1.308   0.001 0.057  1.199  1.270 1.308 1.345 1.421  5518 0.999
    ## beta    0.976   0.011 0.413  0.247  0.694 0.939 1.225 1.886  1491 1.002
    ## 
    ## Samples were drawn using NUTS(diag_e) at Fri Mar  7 13:19:49 2025.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

``` r
#first extract posteriors
posteriors_fit_2 <- rstan::extract(fit_mod2)

#plot
qplot(x = posteriors_fit_2$psi[ ,1], geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = plogis(a_psi[sel]), color = 'red', size = 2) +
  xlab(expression(psi))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_3-1.png)<!-- -->

``` r
qplot(x = posteriors_fit_2$theta[ ,1 ,1], geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = plogis(a_theta[sel]), color = 'red', size = 2) +
  xlab(expression(theta))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_3-2.png)<!-- -->

``` r
qplot(x = posteriors_fit_2$p[ , 1, 1 ,1], geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = plogis(a_p[sel]), color = 'red', size = 2) +
  xlab(expression(p))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_3-3.png)<!-- -->

``` r
qplot(x = posteriors_fit_2$beta, geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = beta_psi, color = 'red', size = 2) +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed", size = 2) +
  xlab(expression(Beta))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_3-4.png)<!-- -->

## 2.3 Same number of sites, fewer replicates

Same number of sites at
![100 \times 2 = 200](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;100%20%5Ctimes%202%20%3D%20200 "100 \times 2 = 200"),
but reduce number of water samples and pcr replicates.

``` r
nag <- 100 # number of ag sites in study
nurban <- 100 # number of urban sites in study
nsite <- nag + nurban # number of sites in study
nsamp <- 2 # number of replicate water samples per site
nrep <- 2 # number of replicate qPCR samples per water sample (per site and date)
nsim <- 1e3 # number of simulated draws from joint prior predictive distribution
```

``` r
nsite
```

    ## [1] 200

``` r
nsite * nsamp
```

    ## [1] 400

``` r
nsite * nsamp * nrep
```

    ## [1] 800

### 2.3.1 Simulate draws

``` r
# Run the prior predictive simulation through likelihood
# Dummy variable for land use
X_land <- c(rep(0, nag), rep(1, nurban))
# Simulate varying effects
gamma_psi <- array(NA, dim = c(nsite, nsim))
gamma_theta <- array(NA, dim = c(nsite, nsim))
gamma_p <- array(NA, dim = c(nsite, nsim))

for (i in 1:nsite){
  gamma_psi[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_psi)
  gamma_theta[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_theta)
  gamma_p[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_p)
}

for(i in 1:nsite){
    psi_i[i, ] <- plogis(a_psi + beta_psi * X_land[i] + gamma_psi[i, ]) # back-transform
    z_i[i, ] <- rbinom(nsim, 1, psi_i[i, ])
  
  for (j in 1:nsamp) {
    theta_ij[i, j, ] <- plogis( a_theta + gamma_theta[i, ]) 
    a_ij[i, j, ] <- rbinom( nsim, 1, z_i[i, ] * theta_ij[i, j, ])
   
    for (k in 1:nrep) {
      p_ijk[i, j, k, ] <- plogis( a_p + gamma_p[i, ])
      y_ijk[i, j, k, ] <- rbinom( nsim, 1, a_ij[ i, j, ] * p_ijk[ i, j, k, ])
      }
    }
  }
```

### 2.3.2 Pick draws at random to analyze

``` r
sel <- sample(1:1000, 1, replace = FALSE)
fake_data <- y_ijk[, , , sel] # draw the fake dataset from 1000 simulations
data.frame(
  a_psi = a_psi[sel],
  beta_psi = beta_psi,
  sigma_gamma_psi = sigma_gamma_psi[[sel]],
  a_theta = a_theta[sel],
  sigma_gamma_theta = sigma_gamma_theta[[sel]],
  a_p = a_p[sel],
  sigma_gamma_p = sigma_gamma_p[[sel]]
) %>% t()
```

    ##                          [,1]
    ## a_psi             0.487509316
    ## beta_psi          1.000000000
    ## sigma_gamma_psi   0.075116356
    ## a_theta           1.186577065
    ## sigma_gamma_theta 0.057954606
    ## a_p               1.589042886
    ## sigma_gamma_p     0.008515509

``` r
# potential combinations of alpha that can lead to all-zero capture history
alpha_potential <- expand.grid(rep(list(c(0, 1)), nsamp))

stan_d <- list(nsite = nsite, 
               nsamp = nsamp, 
               nrep = nrep, 
               X_Land = X_land,
               site = seq(1, nsite, 1),
               y = fake_data, 
               n_possible = 2 ^ nsamp, 
               alpha_potential = alpha_potential)
```

### 2.3.3 Fit Stan model

``` r
fit_mod3 <- sampling(
  object = mocc,
  data = stan_d,
  chains = 4,
  iter = 2000,
  cores = 4,
  thin = 1
  )
```

    ## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
    ## Running the chains for more iterations may help. See
    ## https://mc-stan.org/misc/warnings.html#bulk-ess

### 2.3.4 Print posterior summary

``` r
print(fit_mod3, 
      pars = c("a_psi", "a_theta", "a_p", "beta"), 
      digits_summary = 3)
```

    ## Inference for Stan model: anon_model.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##          mean se_mean    sd   2.5%   25%   50%   75% 97.5% n_eff  Rhat
    ## a_psi   1.421   0.040 0.769  0.422 0.908 1.255 1.750 3.482   368 1.005
    ## a_theta 0.708   0.015 0.288  0.118 0.515 0.723 0.900 1.260   393 1.012
    ## a_p     1.749   0.004 0.196  1.396 1.614 1.740 1.871 2.151  3099 1.002
    ## beta    1.283   0.030 0.942 -0.131 0.660 1.128 1.713 3.730   991 1.004
    ## 
    ## Samples were drawn using NUTS(diag_e) at Fri Mar  7 13:24:02 2025.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

![(0.55/0.5) ^ 2 = 1.21](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%280.55%2F0.5%29%20%5E%202%20%3D%201.21 "(0.55/0.5) ^ 2 = 1.21")
![200 \times 1.7 = 242](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;200%20%5Ctimes%201.7%20%3D%20242 "200 \times 1.7 = 242")

``` r
#first extract posteriors
posteriors_fit_3 <- rstan::extract(fit_mod3)

#plot
qplot(x = posteriors_fit_3$psi[ ,1 ], geom = 'histogram', binwidth = 0.01 ) +
  geom_vline(xintercept = plogis(a_psi[sel]), color = 'red', size = 2) +
  xlab(expression(psi))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_4-1.png)<!-- -->

``` r
qplot(x = posteriors_fit_3$theta[ ,1 ,1], geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = plogis(a_theta[sel]), color = 'red', size = 2) +
  xlab(expression(theta))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_4-2.png)<!-- -->

``` r
qplot(x = posteriors_fit_3$p[ , 1, 1 ,1], geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = plogis(a_p[sel]), color = 'red', size = 2) +
  xlab(expression(p))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_4-3.png)<!-- -->

``` r
qplot(x = posteriors_fit_3$beta, geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = beta_psi, color = 'red', size = 2) +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed", size = 2) +
  xlab(expression(Beta))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_4-4.png)<!-- -->

## 2.4 Same number of sites, lower prior ![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")

Same number of sites at
![100 \times 2 = 200](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;100%20%5Ctimes%202%20%3D%20200 "100 \times 2 = 200"),
but change prior for probability of capturing DNA in water sample,
conditional on occupancy.

``` r
nag <- 100 # number of ag sites in study
nurban <- 100 # number of urban sites in study
nsite <- nag + nurban # number of sites in study
nsamp <- 5 # number of replicate water samples per site
nrep <- 4 # number of replicate qPCR samples per water sample (per site and date)
nsim <- 1e3 # number of simulated draws from joint prior predictive distribution
```

``` r
nsite
```

    ## [1] 200

``` r
nsite * nsamp
```

    ## [1] 1000

``` r
nsite * nsamp * nrep
```

    ## [1] 4000

### 2.4.1 Intercept ![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta") prior

The old prior.

``` r
# prior for log-odds-scale intercept parameter for psi
loc_a_theta = 1.5
scale_a_theta = 0.25
a_theta <- rnorm(nsim, loc_a_theta, scale_a_theta)
```

<img src="multilevel_occupancy_design_scenario_files/figure-gfm/prior_theta_prob_5-1.png" style="display: block; margin: auto;" />

The new prior: reduced probability of capturing in water sample,
conditional on occupancy.

``` r
# prior for log-odds-scale intercept parameter for psi
loc_a_theta = -1
scale_a_theta = 0.5
a_theta <- rnorm(nsim, loc_a_theta, scale_a_theta)
```

<img src="multilevel_occupancy_design_scenario_files/figure-gfm/prior_theta_prob_5b-1.png" style="display: block; margin: auto;" />

## 2.5 Simulate draws

``` r
# Run the prior predictive simulation through likelihood
# Dummy variable for land use
X_land <- c(rep(0, nag), rep(1, nurban))
# Simulate varying effects
gamma_psi <- array(NA, dim = c(nsite, nsim))
gamma_theta <- array(NA, dim = c(nsite, nsim))
gamma_p <- array(NA, dim = c(nsite, nsim))

for (i in 1:nsite){
  gamma_psi[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_psi)
  gamma_theta[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_theta)
  gamma_p[i, ] <- rnorm(nsim, mean = 0, sd = sigma_gamma_p)
}

for(i in 1:nsite){
    psi_i[i, ] <- plogis(a_psi + beta_psi * X_land[i] + gamma_psi[i, ]) # back-transform
    z_i[i, ] <- rbinom(nsim, 1, psi_i[i, ])
  
  for (j in 1:nsamp) {
    theta_ij[i, j, ] <- plogis( a_theta + gamma_theta[i, ]) 
    a_ij[i, j, ] <- rbinom(nsim, 1, z_i[i, ] * theta_ij[i, j, ])
   
    for (k in 1:nrep) {
      p_ijk[i, j, k, ] <- plogis(a_p + gamma_p[i, ])
      y_ijk[i, j, k, ] <- rbinom(nsim, 1, a_ij[ i, j, ] * p_ijk[ i, j, k, ])
      }
    }
  }
```

### 2.5.1 Pick draws at random to analyze

``` r
sel <- sample(1:1000, 1, replace = FALSE)
fake_data <- y_ijk[, , , sel] # draw the fake dataset from 1000 simulations
data.frame(
  a_psi = a_psi[sel],
  beta_psi = beta_psi,
  sigma_gamma_psi = sigma_gamma_psi[[sel]],
  a_theta = a_theta[sel],
  sigma_gamma_theta = sigma_gamma_theta[[sel]],
  a_p = a_p[sel],
  sigma_gamma_p = sigma_gamma_p[[sel]]
) %>% t()
```

    ##                          [,1]
    ## a_psi              0.46534789
    ## beta_psi           1.00000000
    ## sigma_gamma_psi    0.03355722
    ## a_theta           -1.14893773
    ## sigma_gamma_theta  0.01114779
    ## a_p                1.35981291
    ## sigma_gamma_p      0.05752352

``` r
# potential combinations of alpha that can lead to all-zero capture history
alpha_potential <- expand.grid(rep(list(c(0, 1)), nsamp))

stan_d <- list(nsite = nsite, 
               nsamp = nsamp, 
               nrep = nrep, 
               X_Land = X_land,
               site = seq(1, nsite, 1),
               y = fake_data, 
               n_possible = 2 ^ nsamp, 
               alpha_potential = alpha_potential)
```

## 2.6 Fit Stan model

``` r
fit_mod4 <- sampling(
  object = mocc,
  data = stan_d,
  chains = 4,
  iter = 2000,
  cores = 4,
  thin = 1
  )
```

### 2.6.1 Print posterior summary

``` r
print(fit_mod4, 
      pars = c("a_psi", "a_theta", "a_p", "beta"), 
      digits_summary = 3)
```

    ## Inference for Stan model: anon_model.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##           mean se_mean    sd   2.5%    25%    50%    75%  97.5% n_eff  Rhat
    ## a_psi   -0.029   0.007 0.378 -0.677 -0.286 -0.048  0.197  0.782  3030 1.001
    ## a_theta -1.390   0.003 0.129 -1.651 -1.475 -1.387 -1.304 -1.138  2504 1.001
    ## a_p      1.430   0.002 0.123  1.199  1.343  1.428  1.508  1.686  4823 1.000
    ## beta     3.350   0.018 1.014  1.710  2.635  3.239  3.914  5.720  3181 0.999
    ## 
    ## Samples were drawn using NUTS(diag_e) at Fri Mar  7 13:39:06 2025.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

![(1.037/0.5)^2 = 4.3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281.037%2F0.5%29%5E2%20%3D%204.3 "(1.037/0.5)^2 = 4.3")
![200 \times 4.3 = 860](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;200%20%5Ctimes%204.3%20%3D%20860 "200 \times 4.3 = 860")

``` r
#first extract posteriors
posteriors_fit_4 <- rstan::extract(fit_mod4)

#plot
qplot(x = posteriors_fit_4$psi[ ,1], geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = plogis(a_psi[sel]), color = 'red', size = 2) +
  xlab(expression(psi))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_5-1.png)<!-- -->

``` r
qplot(x = posteriors_fit_4$theta[ ,1 ,1], geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = plogis(a_theta[sel]), color = 'red', size = 2) +
  xlab(expression(theta))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_5-2.png)<!-- -->

``` r
qplot(x = posteriors_fit_4$p[ , 1, 1 , 1], geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = plogis(a_p[sel]), color = 'red', size = 2) +
  xlab(expression(p))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_5-3.png)<!-- -->

``` r
qplot(x = posteriors_fit_4$beta, geom = 'histogram', binwidth = 0.01) +
  geom_vline(xintercept = beta_psi, color = 'red', size = 2) +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed", size = 2) +
  xlab(expression(beta))
```

![](multilevel_occupancy_design_scenario_files/figure-gfm/plotting_params_5-4.png)<!-- -->

The “winners curse” describes how, when you infer an effect by chance,
it tends to be over-estimated. The “winners curse” is a commonly cited
phenomenon with relevance to the replication crisis because these types
of noisy findings often get published uncritically, because of
over-reliance on and misunderstanding of p-values. In all of these cases
above, and also in the frequentist mode of inference, we want to look at
the standard error of the effect estimate relative to the known true
effect (or a hypothetical effect size of interest if looking at
estimates from real data and planning future studies). The effect
estimates were too noisy in many of these cases above. That is, relative
to the known true effect size (1). Standard errors decrease proportional
to the square of sample size.

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Nichols_etal_2008" class="csl-entry">

Nichols, James D., Larissa L. Bailey, Allan F. O’Connell Jr., Neil W.
Talancy, Evan H. Campbell Grant, Andrew T. Gilbert, Elizabeth M. Annand,
Thomas P. Husband, and James E. Hines. 2008. “Multi-Scale Occupancy
Estimation and Modelling Using Multiple Detection Methods.” *Journal of
Applied Ecology* 45 (5): 1321–29.
https://doi.org/<https://doi.org/10.1111/j.1365-2664.2008.01509.x>.

</div>

<div id="ref-Schmidt_etal_2013" class="csl-entry">

Schmidt, Benedikt R., Marc Kéry, Sylvain Ursenbacher, Oliver J. Hyman,
and James P. Collins. 2013. “Site Occupancy Models in the Analysis of
Environmental DNA Presence/Absence Surveys: A Case Study of an Emerging
Amphibian Pathogen.” *Methods in Ecology and Evolution* 4 (7): 646–53.
https://doi.org/<https://doi.org/10.1111/2041-210X.12052>.

</div>

</div>
