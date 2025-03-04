NARMS qPCR design: multilevel occupancy
================
Roy Martin
March 04, 2025

# 1 Background and model

For this scenario, imagine that we are looking to estimate the
prevalence of a specific gene variant for 30 sites, visited once every
three weeks for 17 weeks (roughly 1 year). So we have
![t \in 1,...,T=17](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t%20%5Cin%201%2C...%2CT%3D17 "t \in 1,...,T=17")
roughly equally-spaced points in time for which we intend to estimate
the state of occupancy for each of the
![i \in 1,..., S=30](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%20%5Cin%201%2C...%2C%20S%3D30 "i \in 1,..., S=30")
sites. For our water sampling method, we’re also going to utilize
![j \in 1,...,J=5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j%20%5Cin%201%2C...%2CJ%3D5 "j \in 1,...,J=5")
replicate grab samples from a site on each visit; so we’d have to use
![T \times J = 17 \* 5 = 85](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T%20%5Ctimes%20J%20%3D%2017%20%2A%205%20%3D%2085 "T \times J = 17 * 5 = 85")
water bottles over the course of the study for a single site in this
scenario. Then, we’d conduct
![k \in 1,...,K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%20%5Cin%201%2C...%2CK "k \in 1,...,K")
replicate qPCR runs per each of the water samples grabbed. Lets say
we’re going to do
![K=4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%3D4 "K=4")
replicate qPCR analyses for each water sample
![j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j "j"),
where maybe we divide each, say, 500 mL sample equally among
![K=4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%3D4 "K=4")
125 mL samples for filtering and extraction. To summarize, we’d have 4
qPCR replicates, nested in each of 5 replicate water samples, which are
nested in each of 17 days for each of 30 sampling sites. So, our
![N](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;N "N")
as summed across the observation level (qPCR reps) is
![N = K \times J \times T \times S = 4 \times 5 \times 17 \times 30 = 10200](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;N%20%3D%20K%20%5Ctimes%20J%20%5Ctimes%20T%20%5Ctimes%20S%20%3D%204%20%5Ctimes%205%20%5Ctimes%2017%20%5Ctimes%2030%20%3D%2010200 "N = K \times J \times T \times S = 4 \times 5 \times 17 \times 30 = 10200").
This is a huge
![N](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;N "N")
and maybe unfeasible, but it is a place to start for illustrating the
model and getting a general feel for some of the inherent sources of
uncertainty in this type of study.

Next, we’ll use these hypothetical data dimensions to simulate a fake
dataset using a multilevel occupancy model as described in
\[@Schmidt_etal_2013\] and earlier by \[@Nichols_etal_2008\]. At first,
we will ignore any hypothetical covariates (e.g., stream size, land use,
measured water chemistry), but they can easily be included in the model
later. The basic model below consists of a sequence of three coupled
Bernoulli trials for describing the array of nested data,
![y\_{i,t,j,k}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%2Cj%2Ck%7D "y_{i,t,j,k}"),
where we have an observed detection
(![y\_{itjk} = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bitjk%7D%20%3D%201 "y_{itjk} = 1"))
or not
(![y\_{itjk} = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bitjk%7D%20%3D%200 "y_{itjk} = 0"))
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

The model can be described in mathematical notation:

![z\_{i,t} \sim Bernoulli(\psi\_{i,t})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bi%2Ct%7D%20%5Csim%20Bernoulli%28%5Cpsi_%7Bi%2Ct%7D%29 "z_{i,t} \sim Bernoulli(\psi_{i,t})")

![a\_{ij, t} \| z\_{i,t} \sim Bernoulli( z\_{i,t}\theta\_{ij,t} )](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_%7Bij%2C%20t%7D%20%7C%20z_%7Bi%2Ct%7D%20%5Csim%20Bernoulli%28%20z_%7Bi%2Ct%7D%5Ctheta_%7Bij%2Ct%7D%20%29 "a_{ij, t} | z_{i,t} \sim Bernoulli( z_{i,t}\theta_{ij,t} )")

![y\_{ijk,t} \| a\_{ij,t} \sim Bernoulli( a\_{ij,t} p\_{ijk,t} )](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bijk%2Ct%7D%20%7C%20a_%7Bij%2Ct%7D%20%5Csim%20Bernoulli%28%20a_%7Bij%2Ct%7D%20p_%7Bijk%2Ct%7D%20%29 "y_{ijk,t} | a_{ij,t} \sim Bernoulli( a_{ij,t} p_{ijk,t} )")

Here,
![z\_{i,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bi%2Ct%7D "z_{i,t}")
represents a latent random variable defining the
![\textbf{true state of occupancy}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Btrue%20state%20of%20occupancy%7D "\textbf{true state of occupancy}"),
meaning the presence
(![z\_{i,t} = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bi%2Ct%7D%20%3D%201 "z_{i,t} = 1"))
or absence
(![z\_{i,t} = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bi%2Ct%7D%20%3D%200 "z_{i,t} = 0"))
of target DNA in the sampling area of a site
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
at time
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t").
In our hypothetical, we’re assuming this snapshot (within a day) visit
is capable of inferring the
![z](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z "z")
state for that site and day. The latent state
![z](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z "z")
is assumed Bernoulli distributed and is a function of the
![\textbf{occupancy probability}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Boccupancy%20probability%7D "\textbf{occupancy probability}")
parameter
![\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi "\psi").

The next random variable down the hierarchy,
![a\_{ij,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_%7Bij%2Ct%7D "a_{ij,t}"),
is also dependent on
![z\_{i,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bi%2Ct%7D "z_{i,t}"),
and described as the
![\textbf{latent state of availability}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Blatent%20state%20of%20availability%7D "\textbf{latent state of availability}").
It defines whether or not there is DNA in a water sample
![j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j "j")
from site
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")
and time
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")
that would available for detection via qPCR. If the site is occupied on
that sampling visit (i.e,
![z\_{i,t} = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bi%2Ct%7D%20%3D%201 "z_{i,t} = 1"))
with DNA at some concentration, then
![a\_{ij,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_%7Bij%2Ct%7D "a_{ij,t}")
is determined by a Bernoulli trial and according to the parameter
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta"),
or the
![\textbf{availability probability}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bavailability%20probability%7D "\textbf{availability probability}").
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")
could also be described as the probability of capturing observable DNA
in a water sample, conditional on the site being occupied. If
![z\_{i,t}=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bi%2Ct%7D%3D0 "z_{i,t}=0")
then
![a\_{ij,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_%7Bij%2Ct%7D "a_{ij,t}")
would necessarily be
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0"),
as it would be sampled from a Bernoulli distribution with
![\theta = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D%200 "\theta = 0").

The third random variable
![y\_{ijk,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bijk%2Ct%7D "y_{ijk,t}")
is the four-dimensional array containing the observed detects
![y\_{ijk,t}  \in \\0,1\\](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bijk%2Ct%7D%20%20%5Cin%20%5C%7B0%2C1%5C%7D "y_{ijk,t}  \in \{0,1\}")
recorded at the level of the qPCR replicate. These observed states are
also dependent on the latent availability state,
![a\_{ij,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_%7Bij%2Ct%7D "a_{ij,t}"),
in the way that
![a\_{ij,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_%7Bij%2Ct%7D "a_{ij,t}")
is dependent on
![z\_{i,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bi%2Ct%7D "z_{i,t}").
If
![a\_{ij,t} = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_%7Bij%2Ct%7D%20%3D%200 "a_{ij,t} = 0"),
then
![y\_{ijk,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bijk%2Ct%7D "y_{ijk,t}")
is necessarily zero. Likewise, if
![z\_{i,t}=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_%7Bi%2Ct%7D%3D0 "z_{i,t}=0"),
then
![y\_{ijk,t}=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bijk%2Ct%7D%3D0 "y_{ijk,t}=0").
If
![a\_{ij,t} = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_%7Bij%2Ct%7D%20%3D%201 "a_{ij,t} = 1"),
the qPCR observations (binary according to cycle threshold) are assumed
Bernoulli distributed according the the parameter
![p\_{ijk,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_%7Bijk%2Ct%7D "p_{ijk,t}"),
which is the
![\textbf{probability of detection}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bprobability%20of%20detection%7D "\textbf{probability of detection}")
for qPCR rep
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k "k"),
from water sample
![j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j "j"),
in site
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i"),
on day
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t").
The
![\textbf{probability of detection}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7Bprobability%20of%20detection%7D "\textbf{probability of detection}")
could also be described as the probability of detecting DNA via qPCR
given that it is available in the water sample.

Ultimately, we have three types of parameters,
![\psi\_{i,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi_%7Bi%2Ct%7D "\psi_{i,t}"),
![\theta\_{ij,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_%7Bij%2Ct%7D "\theta_{ij,t}"),
and
![p\_{ijk,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_%7Bijk%2Ct%7D "p_{ijk,t}"),
which we would like to estimate by conditioning on the observations,
![y\_{ijk,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bijk%2Ct%7D "y_{ijk,t}"),
using fully Bayesian methods. Again, these are the probability of
occupancy, the probability of availability (given occupied), and the
probability of detection via PCR (given availability). The glaring
assumption of this model is that there are no false positives. This is
maybe a reasonable assumption for qPCR in terms of the measurement
instrument itself working faithfully. It is maybe a reasonable
assumption for the laboratory practice, where we would need to assume no
or very little contamination or cross-reactivity. Outside the
laboratory, false positives can arise in the field due to contamination
during sampling, travel, or transfer. So, aside from building a
different model that relaxes this important assumption, one would need
to closely monitor for contamination in the lab, field, and transport
(e.g., blanks, etc) in order to quantify and minimize the chance for
false positives.

``` r
knitr::opts_chunk$set(include = FALSE)
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

options(mc.cores = 16)
rstan_options(auto_write = TRUE)
options(loo.cores = 16)

options(max.print = 9999)
```

# 2 Prior predictive simulation

Next, we’ll run some simulations from this model using our beliefs about
the values of the parameters
![\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi "\psi"),
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta"),
and
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p").
We’re taking a Bayesian perspective, so we are going to assign each
parameter a probability distribution that reflects these prior beliefs.
After we settle on these priors, we can then do a “push forward” of our
model, which will effectively simulate draws from the *prior predictive
distribution*.

## 2.1 Set up conditions

Lets first set up our hypothetical data dimensions in terms of the
number of dates, water samples, pcr replicates, and simulations.

So, we’re going to sample 30 sites on 17 consecutive days. We’re going
to take 5 water samples on each of those days for each site; and we’re
going to split each of those water samples into 4 PCR replicates.
Finally, we’re going to simulate 1000 draws from the prior predictive
distribution.

Now we need to set up ‘containers’ for the simulation. That is, make
some objects in R for storing for the simulated values of
![z_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_i "z_i"),
![a\_{ij,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_%7Bij%2Ct%7D "a_{ij,t}"),
![y\_{ijk,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bijk%2Ct%7D "y_{ijk,t}"),
etc. Note that these containers all have one more dimension than
mentioned in the model description above. The added dimension is for
holding the simulated draws from the prior predictive distribution. In
this case, that last dimension will always be of size = 1000, which is
the number of draws we decided (above) to use for our simulation. This
number could just as well be
![1e4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1e4 "1e4")
or something otherwise huge, but anything over a few hundred should be
fine for the purposes of getting an understanding the statistical
properties of the posterior.

## 2.2 Models for ![\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi "\psi"), ![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta"), and ![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")

In the following, we’ll choose priors for
![\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi "\psi"),
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta"),
and
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p").

### 2.2.1 ![\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi "\psi")

Our indexing in the introduction indicates that we’ve left room for
extending the models for the parameters with additional linear
structure. In this case, we’re going to employ a linear predictor for
each parameter with a logit link function. Our model for
![\psi\_{i,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi_%7Bi%2Ct%7D "\psi_{i,t}"),
will be:

![logit( \psi\_{i,t} ) = \alpha\_\psi + \gamma\_{\psi_i} + \delta\_{\psi_t} + \nu\_{\psi\_{i,t}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;logit%28%20%5Cpsi_%7Bi%2Ct%7D%20%29%20%3D%20%5Calpha_%5Cpsi%20%2B%20%5Cgamma_%7B%5Cpsi_i%7D%20%2B%20%5Cdelta_%7B%5Cpsi_t%7D%20%2B%20%5Cnu_%7B%5Cpsi_%7Bi%2Ct%7D%7D "logit( \psi_{i,t} ) = \alpha_\psi + \gamma_{\psi_i} + \delta_{\psi_t} + \nu_{\psi_{i,t}}")

where
![\alpha\_\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_%5Cpsi "\alpha_\psi")
is just the intercept for the linear predictor of
![\psi\_{i,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi_%7Bi%2Ct%7D "\psi_{i,t}")
on the log-odds scale. The intercept parameter
![\alpha\_\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_%5Cpsi "\alpha_\psi"),
is back-transformed to the probability scale via the inverse-logit
function. So, equivalently:

![\psi_i = \frac{ 1 }{ 1 + exp( - \alpha\_\psi + \gamma\_{\psi_i} + \delta\_{\psi_t} + \nu\_{\psi\_{i,t}} ) }](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi_i%20%3D%20%5Cfrac%7B%201%20%7D%7B%201%20%2B%20exp%28%20-%20%5Calpha_%5Cpsi%20%2B%20%5Cgamma_%7B%5Cpsi_i%7D%20%2B%20%5Cdelta_%7B%5Cpsi_t%7D%20%2B%20%5Cnu_%7B%5Cpsi_%7Bi%2Ct%7D%7D%20%29%20%7D "\psi_i = \frac{ 1 }{ 1 + exp( - \alpha_\psi + \gamma_{\psi_i} + \delta_{\psi_t} + \nu_{\psi_{i,t}} ) }")

In this parameterization,
![\alpha\_\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_%5Cpsi "\alpha_\psi")
determines the average probability of occupancy across sites and time;
and
![\gamma\_{\psi_i}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cgamma_%7B%5Cpsi_i%7D "\gamma_{\psi_i}"),
![\delta\_{\psi_t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdelta_%7B%5Cpsi_t%7D "\delta_{\psi_t}"),
and
![\nu\_{\psi\_{i,t}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cnu_%7B%5Cpsi_%7Bi%2Ct%7D%7D "\nu_{\psi_{i,t}}")
are varying effects, centered on zero, that determine deviations around
the intercept parameter according to site- and date-specific effects,
and their interaction, respectively. These varying effects further
paremeterized where:

![\gamma\_{\psi_i} \sim N(0, \sigma\_{\gamma\_\psi})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cgamma_%7B%5Cpsi_i%7D%20%5Csim%20N%280%2C%20%5Csigma_%7B%5Cgamma_%5Cpsi%7D%29 "\gamma_{\psi_i} \sim N(0, \sigma_{\gamma_\psi})")

![\delta\_{\psi_t} \sim N(0, \sigma\_{\delta\_\psi})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdelta_%7B%5Cpsi_t%7D%20%5Csim%20N%280%2C%20%5Csigma_%7B%5Cdelta_%5Cpsi%7D%29 "\delta_{\psi_t} \sim N(0, \sigma_{\delta_\psi})")

![\nu\_{\psi\_{i,t}} \sim N(0, \sigma\_{\nu\_\psi})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cnu_%7B%5Cpsi_%7Bi%2Ct%7D%7D%20%5Csim%20N%280%2C%20%5Csigma_%7B%5Cnu_%5Cpsi%7D%29 "\nu_{\psi_{i,t}} \sim N(0, \sigma_{\nu_\psi})")

The standard deviation terms
![\sigma\_{\gamma\_\psi}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cgamma_%5Cpsi%7D "\sigma_{\gamma_\psi}"),
![\sigma\_{\delta\_\psi})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cdelta_%5Cpsi%7D%29 "\sigma_{\delta_\psi})"),
and
![\sigma\_{\nu\_\psi})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cnu_%5Cpsi%7D%29 "\sigma_{\nu_\psi})")
determine the extent of variation around zero for each of the varying
effects; and are conditioned on the observed data.

We must choose priors for all of the parameters to be conditioned on the
observed data, including
![\alpha\_{\psi}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_%7B%5Cpsi%7D "\alpha_{\psi}"),
![\sigma\_{\gamma}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cgamma%7D "\sigma_{\gamma}"),
![\sigma\_{\delta})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cdelta%7D%29 "\sigma_{\delta})"),
and
![\sigma\_{\nu})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cnu%7D%29 "\sigma_{\nu})").
For the intercept parameter, we choose a normal prior, where:

![\alpha\_\psi \sim N( \mu\_{\alpha\_{\psi}}, \sigma\_{\alpha\_{\psi}} )](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_%5Cpsi%20%5Csim%20N%28%20%5Cmu_%7B%5Calpha_%7B%5Cpsi%7D%7D%2C%20%5Csigma_%7B%5Calpha_%7B%5Cpsi%7D%7D%20%29 "\alpha_\psi \sim N( \mu_{\alpha_{\psi}}, \sigma_{\alpha_{\psi}} )")

This is to say that our prior for this intercept parameter is “normally
distributed with location =
![\mu\_{\alpha\_{\psi}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_%7B%5Calpha_%7B%5Cpsi%7D%7D "\mu_{\alpha_{\psi}}"),
and scale =
![\sigma\_{\alpha\_{\psi}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Calpha_%7B%5Cpsi%7D%7D "\sigma_{\alpha_{\psi}}")”.
With regard to location, recall that
![\alpha\_\psi = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_%5Cpsi%20%3D%200 "\alpha_\psi = 0")
would correspond to
![\psi = logit^{-1}( \alpha\_\psi ) =0.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi%20%3D%20logit%5E%7B-1%7D%28%20%5Calpha_%5Cpsi%20%29%20%3D0.5 "\psi = logit^{-1}( \alpha_\psi ) =0.5")
on the probability scale. So if our normal prior is centered such that
![\mu\_{\alpha\_{\psi}}=0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_%7B%5Calpha_%7B%5Cpsi%7D%7D%3D0 "\mu_{\alpha_{\psi}}=0"),
we’re thinking that the most likely true value for
![\psi_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi_i "\psi_i")
is 0.5. Practically, this would suggest that we think that our sampling
sites are occupied with the target DNA for 50% or half of all the
site-day combinations. The intercept is a bit hard to conceptualize in
this context because it requires placing a prior over all the sites and
times, when there is likely a good bit of heterogeneity across those
dimensions. To further flesh out our prior, nevertheless, we’d want to
think about how certain we are about that location and adjust
![\sigma\_{\alpha\_{\psi}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Calpha_%7B%5Cpsi%7D%7D "\sigma_{\alpha_{\psi}}")
accordingly. In practice, this process of choosing useful priors should
be done carefully and critically, employing as much domain expertise as
possible, in order to learn as much as possible from the data/study.

Lets now choose the priors for the intercept parameter, which requires
choices for
![\mu\_{\alpha\_{\psi}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_%7B%5Calpha_%7B%5Cpsi%7D%7D "\mu_{\alpha_{\psi}}")
and
![\sigma\_{\alpha\_{\psi}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Calpha_%7B%5Cpsi%7D%7D "\sigma_{\alpha_{\psi}}").
Then we can simulate the 1000 draws from this prior.

So, we’ve chosen
![\mu\_{\alpha\_{\psi}} = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_%7B%5Calpha_%7B%5Cpsi%7D%7D%20%3D%200 "\mu_{\alpha_{\psi}} = 0")
and
![\sigma\_{\alpha\_{\psi}} = 1.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Calpha_%7B%5Cpsi%7D%7D%20%3D%201.5 "\sigma_{\alpha_{\psi}} = 1.5").
Therefore, our prior is centered over
![\alpha\_\psi = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_%5Cpsi%20%3D%200 "\alpha_\psi = 0")
on the log-odds scale, or
![logit^{-1}( \alpha\_\psi ) = 0.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;logit%5E%7B-1%7D%28%20%5Calpha_%5Cpsi%20%29%20%3D%200.5 "logit^{-1}( \alpha_\psi ) = 0.5")
on the probability scale. Along with the scale we’ve chosen, this prior
places 90% of the density of
![\alpha\_\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_%5Cpsi "\alpha_\psi")
between about -2.467 and -2.467 on the log-odds scale, or 0.078 and
0.922 on the probability scale. Here’s a vizualization of the prior on
the log-odds scale.

And on the probability scale.

And, finally, here are our simulated draws for this prior on the
probability scale.

Next, we’ll choose the priors standard deviation parameters for the
varying effects. For these standard deviation parameters, the values
must be non-zero. We use half-normal priors, which again require a
choice of location and scale.

First, the prior for
![\sigma\_{\gamma\_{\psi}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cgamma_%7B%5Cpsi%7D%7D "\sigma_{\gamma_{\psi}}"),
which is the parameter determining the amount of site to site variation
in
![\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi "\psi").
Again, this variation is implemented in the model structure above as
variation around the mean probability of occupancy,
![\alpha\_{\psi}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_%7B%5Cpsi%7D "\alpha_{\psi}").

We’ve chosen the location parameter for the prior on
![\sigma\_{\gamma\_{\psi}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cgamma_%7B%5Cpsi%7D%7D "\sigma_{\gamma_{\psi}}")
to be
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0")
and the scale parameter to be
![0.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.5 "0.5").
Therefore, our prior is centered over
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0").
Along with the scale we’ve chosen, this prior places 90% of the density
of this standard deviation parameter between about 0.031 and 0.98. Below
is a visualization of the prior.

Here are our simulated draws for this prior on the probability scale.

Next, the prior for
![\sigma\_{\delta\_{\psi}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cdelta_%7B%5Cpsi%7D%7D "\sigma_{\delta_{\psi}}"),
the parameter controlling the amount of day to day variation.

We’ve chosen the location parameter for the prior on
![\sigma\_{\delta\_{\psi}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cdelta_%7B%5Cpsi%7D%7D "\sigma_{\delta_{\psi}}")
to be
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0")
and the scale parameter to be
![0.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.5 "0.5").
Our prior places most of the probability over
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0").
Along with the scale we’ve chosen, this prior places 90% of the density
of this standard deviation parameter between about 0.031 and 0.98. Below
is a visualization of the prior.

Here are our simulated draws for this prior on the probability scale.

The last parameter that needs a prior for the
![\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi "\psi")
portion of the model is
![\sigma\_{\nu\_{\psi}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cnu_%7B%5Cpsi%7D%7D "\sigma_{\nu_{\psi}}"),
which is the parameter controlling the amount of variation in the
interaction between the site and day effects. This effectively
determines how much between-site variation there is in day to day
variability. The all-encompassing term might be “spatio-temporal”
variation.

We’ve chosen the location parameter for the prior on
![\sigma\_{\nu\_{\psi}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cnu_%7B%5Cpsi%7D%7D "\sigma_{\nu_{\psi}}")
to be
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0")
and the scale parameter to be
![0.25](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.25 "0.25").
Our prior places most of the probability over
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0").
Along with the scale we’ve chosen, this prior places 90% of the density
of this standard deviation parameter between about 0.016 and 0.49. Below
is a visualization of the prior.

Here are our simulated draws for this prior on the probability scale.

### 2.2.2 ![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")

Now lets run through the same process for
![\theta\_{ij,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_%7Bij%2Ct%7D "\theta_{ij,t}"),
the availability probability. Again, we’ll use the logit-link linear
predictor with the intercept parameter, where:

![logit( \theta\_{ij,t} ) = \alpha\_\theta + \gamma\_{\theta_i} + \delta\_{\theta_t} + \nu\_{\theta\_{i,t}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;logit%28%20%5Ctheta_%7Bij%2Ct%7D%20%29%20%3D%20%5Calpha_%5Ctheta%20%2B%20%5Cgamma_%7B%5Ctheta_i%7D%20%2B%20%5Cdelta_%7B%5Ctheta_t%7D%20%2B%20%5Cnu_%7B%5Ctheta_%7Bi%2Ct%7D%7D "logit( \theta_{ij,t} ) = \alpha_\theta + \gamma_{\theta_i} + \delta_{\theta_t} + \nu_{\theta_{i,t}}")

We, again, further describe the varying effects as:

![\gamma\_{\theta\_{i}} \sim N(0, \sigma\_{\gamma\_\theta})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cgamma_%7B%5Ctheta_%7Bi%7D%7D%20%5Csim%20N%280%2C%20%5Csigma_%7B%5Cgamma_%5Ctheta%7D%29 "\gamma_{\theta_{i}} \sim N(0, \sigma_{\gamma_\theta})")

![\delta\_{\theta\_{t}} \sim N(0, \sigma\_{\delta\_\theta})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdelta_%7B%5Ctheta_%7Bt%7D%7D%20%5Csim%20N%280%2C%20%5Csigma_%7B%5Cdelta_%5Ctheta%7D%29 "\delta_{\theta_{t}} \sim N(0, \sigma_{\delta_\theta})")

![\nu\_{\theta\_{j,t}} \sim N(0, \sigma\_{\nu\_\theta})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cnu_%7B%5Ctheta_%7Bj%2Ct%7D%7D%20%5Csim%20N%280%2C%20%5Csigma_%7B%5Cnu_%5Ctheta%7D%29 "\nu_{\theta_{j,t}} \sim N(0, \sigma_{\nu_\theta})")

So, well choose prior values for
![\mu\_{\alpha\_\theta}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_%7B%5Calpha_%5Ctheta%7D "\mu_{\alpha_\theta}")
and
![\sigma\_{\alpha\_\theta}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Calpha_%5Ctheta%7D "\sigma_{\alpha_\theta}").

We’ve chosen
![\mu\_{\alpha\_\theta} = -0.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_%7B%5Calpha_%5Ctheta%7D%20%3D%20-0.5 "\mu_{\alpha_\theta} = -0.5")
and
![\sigma\_{\alpha\_\theta} = 0.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Calpha_%5Ctheta%7D%20%3D%200.5 "\sigma_{\alpha_\theta} = 0.5").
Therefore, our prior is centered over
![\alpha\_\theta = -0.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_%5Ctheta%20%3D%20-0.5 "\alpha_\theta = -0.5")
on the log-odds scale, or
![logit^{-1}( \alpha\_\theta ) = 0.378](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;logit%5E%7B-1%7D%28%20%5Calpha_%5Ctheta%20%29%20%3D%200.378 "logit^{-1}( \alpha_\theta ) = 0.378")
on the probability scale. This prior places 90% of the density of
![\alpha\_\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_%5Ctheta "\alpha_\theta")
between about -1.322 and -1.322 on the log-odds scale, or 0.21 and 0.58
on the probability scale. Here’s a visualization of the prior on the
log-odds scale.

On the probability scale.

And, finally, our simulated draws from this prior on the probability
scale.

The prior for
![\sigma\_{\gamma\_{\theta}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cgamma_%7B%5Ctheta%7D%7D "\sigma_{\gamma_{\theta}}"),
which is the parameter determining the amount of site to site variation
in
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta").

We’ve chosen the location parameter for the prior on
![\sigma\_{\gamma\_{\theta}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cgamma_%7B%5Ctheta%7D%7D "\sigma_{\gamma_{\theta}}")
to be
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0")
and the scale parameter to be
![0.25](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.25 "0.25").
Therefore, our prior is centered over
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0").
Along with the scale we’ve chosen, this prior places 90% of the density
of this standard deviation parameter between about 0.016 and 0.49. Below
is a visualization of the prior.

Here are our simulated draws for this prior on the probability scale.

Next, the prior for
![\sigma\_{\delta\_{\theta}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cdelta_%7B%5Ctheta%7D%7D "\sigma_{\delta_{\theta}}"),
the parameter controlling the amount of day to day variation.

We’ve chosen the location parameter for the prior on
![\sigma\_{\delta\_{\theta}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cdelta_%7B%5Ctheta%7D%7D "\sigma_{\delta_{\theta}}")
to be
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0")
and the scale parameter to be
![0.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.5 "0.5").
Our prior places most of the probability over
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0").
Along with the scale we’ve chosen, this prior places 90% of the density
of this standard deviation parameter between about 0.016 and 0.49. Below
is a visualization of the prior.

Here are our simulated draws for this prior on the probability scale.

The last parameter that needs a prior for the
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")
portion of the model is
![\sigma\_{\nu\_{\theta}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cnu_%7B%5Ctheta%7D%7D "\sigma_{\nu_{\theta}}"),
which is the parameter controlling the amount of “spatio-temporal”
variation.

We’ve chosen the location parameter for the prior on
![\sigma\_{\nu\_{\theta}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cnu_%7B%5Ctheta%7D%7D "\sigma_{\nu_{\theta}}")
to be
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0")
and the scale parameter to be
![0.1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.1 "0.1").
Our prior places most of the probability over
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0").
Along with the scale we’ve chosen, this prior places 90% of the density
of this standard deviation parameter between about 0.006 and 0.196.
Below is a visualization of the prior.

Here are our simulated draws for this prior on the probability scale.

### 2.2.3 ![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")

Finally, lets look at priors for
![p\_{ijk,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_%7Bijk%2Ct%7D "p_{ijk,t}").
Again, we’ll use the logit-link linear predictor and an intercept-only
model with the intercept parameter
![\alpha_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_p "\alpha_p"),
where:

![logit( p\_{ijk,t} ) = \alpha_p + \gamma\_{p_i} + \delta\_{p_t} + \nu\_{p\_{i,t}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;logit%28%20p_%7Bijk%2Ct%7D%20%29%20%3D%20%5Calpha_p%20%2B%20%5Cgamma_%7Bp_i%7D%20%2B%20%5Cdelta_%7Bp_t%7D%20%2B%20%5Cnu_%7Bp_%7Bi%2Ct%7D%7D "logit( p_{ijk,t} ) = \alpha_p + \gamma_{p_i} + \delta_{p_t} + \nu_{p_{i,t}}")

We, again, further describe the varying effects as:

![\gamma\_{p_i} \sim N(0, \sigma\_{\gamma_p})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cgamma_%7Bp_i%7D%20%5Csim%20N%280%2C%20%5Csigma_%7B%5Cgamma_p%7D%29 "\gamma_{p_i} \sim N(0, \sigma_{\gamma_p})")

![\delta\_{p_t} \sim N(0, \sigma\_{\delta_p})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdelta_%7Bp_t%7D%20%5Csim%20N%280%2C%20%5Csigma_%7B%5Cdelta_p%7D%29 "\delta_{p_t} \sim N(0, \sigma_{\delta_p})")

![\nu\_{p\_{i,t}} \sim N(0, \sigma\_{\nu_p})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cnu_%7Bp_%7Bi%2Ct%7D%7D%20%5Csim%20N%280%2C%20%5Csigma_%7B%5Cnu_p%7D%29 "\nu_{p_{i,t}} \sim N(0, \sigma_{\nu_p})")

Next, choose prior values for
![\mu\_{\alpha_p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_%7B%5Calpha_p%7D "\mu_{\alpha_p}")
and
![\sigma\_{\alpha_p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Calpha_p%7D "\sigma_{\alpha_p}").

We’ve chosen
![\mu\_{\alpha_p} = 2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_%7B%5Calpha_p%7D%20%3D%202 "\mu_{\alpha_p} = 2")
and
![\sigma\_{\alpha_p} = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Calpha_p%7D%20%3D%201 "\sigma_{\alpha_p} = 1").
Therefore, our prior is centered over
![\alpha_p = 2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_p%20%3D%202 "\alpha_p = 2")
on the log-odds scale, or
![logit^{-1}( \alpha_p ) = 0.881](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;logit%5E%7B-1%7D%28%20%5Calpha_p%20%29%20%3D%200.881 "logit^{-1}( \alpha_p ) = 0.881")
on the probability scale. This prior places 90% of the density of
![\alpha_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha_p "\alpha_p")
between about 0.355 and 0.355 on the log-odds scale, or 0.588 and 0.975
on the probability scale. Here’s a vizualization of the prior on the
log-odds scale.

On the probability scale.

And, finally, our simulated draws from this prior on the probability
scale.

The prior for
![\sigma\_{\gamma_p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cgamma_p%7D "\sigma_{\gamma_p}"),
which is the parameter determining the amount of site to site variation
in
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p").

We’ve chosen the location parameter for the prior on
![\sigma\_{\gamma_p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cgamma_p%7D "\sigma_{\gamma_p}")
to be
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0")
and the scale parameter to be
![0.1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.1 "0.1").
Therefore, our prior is centered over
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0").
Along with the scale we’ve chosen, this prior places 90% of the density
of this standard deviation parameter between about 0.006 and 0.196.
Below is a visualization of the prior.

Here are our simulated draws for this prior on the probability scale.

Next, the prior for
![\sigma\_{\delta_p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cdelta_p%7D "\sigma_{\delta_p}"),
the parameter controlling the amount of day to day variation in the
probability of detection via qPCR conditional on it being available in
the water sample.

We’ve chosen the location parameter for the prior on
![\sigma\_{\delta_p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cdelta_p%7D "\sigma_{\delta_p}")
to be
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0")
and the scale parameter to be
![0.1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.1 "0.1").
Our prior places most of the probability over
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0").
Along with the scale we’ve chosen, this prior places 90% of the density
of this standard deviation parameter between about 0.006 and 0.196.
Below is a visualization of the prior.

Here are our simulated draws for this prior on the probability scale.

The last parameter that needs a prior for the
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
portion of the model is
![\sigma\_{\nu_p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cnu_p%7D "\sigma_{\nu_p}"),
which is the parameter controlling the amount of “spatio-temporal”
variation.

We’ve chosen the location parameter for the prior on
![\sigma\_{\nu_p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma_%7B%5Cnu_p%7D "\sigma_{\nu_p}")
to be
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0")
and the scale parameter to be
![0.05](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0.05 "0.05").
Our prior places most of the probability over
![0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;0 "0").
Along with the scale we’ve chosen, this prior places 90% of the density
of this standard deviation parameter between about 0.003 and 0.098.
Below is a visualization of the prior.

Here are our simulated draws for this prior on the probability scale.

## 2.3 Simulate draws

Now we can actually simulate draws from the joint prior predictive
distribution. We’ll use our draws from the priors above above and then
push forward or loop through the likelihood.

### 2.3.1 Freqency of occupancy

We can then use these draws from the prior predictive distribution to
examine some of the practical implications of our priors and the model.
For example, we can plot the prior predictive distribution of the number
of days the sampling site is expected to be occupied.

Or the number of sites occupied on any particular day.

### 2.3.2 Frequency of detections

Or we can plot the prior predictive distribution of the total number of
qPCR runs returning a positive result (out of 1.02^{4} possible for each
simulation).

### 2.3.3 Replicate datasets

We could look at it another way by, say, visualizing replicated
hypothetical datasets from the prior predictive distribution. Lets draw
several replicate ‘datasets’ at random from the prior predictive
distribution for the first date and plot them.

Specically, we’ll plot 20 of them, along with the values of the
parameters they were generated from. Note that the parameters
(![\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi "\psi"),
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta"),
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p"))
are fixed for each of the ‘datasets’ in the plot, but each of those
panels is only one probabilistic realization from the potential set
generated from that particular configuration of the parameters.

It may be informative to generate multiple realizations from a single
configuration of the paramters. We’ll now try that by taking one of the
draws from the prior for the three parameters and plotting some
hypothetical datasets based on that configuration. Another way to think
of it might be to think of the prior for each parameter as just being a
point value (or very strong).

So, lets just take the last configuration of parameters from our
previous plot, where we have
![\psi =](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi%20%3D "\psi =")
0.67,
![\theta =](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta%20%3D "\theta =")
0.4, and
![p =](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%20%3D "p =")
0.95.

# 3 A model in Stan

Now let use this model simulation to check our computational algorithm
that we would use for Bayesian inference given real data (i.e., Stan).
We’ll code the Stan model, then we’ll select a “dataset” from the
simulation above and see if our Stan program recovers the known
parameters.

We’ll first code up our Stan model.

## 3.1 Stan code

Explanation of this code may come later.

## 3.2 Pick fake dataset for Stan

Lets just pick the last dataset from the simulation with known
parameters above. We have to creat a variable called ‘alpha_potential’
to account for all the potential combinations of
![a\_{ij}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_%7Bij%7D "a_{ij}")
that could lead to an all-zero observation history. This is in order to
integrate out the latent discrete parameters (i.e.,
![z_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;z_i "z_i"),
and
![a\_{ij}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a_%7Bij%7D "a_{ij}"))
for our Stan program because Stan doesn’t allow discrete paramters. This
marginalization is actually optimal even in BUGS, but we’ll save the
explanation for later.

In any case, we also create a data list for our Stan program.

## 3.3 Fit Stan model

Now we can fit the model to the data with Stan in R.

## 3.4 Print posterior summary

And we’ll print the posterior estimates. Note the indexing on “psi”,
“theta” and “p”, We chose just the first index for each of these because
all the others are the same due to this being an “intercept only” model.
With a more complex linear predictor, each of these parameters could
vary among sites, samples and/or reps; and in that case we might choose
to print more dimensions. More likely, though, we’d be more interested
in the logit scale parameters in that case. Nevertheless, for this
simpler model we’ll print both the logit scale intercept parameters and
each of
![\psi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpsi "\psi"),
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta"),
and
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p").

We can also plot the posteriors and overlay the known parameters to
assess recovery.
