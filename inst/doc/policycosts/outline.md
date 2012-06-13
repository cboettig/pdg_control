Introduction
------------

Optimal management of ecosystems and natural resources often involves
frequent adjustment of a control variable in response to the observed
state of a system. When this control is set in a policy-making process,
such as determining a fishing quota, it may be more costly to change the
policy in response to new information than to continue with the status
quo (Bohm 1974) (Xepapadeas 1992).

### Types of real world policy costs

-   Capital adjustment costs, smoothing Singh, Weninger, and Doyle
    (2006)

### Stochastic fisheries model

-   Model of Reed (1979)

Methods
-------

### Model setup

-   Fish population dynamics / state equation
-   Fishing profit function
-   Policy cost function
-   Stochastics shocks function
-   Discounting
-   Model Parameters

### Solution method: Value iteration SDP

-   Finite time horizon

### Formulating the penalty functions

-   [L1](https://github.com/cboettig/pdg_control/blob/master/inst/examples/policycosts/L1.md)
-   [L2](https://github.com/cboettig/pdg_control/blob/master/inst/examples/policycosts/L2.md)
-   assymetric ([free to increase
    quota](https://github.com/cboettig/pdg_control/blob/master/inst/examples/policycosts/free_increase.md)
    vs [free to decrease
    quota](https://github.com/cboettig/pdg_control/blob/master/inst/examples/policycosts/free_decrease.md))
-   [fixed
    costs](https://github.com/cboettig/pdg_control/blob/master/inst/examples/policycosts/fixed.md)
-   quadratic smoothing

### Apples to Apples comparisons

We compare the impact of the different functional forms for the penalty
function at a value of the penalty scaling parameter in each model that
induces an equivalent net present value to the stock when optimally
managed under this penalty. This value represents the expected net
present value before the costs of policy adjustment are paid.

[plot](http://farm8.staticflickr.com/7217/7258601130_c2fc0bcfa4_o.png)

Results
-------

### L1

[single replicate
example](http://farm8.staticflickr.com/7096/7258516896_5c89f034d5_o.png)

### L2

[example](http://farm8.staticflickr.com/7214/7258563112_2f5f9ffecd_o.png)

### Fixed fee

[example](http://farm8.staticflickr.com/7093/7258506664_d6235e5f8e_o.png)

Asymmetric
----------

[example](http://farm8.staticflickr.com/7076/7258432026_d6f8179f54_o.png)

How do we visualize these policies?
-----------------------------------

-   As single replicates
-   As distribution of replicates
-   As distribution of profits/cost
-   Visualize policy as
    [colormap1](https://a248.e.akamai.net/camo.github.com/c08160f9c375b916507740264bcd8be87259815e/687474703a2f2f6661726d382e737461746963666c69636b722e636f6d2f373233352f373235383531383739325f326330306365326165655f6f2e706e67)
    or
    [colormap2](https://a248.e.akamai.net/camo.github.com/b0754b0e121edb0a8b3a79d7e46c951b84d52f48/687474703a2f2f6661726d382e737461746963666c69636b722e636f6d2f373231312f373235383531393235365f663833623537363264635f6f2e706e67)
    or nothing?

### Summary values

-   mean/sd profit from harvest
-   mean/sd cost of solution % reduction in net present value from the
    cost (Induced cost)
-   mean/sd cost of changing policies (Direct costs)

### Other comparisons

-   Plot of variance vs magnitude of penalty, by penalty function.
    [variance
    trends](http://farm8.staticflickr.com/7224/6850042286_ef81b74acc_o.png)

-   Plot of autocorrelation vs magnitude of penalty, by penalty
    function. [autocorrelation
    trends](http://farm8.staticflickr.com/7248/6996165783_41c9894bdb_o.png)

Discussion
----------

-   Comparison of the optimal policy solution to the Reed model
-   Comparison between the different functional forms -- smoothing (L2)
    non-smoothing (L1), and destablizing (fixed).
-   Induced costs (relative to free adjustment) vs direct costs -- why
    are the induced costs much larger?
-   Impact of symmetric vs assymetric costs
-   Comparison to deterministic results
-   Contrast steady-state results to dynamic solutions under stochastic
    shocks.

References
==========

Bohm, Peter. 1974. “On the Effects of Policy Costs.” *The Swedish
Journal of Economics* 76 (mar): 104. doi:10.2307/3439361.
[http://www.jstor.org/stable/10.2307/3439361](http://www.jstor.org/stable/10.2307/3439361 "http://www.jstor.org/stable/10.2307/3439361").

Reed, William J. 1979. “Optimal escapement levels in stochastic and
deterministic harvesting models.” *Journal of Environmental Economics
and Management* 6 (dec): 350–363. doi:10.1016/0095-0696(79)90014-7.
[http://www.sciencedirect.com/science/article/pii/0095069679900147
http://linkinghub.elsevier.com/retrieve/pii/0095069679900147](http://www.sciencedirect.com/science/article/pii/0095069679900147 http://linkinghub.elsevier.com/retrieve/pii/0095069679900147 "http://www.sciencedirect.com/science/article/pii/0095069679900147 http://linkinghub.elsevier.com/retrieve/pii/0095069679900147").

Singh, R., Q. Weninger, and M. Doyle. 2006. “Fisheries management with
stock growth uncertainty and costly capital adjustment.” *Journal of
Environmental Economics and Management* 52 (sep): 582–599.
doi:10.1016/j.jeem.2006.02.006.
[http://linkinghub.elsevier.com/retrieve/pii/S009506960600043X](http://linkinghub.elsevier.com/retrieve/pii/S009506960600043X "http://linkinghub.elsevier.com/retrieve/pii/S009506960600043X").

Xepapadeas, A. P. 1992. “Environmental policy, adjustment costs, and
behavior of the firm.” *Journal of Environmental Economics and
Management*: 258–275.
[http://www.sciencedirect.com/science/article/pii/009506969290004G](http://www.sciencedirect.com/science/article/pii/009506969290004G "http://www.sciencedirect.com/science/article/pii/009506969290004G").
