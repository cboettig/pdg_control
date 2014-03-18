
optimal control management stochastic dynamic programming fisheries

Introduction
============

*Could use input from those most familiar with this literature to flush
this out more.*

-   Ecosystem management frequently set in terms of policy of quotas.

-   Policy relatively static and costly to change

-   Meanwhile, the natural world is variable

-   Optimal solutions typically track these shocks, resulting in
    impractical management recommendations in face of policy costs

(See Figure [fig:1] which motivates the importance of considering policy
costs by comparing a theoretically optimal policy over time to a typical
real world policy over time. Possible candidates - bluefin tuna vs.
halibut comparison to show variation).


[fig:1]

Types of real world policy costs
--------------------------------

-   A discussion of capital adjustment costs, smoothing, *e.g.*
    [@Singh2006].

-   Historical context, *e.g.* discussion of
    [@Bohm1974; @Reed1979; @Xepapadeas1992].

Stochastic fisheries model
--------------------------

-   Historical context of the Reed (1979) model.

-   Discuss the importance and relevance of stochastic, discrete time
    models, which frequently yield highly variable solutions since the
    policy tracks the shock.

Methods
=======

Model setup
-----------

Fish population dynamics / state equation. We will assume Beverton-Holt
dynamics,

$$X_{t+1} = Z_t \frac{A X_t}{1 + B X_t},$$

where $Z_t$ gives the stochastic shocks. $Z_t$ may be distributed
log-normally [^1]

$$\Pi_0(x,h) = p h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \label{profit}$$

where $h$ is the harvest level, $x$ the stock size, $h/x$ represents
fishing effort, $p$ the price per unit harvest. The coefficient $c_1$
introduces a quadratic cost to effort, a typical way to introduce
smoothing [@Singh2006]. For simplicity, we will
consider $c_0 = c_1 = 0$.

Policy cost function: L1, L2, fixed fee, asymmetric, are introduced as
modifications to the cost function that depend on the action taken in
the previous time-step. (This makes the previous action part of the
state space). All cost functions are characterized in terms of the
coupling coefficient $c_2$. (Because $c_2$ takes different units and
strength of interaction under the different functional forms, it is
necessary to calibrate the choice of this coefficient such that these
penalties can be compared directly, as described in the next section.)

$$\begin{aligned}
    \Pi_{L_1}(x_t,h_t, h_{t-1}) &= \Pi_0 + c_2 \operatorname{abs}\left( h_t - h_{t-1} \right) \label{L1} \\
    \Pi_{L_2}(x_t,h_t, h_{t-1}) &= \Pi_0 + c_2 \left( h_t - h_{t-1} \right)^2 \label{L2} \\
    \Pi_{\textrm{fixed}}(x_t,h_t, h_{t-1}) &= \Pi_0 + c_2 \mathbb{I}(h_t, h_{t-1})  \label{fixed_fee} \\
    \Pi_{\textrm{asym}}(x_t,h_t, h_{t-1}) &= \Pi_0 + c_2 \operatorname{max}\left( h_t - h_{t-1}, 0 \right) \label{asym}
  \end{aligned}$$

Where $\mathbb{I}(a,b) = 0$ for $a \neq b$ and $\mathbb{I}(a,b) = 1$ for
$a = b$.

Economic discounting, boundary conditions, constraints, Bellman
equation, SDP solution method on finite time horizon.

Choice of model parameters

Apples to apples comparisons
----------------------------

We compare the impact of the different functional forms for the penalty
function at a value of the penalty scaling parameter in each model that
induces an equivalent net present value to the stock when optimally
managed under this penalty. This value represents the expected net
present value before the costs of policy adjustment are paid,
Figure [fig:apples]

![image](figure/npv0.png) [fig:apples]

Results
=======

Characterize and compare functional forms of costs
--------------------------------------------------

Figure [fig:shapes] shows the impact the adjustment cost has on the
optimal policy. Each panel is generated against the same sequence of
random shocks so that they can be compared directly. In each case, the
optimal solution without any adjustment cost is shown in black and the
policy induced by optimization under the given adjustment penalty is
overlaid in blue. The $L_1$ panel of Figure [fig:shapes] shows a pattern
typical of the $L_1$-norm, which tends to avoid small adjustments to the
policy, resulting in long periods of constant policy followed by sudden
bursts of adjustment. This results in a relatively step-like policy
pattern. In contrast, the $L_2$-norm tracks all of the changes made by
the cost-free policy, but with smaller magnitude. This results in a
smoother curve that undershoots the larger oscillations seen in the
cost-free optimum, but results in a policy that changes incrementally
each interval. The fixed fee only makes large adjustments $L_1$, as
expected where the adjustment magnitude does not impact the cost. Where
the $L_1$ norm attempts to accomidate the large oscillations at the end
with a constant quota in between the extremes, the fixed cost takes the
opposite approach of one very large swing. Asymmetric ($L_1$) adjustment
costs attempt to track the optimal solution, and can end up with larger
fees than anticipated. * Not quite sure why we get this pattern in the
asymmetrics!*

The resulting stock dynamics are more similar to the optimal solution
(Figure S1).

* In terms of a policy set on escapement instead of harvest, our
solutions are actually more volatile than [@Reed1979]. Maybe this
leads to a discussion of what the control policy should be defined in
terms of (i.e. catch vs escapement)? Surely the simplicity of the
constant-escapement rule has led to it’s popularity in some part due
because that simplicity avoids adjustment costs? We should consider a
parameter range where the solution has less chatter? Parameters that
actually satisify the [@Reed1979] self-sustaining condition?*

![image](figure/2012-12-04-8f1e7b92fa-Figure3.png)

[fig:shapes]

* Figure 3 stays largely as is. But we also need to see the
corresponding x trajectories in a way that makes for easy comparisons
but not cluttered axes *

Responses to increasing costs
-----------------------------

* A version of Figure 4 to go in. But limit the range of the fraction on
the horizontal axis more. We need to see the corresponding variance and
autocorrelation for x as well as h. Also suggested we might want a panel
showing cross-correlation statistic as well x to h variation although
depending on what it shows it could perhaps go to SI. *

The comparisons shown in Figure [fig:shapes] are made at a fixed
magnitude of the policy cost equal to 25% of the net present value
derived from the cost-free optimal strategy. The different functional
forms may also respond differently to an increasing magnitude of the
adjustment cost. To quantify the extent to which the policy cost
decreases the volatility of the proposed policy, we calculate the
variance in the fishing quota over time under each policy. The variance
is averaged over 50 replicates, and recalculated at each of 20 different
intensities ranging from 0 (no cost) to 1 (NPV reduced to zero by the
adjustment costs)The variance is averaged over 50 replicates, and
recalculated at each of 20 different intensities ranging from 0 (no
cost) to 1 (NPV reduced to zero, by the adjustment costs), shown in
Figure [fig:summarystats]. * Correction: not actually zero, just to max
reduction (usually what is obtained by constant quota). Should indicate
where this asymptote has occurred)*.

The variance is always smallest in $L_2$, which creates relatively
smooth policies, and largest in the fixed cost, which creates highly
volatile adjustments (Figure [fig:var]). Stronger penalties further
shrink the variance of the $L_2$ norm and the smoothing quadratic costs
on effort, while provoking even higher variance in the fixed cost
solution. The $L_1$ solution does not become either more or less
variable as costs increase.

The autocorrelation of the policies also captures the qualitative
pattern observed in the individual replicates (Figure [fig:acor]). The
$L_2$ norm shows postively auto-correlated patterns, since a series of
small increases is often required in lieu of a single large adjustment
in the cost-free solution. The other penalties tend to have more
oscillatory solutions, leading to a negative autocorrelation.

* Add additional plot: Is value a linear function of policy cost, or is
there a sharp transition? * * Must be nonlinear somewhere when
adjustment costs exceed profits? *

[fig:summarystats]

Consequences of ignoring policy adjustment costs
------------------------------------------------

[fig:profits~c~osts]

-   Distribution of profits and costs under a policy-cost scenario vs
    profits and costs under Reed solution, by penalty type.
    Figure [fig:profits~c~osts].

-   Distribution of net present value under the policy-cost scenario vs
    the Reed solution, Figure [fig:npv~d~ist].

[fig:npv~d~ist]

Add two-by-two table of net present value accounting for: “Adjustment
cost in policy calculation” by “adjustment cost in reality”?

Discussion
==========

* Some possible things to discuss: could add others, could remove some.*

-   From Meeting 3: Interest in how the inclusion of these costs would
    operationally put a population at greater risk of extinction (e.g.
    if including the transaction costs increases the variance and
    autocorrelation of N when subject to the optimal policy). Should
    expect this, because (AH) as the transaction costs become infinite,
    you would never change your policy which puts you in the world of
    open loop control; TAC fixed for all time, which we know puts the
    stock at greater risk of extinction.

-   Alan on the Table looking at the impacts of $\sigma_{\epsilon}$
    should flag that might also be good to look at autocorrelated
    environmental noise and how that affects thing. Here we look at how
    the optimal policy induces more color in the stock size EVEN though
    the environmental variability going into it is white noise.

-   (JL) Interpreting (1), (2), (3), and (4) in terms of how much you
    would have available to pay the Japanese and Canadians to convince
    them to participate in an expedited process to change TACs.

-   Comparison of the optimal policy solution to the Reed model,
    traditional economic smoothing penalty - which now appears in SI

-   Comparison between the different functional forms – smoothing (L2)
    non-smoothing (L1), and de-stablizing (fixed).

-   Induced costs (relative to free adjustment) vs direct costs (paid
    for adjusting)

-   Contrast steady-state results to dynamic solutions under stochastic
    shocks.

**Key conclusions:**

1.  **Lesson for modelers**: Ignoring the reality of policy costs leads
    to less effective management (decreased value) (increased extinction
    risk?)

2.  **Lesson for policy makers**: Decreasing adjustment costs not only
    saves direct costs, but permits higher-value solutions.

References
==========

Bovenberg, AL, Goulder LH, and Gurney DJ (2005). “Efficiency Costs of
Meeting Industry-Distributional Constraints Under Environmental Permits
and Taxes” RAND Journal of Economics. 36(4): 951- 971.

Bovenberg, AL, Goulder LH, and Jacobsen MR (2008). “Costs of Alternative
Policy Instruments in the Presence of Industry Compensation
Requirements” Journal of Public Economics. 92: 1236-1253.

[^1]: though this violates the self-sustaining property
    of [@Reed1979], such that the optimal constant-escapement level
    $S$ in the stochastic model is greater than the optimal escapement
    in the deterministic scenario.
