% Department of Ecology and Evolutionary Biology, University of
  Tennessee, Knoxville, TN 37996, USA
% P. R Armsworth
% 

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

[ht]

[fig:1]

Problem description resolved in Meeting 3 - note the list of possible
causal mechanisms here (reaching agreement to processor plants) was
negotiated with some not being acceptable suggestions:

Empirics suggest not always as responsive as Reed solution. Explore
conditions under which that makes sense. Idea is that regulator’s
objective sometimes may not be ‘raw’ NPV, but rather NPV with some
penalties on fast changes. Possible reasons pure administrative
transaction cost of changing policies (reaching agreement), fishermens’
preferences for less variable quotas, processing plants’ preferences for
less variable quotas especially if tied into contracts for canning and
so forth. But we do not know what the functional form is and will be a
lot of work to estimate it. Here we explore several plausible candidates
that each reflect different ways these costs could work. Would be good
to know if differences between functional forms is particularly
important.

Types of real world policy costs
--------------------------------

-   A discussion of capital adjustment costs, smoothing, *e.g.*
    \citet{Singh2006}.

-   Historical context, *e.g.* discussion of
    \citep{Bohm1974, Reed1979, Xepapadeas1992}.

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

Resolved in Meeting 3: Base case for all analyses $ph-c_0E$ where
$h=qEx$ and $c_0>0$.

Resolved in Meeting 3: The economic smoothing formulation:
$ph - (c_0 E+c_4 E^2)$ will ONLY appear in SI with some figures for
comparison. We do not make a big play out of comparing one with the
other

$$\Pi_0(x,h) = p h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \label{profit}$$

where $h$ is the harvest level, $x$ the stock size, $h/x$ represents
fishing effort, $p$ the price per unit harvest. The coefficient $c_1$
introduces a quadratic cost to effort, a typical way to introduce
smoothing \citep[\emph{e.g.}][]{Singh2006}. For simplicity, we will
consider $c_0 = c_1 = 0$.

Policy cost function: L1, L2, fixed fee, asymmetric, are introduced as
modifications to the cost function that depend on the action taken in
the previous time-step. (This makes the previous action part of the
state space). All cost functions are characterized in terms of the
coupling coefficient $c_2$. (Because $c_2$ takes different units and
strength of interaction under the different functional forms, it is
necessary to calibrate the choice of this coefficient such that these
penalties can be compared directly, as described in the next section.)

Resolved in Meeting 3: We will only do comparisons described last time
as $L_1$, $L_2$ and fixed. They have new and more descriptive names now.
$c_0>0$ in all comparisons please. We do not present the asymmetric
comparison.

Resolved in Meeting 3: penalties to apply to h(t) NOT E(t)

$$\begin{aligned}
    \Pi_{L_1}(x_t,h_t, h_{t-1}) &= \Pi_0 + c_2 \operatorname{abs}\left( h_t - h_{t-1} \right)  \\
    \Pi_{L_2}(x_t,h_t, h_{t-1}) &= \Pi_0 + c_2 \left( h_t - h_{t-1} \right)^2  \\
    \Pi_{\textrm{fixed}}(x_t,h_t, h_{t-1}) &= \Pi_0 + c_2 \mathbb{I}(h_t, h_{t-1})  \\
    \Pi_{\textrm{asym}}(x_t,h_t, h_{t-1}) &= \Pi_0 + c_2 \operatorname{max}\left( h_t - h_{t-1}, 0 \right) 
  \end{aligned}$$

Where $\mathbb{I}(a,b) = 0$ for $a \neq b$ and $\mathbb{I}(a,b) = 1$ for
$a = b$.

Economic discounting, boundary conditions, constraints, Bellman
equation, SDP solution method on finite time horizon.

Choice of model parameters

Suggested from Meeting 3:

Give them more meaningful names. we suggested

$L_1$ - Variable cost with extent of policy change, case 1 / linear

$L_2$ - Variable cost with extent of policy change, case 2 / quadratic

fixed - Fixed cost associated with any policy change

My own suggestion is that rather than calling ALL of the relevant
parameters $c_2$, that you than call them $c_1$, $c_2$, $c_3$ and when
seeking to refer to them collectively you call them $c_i$

That is why the suggested parameter on $E^2$ penalties for the SI was
$c_4$ This section needs also to introduce the NPV function including
discounting. We also resolved in Meeting 3 that we should use a positive
discount rate.

Apples to apples comparisons
----------------------------

Resolved in meeting 3.

Let $h_0^*$ be the optimal control path for the base case $\Pi_0$
version of the NPV function, $NPV_0$.

Let $h_i^*$ be the optimal control path for one comparator case $\Pi_i$
version of the NPV function, $NPV_i$, e.g. $\Pi_i$ might now be what we
called the L1 penalty.

Then our independent variable in Fig. 4 for example is the fraction:
$1 - NPV_i(h_i^*) / NPV_0(h_0^*)$

By solving that and using the approach outlined in Figure 2, we obtain
the $c_i$ value that we use for a given comparison.

Given what is shown in Figure 2, we should make the focal range of
fraction be 0-30% of the $NPV_0(h_0^*)$ value to avoid a blow-up on
$c_i$ for some penalty function forms

Resolved in Meeting 3: we would do this instead of a dimensionality
argument OR trying to calibrate on costs arising only from induced
changes to controls.

The fact that there is no uniform metric for comparing different
functional forms of the penalty function complicates identifying the
implications of various penalty functions on optimal management. Put
another way, it is unclear what parameter values should be used for each
functional form in order to compare the relative effect of each penalty
function on optimal management. To address this issue, consider the
following rule for comparing functional forms inspired by Bovenberg,
Goulder and Gurnery (2005) and Bovenberg, Goulder and Jacobsen (2008).
By definition, adding a constraint, in the form of a penalty function,
to the social planner’s objective function will reduce the NPV of the
resource in addition to affecting the optimal management policy. For
example, each individual penalty function equation (3)-(6) is multiplied
by a constant C2. As the magnitude of the penalty function increases
(e.g., C2 increases), the NPV value of the resource decreases.
Importantly, though, the rate at which the NPV is affected by a given
increase in magnitude of the penalty function is different for different
penalty functions. One way to compare the relative impact of penalty
functions on optimal management is to choose parameters for the penalty
function that imply the same level of NPV conditional on optimal
management.

Figure 2 shows this comparison criterion graphically. Each curve plots
the percentage change in NPV relative to the unconstrained problem for a
given penalty function over different magnitudes of that penalty
function, C2, given optimal management. There are two important features
in Figure 2. First, the optimal management policy is not constant across
penalty function. This will be discussed in detail below. Second, the
rate of change in NPV as a function of the magnitude of the penalty
function varies both within and across policies. For example, the
existence of a relatively small L2 penalty function dramatically affects
the resource’s value initially but there is little marginal effect on
NPV for relatively large L2 penalty functions. Conversely, the rate of
change in NPV for the fixed penalty function is roughly constant over
all ranges of penalty function magnitudes.

To compare the impact of penalty functions on optimal management, we
select penalty magnitudes that make the resource worth 75% of its
unconstrained value. The dashed vertical lines in Figure 2 map the
needed penalty function magnitude, C2, to each penalty function such
that the resource is worth 75% of its unconstrained value when optimally
managed. We choose 75% somewhat arbitrarily but the method is applicable
to any NPV level and the qualitative effects on management are robust
for other levels.

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
fees than anticipated. Not quite sure why we get this pattern in the
asymmetrics!

The resulting stock dynamics are more similar to the optimal solution
(Figure S1).

In terms of a policy set on escapement instead of harvest, our solutions
are actually more volatile than \citet{Reed1979}. Maybe this leads to a
discussion of what the control policy should be defined in terms of
(i.e. catch vs escapement)? Surely the simplicity of the
constant-escapement rule has led to it’s popularity in some part due
because that simplicity avoids adjustment costs? We should consider a
parameter range where the solution has less chatter? Parameters that
actually satisify the \citet{Reed1979} self-sustaining condition?

![image](figure/2012-12-04-8f1e7b92fa-Figure3.png)

[fig:shapes]

Figure 3 stays largely as is. But we also need to see the corresponding
x trajectories in a way that makes for easy comparisons but not
cluttered axes

Responses to increasing costs
-----------------------------

A version of Figure 4 to go in. But limit the range of the fraction on
the horizontal axis more. We need to see the corresponding variance and
autocorrelation for x as well as h. Also suggested we might want a panel
showing cross-correlation statistic as well x to h variation although
depending on what it shows it could perhaps go to SI.

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
Figure [fig:summarystats]. Correction: not actually zero, just to max
reduction (usually what is obtained by constant quota). Should indicate
where this asymptote has occurred).

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

Add additional plot: Is value a linear function of policy cost, or is
there a sharp transition? Must be nonlinear somewhere when adjustment
costs exceed profits?

[fig:summarystats]

Consequences of ignoring policy adjustment costs
------------------------------------------------

**Option 1.**

We don’t include do a section 3.3. I’m not convinced the paper needs
one. Some of this I think is a tension between a focus in the paper on :

A. “What are the effects of including policy adjustment costs on optimal
policy recommendations at all?”

Vs.

B. “Given that there will be costs to adjusting policies but we don’t
know just what form these will take, what are the consequences of
different cost formulations”

If we focus the paper on question B there is less need for a section 3.3
at least as Carl has that section currently entitled (“Consequences of
ignoring policy adjustment costs” – which is not really a statement
about functional forms).

**Option 2.**

We include a version of the results Carl currently shows in his Fig 5
and 6 for the L1 penalty, but we include it also for the other penalty
function forms - e.g. using Carl’s suggestion of a box-whisker plot.

Reminder - The two figures compare the performance of the optimal policy
when accounting for the policy adjustment costs (maximizing $\Pi_{L_1}$)
with what the optimal policy recommendation one would arrive if ignoring
these costs (maximizing $\Pi_0$) WHEN scoring the performance of the two
control recommendations against $\Pi_{L_1}$.

Fig. 5 shows the breakdown into how the gross yield (ph) and adjustment
cost ($c_2 abs(h_t-h_{t-1})$) separately. The policy ignoring costs of
policy adjustment gets higher gross yields by tracking the variations
more closely but also incurs (even) higher policy adjustment costs. Fig
6 combines the two bits to show the effect on net present value.
Distribution of profits and costs under a policy-cost scenario vs
profits and costs under Reed solution, by penalty type.
Figure [fig:profits~c~osts].

[fig:profits~c~osts]

Keep Figure 5 or its equivalent for one penalty function formulation and
talk through what it means and how to interpret it.

Underneath lay out

$$NPV_i(h_{it}^*) = \sum \frac{p h^*_{it}-c_0E^*_{it}-c(h^*_{it})}{(1+\delta)^t}$$

ALL MINUS

$$NPV_i(h_{0t}^*) = \sum \frac{p h^*_{0t}-c_0E^*_{0t}-c(h^*_{0t})}{(1+\delta)^t}$$

Use labelling on the equations to label:

$p h^*_{it}-c_0E^*_{it} $ as contribution (1 - red in Figure 5) and
$c(h^*_{it})$ as contribution (3 - green in Figure 5) in the first
expression

and then

$p h^*_{0t}-c_0E^*_{0t}$ as contribution (2 - blue ) $c(h^*_{0t}$ as
contribution (4 - mauve) in the second expression

Then cost of assuming penalties when they are not actually there, i.e.
induced cost of controls is (2) - (1), blue - red cost of assuming no
penalties when they are present is [1-3] - [2-4]. Finally tabulate
separately for each of $\Pi_1$, $\Pi_2$, $\Pi_3$

rows: $\sigma_{\epsilon}$ low, medium, high

columns: $c_i$ values corresponding to fraction NPV or 90%, 80%, 70%

Cells: (2)-(1) AND [1-3]-[2-4] all divided by $NPV_0(h^*_{0t})$ to make
them dimensionless

(alternatively just draw three versions of Fig. 5 if not much is going
on as we scope parameter space in this way.)

This option would let one answer:

​i) What is the consequence of ignoring policy adjustment costs when
they are present?

​ii) It also answers Jim’s suggested reversal of that question, I think,
namely: what is the consequence of assuming adjustment costs are present
when in fact they are absent? (I think this is just the difference
between the gross yields / the turquoise and red shaded areas in Fig 5).

​iii) And then how do these things vary for each of the different
functional forms?

**Option 3**

We populate a table that organizes columns by the “Actual functional
form of adjustment costs” and the rows into “Assumed functional form of
adjustment costs.” The types of functional form in both columns and rows
would then run “None, $L_1$, $L_2$ etc.”

IF we denote the optimal control for objective function $NPV_i$ as being
$h^*_i$, then I would populate the Table with the proportions
$NPV_i(h^*_j) / NPV_i(h^*_i)$ (i.e. what proportion of the available NPV
you still get if assuming the wrong functional form for the adjustment
costs).

The inclusion of the comparison against None means that all of what was
being talked about for Option 2 would be captured. In addition it would
provide a robustness check on a policy recommendation if uncertain of
either the presence of policy adjustment costs or functional form that
they should adopt. Specifically, by comparing across rows one could
ascertain whether assuming one choice of functional form was likely to
lead to large efficiency costs should other assumptions about adjustment
costs prove correct. Similarly one could compare across columns to see
if one form were particularly problematic in reality.

I don’t know if those are clear – you need to get mired into the
document once again a bit to parse it out. I’ll happily discuss with
anyone in a bid to shed more light.

Discussion of these options and alternative suggestions invited.

-   Distribution of profits and costs under a policy-cost scenario vs
    profits and costs under Reed solution, by penalty type.
    Figure [fig:profits~c~osts].

-   Distribution of net present value under the policy-cost scenario vs
    the Reed solution, Figure [fig:npv~d~ist].

[fig:npv~d~ist]

Add two-by-two table of net present value accounting for: “Adjustment
cost in policy calculation” by “adjustment cost in reality”?

**Further things to look at?**

-   How are costs partitioned between policy adjustments and deviation
    from optimality?

-   Risk sensitivity: does policy inertia increase or decrease
    extinction risk? Does accounting for the policy cost increase or
    decrease extinction risk? (only in cases where extinction is
    possible: *e.g.* in large variance, outside of Reed’s
    self-sustaining range. Or Allee model)

Discussion
==========

Some possible things to discuss: could add others, could remove some.

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

[^1]: though this violates the self-sustaining property
    of \citet{Reed1979}, such that the optimal constant-escapement level
    $S$ in the stochastic model is greater than the optimal escapement
    in the deterministic scenario.
