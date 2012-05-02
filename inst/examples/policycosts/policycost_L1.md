




# L1 Policy Costs 
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

## Setup the system



This example illustrates the impact of adding a cost to changing the harvest level between years 

### Define all parameters 



we'll use log normal noise functions




Chose the state equation / population dynamics function



and we use a harvest-based profit function with default parameters



Set up the discrete grids for stock size and havest levels




### Calculate the stochastic transition matrix
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations. Note that this only includes uncertainty in the growth rate (projected stock next year). 
<pre class="knitr"><div class="source">    <span class="symbol">SDP_Mat</span> <span class="assignement">&lt;-</span> <span class="functioncall">determine_SDP_matrix</span><span class="keyword">(</span><span class="symbol">f</span><span class="keyword">,</span> <span class="symbol">pars</span><span class="keyword">,</span> <span class="symbol">x_grid</span><span class="keyword">,</span> <span class="symbol">h_grid</span><span class="keyword">,</span> <span class="symbol">sigma_g</span> <span class="keyword">)</span>
    <span class="symbol">opt</span> <span class="assignement">&lt;-</span> <span class="functioncall">find_dp_optim</span><span class="keyword">(</span><span class="symbol">SDP_Mat</span><span class="keyword">,</span> <span class="symbol">x_grid</span><span class="keyword">,</span> <span class="symbol">h_grid</span><span class="keyword">,</span> <span class="symbol">OptTime</span><span class="keyword">,</span> <span class="symbol">xT</span><span class="keyword">,</span>
                     <span class="symbol">profit</span><span class="keyword">,</span> <span class="symbol">delta</span><span class="keyword">,</span> <span class="argument">reward</span><span class="argument">=</span><span class="symbol">reward</span><span class="keyword">)</span>
</div></pre>

### Find the optimum by dynamic programming 

I've updated the algorithm to allow an arbitrary penalty function. Must be a function of the harvest and previous harvest. 
<pre class="knitr"><div class="source"><span class="symbol">c2</span> <span class="assignement">&lt;-</span> <span class="number">4</span>
<span class="symbol">L1</span> <span class="assignement">&lt;-</span> <span class="keyword">function</span><span class="keyword">(</span><span class="formalargs">c2</span><span class="keyword">)</span> <span class="keyword">function</span><span class="keyword">(</span><span class="formalargs">h</span><span class="keyword">,</span> <span class="formalargs">h_prev</span><span class="keyword">)</span>  <span class="symbol">c2</span> <span class="keyword">*</span> <span class="functioncall">abs</span><span class="keyword">(</span><span class="symbol">h</span> <span class="keyword">-</span> <span class="symbol">h_prev</span><span class="keyword">)</span>
<span class="symbol">policycost</span> <span class="assignement">&lt;-</span> <span class="functioncall">optim_policy</span><span class="keyword">(</span><span class="symbol">SDP_Mat</span><span class="keyword">,</span> <span class="symbol">x_grid</span><span class="keyword">,</span> <span class="symbol">h_grid</span><span class="keyword">,</span> <span class="symbol">OptTime</span><span class="keyword">,</span> <span class="symbol">xT</span><span class="keyword">,</span>
                    <span class="symbol">profit</span><span class="keyword">,</span> <span class="symbol">delta</span><span class="keyword">,</span> <span class="symbol">reward</span><span class="keyword">,</span> <span class="argument">penalty</span> <span class="argument">=</span> <span class="functioncall">L1</span><span class="keyword">(</span><span class="symbol">c2</span><span class="keyword">)</span><span class="keyword">)</span>
</div><div class="error">Error: object 'x_grid' not found
</div></pre>



### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above.  We use a modified simulation function that can simulate an alternate policy (the Reed optimum, where policy costs are zero, `opt$D` ) and a focal policy, `policycost$D`

<pre class="knitr"><div class="source"><span class="symbol">sims</span> <span class="assignement">&lt;-</span> <span class="functioncall">lapply</span><span class="keyword">(</span><span class="number">1</span><span class="keyword">:</span><span class="number">100</span><span class="keyword">,</span> <span class="keyword">function</span><span class="keyword">(</span><span class="formalargs">i</span><span class="keyword">)</span>
  <span class="functioncall">simulate_optim</span><span class="keyword">(</span><span class="symbol">f</span><span class="keyword">,</span> <span class="symbol">pars</span><span class="keyword">,</span> <span class="symbol">x_grid</span><span class="keyword">,</span> <span class="symbol">h_grid</span><span class="keyword">,</span> <span class="symbol">x0</span><span class="keyword">,</span> <span class="symbol">policycost</span><span class="keyword">$</span><span class="symbol">D</span><span class="keyword">,</span> <span class="symbol">z_g</span><span class="keyword">,</span> <span class="symbol">z_m</span><span class="keyword">,</span> <span class="symbol">z_i</span><span class="keyword">,</span> <span class="symbol">opt</span><span class="keyword">$</span><span class="symbol">D</span><span class="keyword">,</span> <span class="argument">profit</span><span class="argument">=</span><span class="symbol">profit</span><span class="keyword">,</span> <span class="argument">penalty</span><span class="argument">=</span><span class="functioncall">L1</span><span class="keyword">(</span><span class="symbol">c2</span><span class="keyword">)</span><span class="keyword">)</span>
  <span class="keyword">)</span>
</div><div class="error">Error: object 'policycost' not found
</div></pre>



Make data tidy (melt), fast (data.tables), and nicely labeled.
<pre class="knitr"><div class="source"><span class="symbol">dat</span> <span class="assignement">&lt;-</span> <span class="functioncall">melt</span><span class="keyword">(</span><span class="symbol">sims</span><span class="keyword">,</span> <span class="argument">id</span><span class="argument">=</span><span class="functioncall">names</span><span class="keyword">(</span><span class="symbol">sims</span><span class="keyword">[[</span><span class="number">1</span><span class="keyword">]</span><span class="keyword">]</span><span class="keyword">)</span><span class="keyword">)</span>
<span class="symbol">dt</span> <span class="assignement">&lt;-</span> <span class="functioncall">data.table</span><span class="keyword">(</span><span class="symbol">dat</span><span class="keyword">)</span>
<span class="functioncall">setnames</span><span class="keyword">(</span><span class="symbol">dt</span><span class="keyword">,</span> <span class="string">"L1"</span><span class="keyword">,</span> <span class="string">"reps"</span><span class="keyword">)</span> <span class="comment"># names are nice</span>
</div></pre>


### Plots 

A single replicate, alternate dynamics should show the Reed optimum, while harvest/fishstock should show the impact of having policy costs. 
<pre class="knitr"><div class="source"><span class="functioncall">ggplot</span><span class="keyword">(</span><span class="functioncall">subset</span><span class="keyword">(</span><span class="symbol">dt</span><span class="keyword">,</span><span class="symbol">reps</span>==<span class="number">1</span><span class="keyword">)</span><span class="keyword">)</span> <span class="keyword">+</span>
  <span class="functioncall">geom_line</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="symbol">time</span><span class="keyword">,</span> <span class="symbol">alternate</span><span class="keyword">)</span><span class="keyword">)</span> <span class="keyword">+</span>
  <span class="functioncall">geom_line</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="symbol">time</span><span class="keyword">,</span> <span class="symbol">fishstock</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">col</span><span class="argument">=</span><span class="string">"darkblue"</span><span class="keyword">)</span> <span class="keyword">+</span>
  <span class="functioncall">geom_line</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="symbol">time</span><span class="keyword">,</span> <span class="symbol">harvest</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">col</span><span class="argument">=</span><span class="string">"purple"</span><span class="keyword">)</span> <span class="keyword">+</span>
  <span class="functioncall">geom_line</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="symbol">time</span><span class="keyword">,</span> <span class="symbol">harvest_alt</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">col</span><span class="argument">=</span><span class="string">"darkgreen"</span><span class="keyword">)</span>
</div><img src="http://farm8.staticflickr.com/7072/6987838346_6ef5dfd8de_o.png" class="plot" />
</pre>


A second replicate

<pre class="knitr"><div class="source"><span class="functioncall">ggplot</span><span class="keyword">(</span><span class="functioncall">subset</span><span class="keyword">(</span><span class="symbol">dt</span><span class="keyword">,</span><span class="symbol">reps</span>==<span class="number">2</span><span class="keyword">)</span><span class="keyword">)</span> <span class="keyword">+</span>
  <span class="functioncall">geom_line</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="symbol">time</span><span class="keyword">,</span> <span class="symbol">alternate</span><span class="keyword">)</span><span class="keyword">)</span> <span class="keyword">+</span>
  <span class="functioncall">geom_line</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="symbol">time</span><span class="keyword">,</span> <span class="symbol">fishstock</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">col</span><span class="argument">=</span><span class="string">"darkblue"</span><span class="keyword">)</span> <span class="keyword">+</span>
  <span class="functioncall">geom_line</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="symbol">time</span><span class="keyword">,</span> <span class="symbol">harvest</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">col</span><span class="argument">=</span><span class="string">"purple"</span><span class="keyword">)</span> <span class="keyword">+</span>
  <span class="functioncall">geom_line</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="symbol">time</span><span class="keyword">,</span> <span class="symbol">harvest_alt</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">col</span><span class="argument">=</span><span class="string">"darkgreen"</span><span class="keyword">)</span>
</div><img src="figure/rep2.png" class="plot" />
</pre>



We can visualize the equilibrium policy for each possible harvest:

<pre class="knitr"><div class="source"><span class="symbol">policy</span> <span class="assignement">&lt;-</span> <span class="functioncall">sapply</span><span class="keyword">(</span><span class="number">1</span><span class="keyword">:</span><span class="functioncall">length</span><span class="keyword">(</span><span class="symbol">h_grid</span><span class="keyword">)</span><span class="keyword">,</span> <span class="keyword">function</span><span class="keyword">(</span><span class="formalargs">i</span><span class="keyword">)</span> <span class="symbol">policycost</span><span class="keyword">$</span><span class="symbol">D</span><span class="keyword">[[</span><span class="symbol">i</span><span class="keyword">]</span><span class="keyword">]</span><span class="keyword">[</span><span class="keyword">,</span><span class="number">1</span><span class="keyword">]</span><span class="keyword">)</span>
<span class="functioncall">ggplot</span><span class="keyword">(</span><span class="functioncall">melt</span><span class="keyword">(</span><span class="symbol">policy</span><span class="keyword">)</span><span class="keyword">)</span> <span class="keyword">+</span>
  <span class="functioncall">geom_point</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="symbol">h_grid</span><span class="keyword">[</span><span class="symbol">Var2</span><span class="keyword">]</span><span class="keyword">,</span> <span class="keyword">(</span><span class="symbol">x_grid</span><span class="keyword">[</span><span class="symbol">Var1</span><span class="keyword">]</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">col</span><span class="argument">=</span><span class="symbol">h_grid</span><span class="keyword">[</span><span class="symbol">value</span><span class="keyword">]</span><span class="keyword">-</span><span class="symbol">h_grid</span><span class="keyword">[</span><span class="symbol">Var2</span><span class="keyword">]</span><span class="keyword">)</span><span class="keyword">)</span> <span class="keyword">+</span>
    <span class="functioncall">labs</span><span class="keyword">(</span><span class="argument">x</span> <span class="argument">=</span> <span class="string">"prev harvest"</span><span class="keyword">,</span> <span class="argument">y</span> <span class="argument">=</span> <span class="string">"fishstock"</span><span class="keyword">)</span> <span class="keyword">+</span>
      <span class="functioncall">scale_colour_gradientn</span><span class="keyword">(</span><span class="argument">colours</span> <span class="argument">=</span> <span class="functioncall">rainbow</span><span class="keyword">(</span><span class="number">4</span><span class="keyword">)</span><span class="keyword">)</span>
</div><img src="http://farm9.staticflickr.com/8160/6987838596_5e98df06b7_o.png" class="plot" />
</pre>


Here we plot previous harvest against the recommended harvest, coloring by stocksize.  Note this swaps the y axis from above with the color density.  Hence each x-axis value has all possible colors, but they map down onto a subset of optimal harvest values (depending on their stock). 
<pre class="knitr"><div class="source"><span class="symbol">policy</span> <span class="assignement">&lt;-</span> <span class="functioncall">sapply</span><span class="keyword">(</span><span class="number">1</span><span class="keyword">:</span><span class="functioncall">length</span><span class="keyword">(</span><span class="symbol">h_grid</span><span class="keyword">)</span><span class="keyword">,</span> <span class="keyword">function</span><span class="keyword">(</span><span class="formalargs">i</span><span class="keyword">)</span> <span class="symbol">policycost</span><span class="keyword">$</span><span class="symbol">D</span><span class="keyword">[[</span><span class="symbol">i</span><span class="keyword">]</span><span class="keyword">]</span><span class="keyword">[</span><span class="keyword">,</span><span class="number">1</span><span class="keyword">]</span><span class="keyword">)</span>
<span class="functioncall">ggplot</span><span class="keyword">(</span><span class="functioncall">melt</span><span class="keyword">(</span><span class="symbol">policy</span><span class="keyword">)</span><span class="keyword">)</span> <span class="keyword">+</span>
  <span class="functioncall">geom_point</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="symbol">h_grid</span><span class="keyword">[</span><span class="symbol">Var2</span><span class="keyword">]</span><span class="keyword">,</span> <span class="keyword">(</span><span class="symbol">h_grid</span><span class="keyword">[</span><span class="symbol">value</span><span class="keyword">]</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">col</span> <span class="argument">=</span> <span class="symbol">x_grid</span><span class="keyword">[</span><span class="symbol">Var1</span><span class="keyword">]</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">position</span><span class="argument">=</span><span class="functioncall">position_jitter</span><span class="keyword">(</span><span class="argument">w</span><span class="argument">=</span><span class="number">.005</span><span class="keyword">,</span><span class="argument">h</span><span class="argument">=</span><span class="number">.005</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">alpha</span><span class="argument">=</span><span class="number">.5</span><span class="keyword">)</span> <span class="keyword">+</span>
    <span class="functioncall">labs</span><span class="keyword">(</span><span class="argument">x</span> <span class="argument">=</span> <span class="string">"prev harvest"</span><span class="keyword">,</span> <span class="argument">y</span> <span class="argument">=</span> <span class="string">"harvest"</span><span class="keyword">)</span> <span class="keyword">+</span>
      <span class="functioncall">scale_colour_gradientn</span><span class="keyword">(</span><span class="argument">colours</span> <span class="argument">=</span> <span class="functioncall">rainbow</span><span class="keyword">(</span><span class="number">4</span><span class="keyword">)</span><span class="keyword">)</span>
</div><img src="http://farm8.staticflickr.com/7064/6987838958_199e6c3018_o.png" class="plot" />
</pre>



### Profits
<pre class="knitr"><div class="source"><span class="symbol">dt</span> <span class="assignement">&lt;-</span> <span class="functioncall">data.table</span><span class="keyword">(</span><span class="symbol">dt</span><span class="keyword">,</span> <span class="argument">id</span><span class="argument">=</span><span class="number">1</span><span class="keyword">:</span><span class="functioncall">dim</span><span class="keyword">(</span><span class="symbol">dt</span><span class="keyword">)</span><span class="keyword">[</span><span class="number">1</span><span class="keyword">]</span><span class="keyword">)</span>
<span class="symbol">profits</span> <span class="assignement">&lt;-</span> <span class="symbol">dt</span><span class="keyword">[</span><span class="keyword">,</span> <span class="functioncall">profit</span><span class="keyword">(</span><span class="symbol">fishstock</span><span class="keyword">,</span> <span class="symbol">harvest</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">by</span><span class="argument">=</span><span class="symbol">id</span><span class="keyword">]</span>
</div></pre>


Merge in profits to data.table (should be a way to avoid having to do these joins?)
<pre class="knitr"><div class="source"><span class="functioncall">setkey</span><span class="keyword">(</span><span class="symbol">dt</span><span class="keyword">,</span> <span class="symbol">id</span><span class="keyword">)</span>
<span class="functioncall">setkey</span><span class="keyword">(</span><span class="symbol">profits</span><span class="keyword">,</span> <span class="symbol">id</span><span class="keyword">)</span>
<span class="symbol">dt</span> <span class="assignement">&lt;-</span> <span class="symbol">dt</span><span class="keyword">[</span><span class="symbol">profits</span><span class="keyword">]</span>
<span class="functioncall">setnames</span><span class="keyword">(</span><span class="symbol">dt</span><span class="keyword">,</span> <span class="string">"V1"</span><span class="keyword">,</span> <span class="string">"profits"</span><span class="keyword">)</span>
</div></pre>


merge in total profits to data.table
<pre class="knitr"><div class="source"><span class="symbol">total_profit</span> <span class="assignement">&lt;-</span> <span class="symbol">dt</span><span class="keyword">[</span><span class="keyword">,</span><span class="functioncall">sum</span><span class="keyword">(</span><span class="symbol">profits</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">by</span><span class="argument">=</span><span class="symbol">reps</span><span class="keyword">]</span>
<span class="functioncall">setkey</span><span class="keyword">(</span><span class="symbol">total_profit</span><span class="keyword">,</span> <span class="symbol">reps</span><span class="keyword">)</span>
<span class="functioncall">setkey</span><span class="keyword">(</span><span class="symbol">dt</span><span class="keyword">,</span> <span class="symbol">reps</span><span class="keyword">)</span>
<span class="symbol">dt</span> <span class="assignement">&lt;-</span> <span class="symbol">dt</span><span class="keyword">[</span><span class="symbol">total_profit</span><span class="keyword">]</span>
<span class="functioncall">setnames</span><span class="keyword">(</span><span class="symbol">dt</span><span class="keyword">,</span> <span class="string">"V1"</span><span class="keyword">,</span> <span class="string">"total.profit"</span><span class="keyword">)</span>
</div></pre>


<pre class="knitr"><div class="source"><span class="functioncall">ggplot</span><span class="keyword">(</span><span class="symbol">dt</span><span class="keyword">,</span> <span class="functioncall">aes</span><span class="keyword">(</span><span class="symbol">total.profit</span><span class="keyword">)</span><span class="keyword">)</span> <span class="keyword">+</span> <span class="functioncall">geom_histogram</span><span class="keyword">(</span><span class="argument">alpha</span><span class="argument">=</span><span class="number">.8</span><span class="keyword">)</span>
</div><img src="http://farm8.staticflickr.com/7178/7133923039_0363e7c02d_o.png" class="plot" />
</pre>


<pre class="knitr"><div class="source"><span class="functioncall">save</span><span class="keyword">(</span><span class="argument">list</span><span class="argument">=</span><span class="functioncall">ls</span><span class="keyword">(</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">file</span><span class="argument">=</span><span class="string">"L1.rda"</span><span class="keyword">)</span>
</div></pre>


The mean dynamics of the state
<pre class="knitr"><div class="source"><span class="symbol">stats</span> <span class="assignement">&lt;-</span> <span class="symbol">dt</span><span class="keyword">[</span> <span class="keyword">,</span> <span class="functioncall">mean_sdl</span><span class="keyword">(</span><span class="symbol">fishstock</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">by</span> <span class="argument">=</span> <span class="symbol">time</span><span class="keyword">]</span>
<span class="functioncall">ggplot</span><span class="keyword">(</span><span class="symbol">stats</span><span class="keyword">)</span> <span class="keyword">+</span>   <span class="functioncall">geom_ribbon</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="argument">x</span> <span class="argument">=</span> <span class="symbol">time</span><span class="keyword">,</span> <span class="argument">ymin</span> <span class="argument">=</span> <span class="symbol">ymin</span><span class="keyword">,</span> <span class="argument">ymax</span> <span class="argument">=</span> <span class="symbol">ymax</span><span class="keyword">)</span><span class="keyword">,</span>
                <span class="argument">fill</span> <span class="argument">=</span> <span class="string">"darkblue"</span><span class="keyword">,</span> <span class="argument">alpha</span> <span class="argument">=</span> <span class="number">0.2</span><span class="keyword">,</span> <span class="argument">dat</span><span class="argument">=</span><span class="symbol">stats</span><span class="keyword">)</span> <span class="keyword">+</span>
                <span class="functioncall">geom_line</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="argument">x</span><span class="argument">=</span><span class="symbol">time</span><span class="keyword">,</span> <span class="argument">y</span><span class="argument">=</span><span class="symbol">y</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">lwd</span><span class="argument">=</span><span class="number">1</span><span class="keyword">)</span>
</div><img src="http://farm8.staticflickr.com/7080/7133923325_d7e4a88e85_o.png" class="plot" />
</pre>


The mean dynamics of the control
<pre class="knitr"><div class="source"><span class="symbol">stats</span> <span class="assignement">&lt;-</span> <span class="symbol">dt</span><span class="keyword">[</span> <span class="keyword">,</span> <span class="functioncall">mean_sdl</span><span class="keyword">(</span><span class="symbol">harvest</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">by</span> <span class="argument">=</span> <span class="symbol">time</span><span class="keyword">]</span>
<span class="functioncall">ggplot</span><span class="keyword">(</span><span class="symbol">stats</span><span class="keyword">)</span> <span class="keyword">+</span>  <span class="functioncall">geom_ribbon</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="argument">x</span> <span class="argument">=</span> <span class="symbol">time</span><span class="keyword">,</span> <span class="argument">ymin</span> <span class="argument">=</span> <span class="symbol">ymin</span><span class="keyword">,</span> <span class="argument">ymax</span> <span class="argument">=</span> <span class="symbol">ymax</span><span class="keyword">)</span><span class="keyword">,</span>
                <span class="argument">fill</span> <span class="argument">=</span> <span class="string">"darkblue"</span><span class="keyword">,</span> <span class="argument">alpha</span> <span class="argument">=</span> <span class="number">0.2</span><span class="keyword">)</span> <span class="keyword">+</span>
                <span class="functioncall">geom_line</span><span class="keyword">(</span><span class="functioncall">aes</span><span class="keyword">(</span><span class="argument">x</span><span class="argument">=</span><span class="symbol">time</span><span class="keyword">,</span> <span class="argument">y</span><span class="argument">=</span><span class="symbol">y</span><span class="keyword">)</span><span class="keyword">,</span> <span class="argument">lwd</span><span class="argument">=</span><span class="number">1</span><span class="keyword">)</span>
</div><img src="http://farm8.staticflickr.com/7213/7133923541_62efc05029_o.png" class="plot" />
</pre>

