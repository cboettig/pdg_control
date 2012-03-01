function [xt_h,xt,x_ph] = ForwardSimulate(x0,pars,D,dev,n,H)
    OptTime = size(D,2); 
  %% FIXME these should be initialized as vectors!
    xt_h = x0;  
    xt = x0; 
    x_ph = x0;% Initial conditions
  for t = 1:OptTime-1
      [dummy,St] = min(abs(n-xt_h(end))); % Current state
      %xt_h[end] shouldn't really be xt_h[t]; 
      %this is correct only if we didn't initialize
      ht = H(D(St,t+1)); % The corresponding optimal harvest
      zt = lognrnd(1,dev); % Realised stochasticity this year
      x_ph(t) = xt_h(t)-ht; % The population size after harvesting takes place
      xt_h(t+1) = zt*f(x_ph(t), pars); % Population dynamics with harvest
      xt(t+1) = zt*f(xt_h(t), pars); % What would have happened without harvest
  end
end

