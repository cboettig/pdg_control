function [D,SDP_Mat,xt_h] = Reed_SDP(DISCRETISE,SIM)
% REED_SDP implements a numerical version of the SDP described in:
%  Reed, W.J., 1979. Optimal Escapement Levels in Stochastic and
%  Deterministic Harvesting Models. Journal of Environmental Economics
%  and Management. 6: 350-363.
%
%  Updated: November 8 2011
%       
%  [D,SDP_Mat,xt_h,xt,x_ph] = REED_SDP(DISCRETISE,SIM)
%   
%  Inputs are:
%   DISCRETISE  - the resolution at which the state space (i.e., the 
%                 number of fish in the population) is discretised.
%   SIM         - a binary stating whether a stoch forward simulation
%                 of the population should be run and outputted. 
%
%  Outputs are: 
%    D          - The States x Tmax state-dependent optimal decision
%                   matrix
%    SDP_Mat     - The Markovian state transition matrix
%    xt_h        - The population level through time
% 
  if nargin < 2
    DISCRETISE = 100; 
    SIM = 1;
  end

  delta = 0.1;  % The economic discount rate
  dev = 0.4;    % Standard deviation of the noise process
  OptTime = 25; % Stopping time 

  pars = [2,4]; % Beverton-Holt parameters
  K = (pars(1)-1)/pars(2); % K is unharvested equilib 
 
  % Create the grid 
  n_vec = linspace(0,2*K,DISCRETISE); % Vector of state-spaces
  HVec = n_vec; 
 
  function out = profit(h) 
    % Define a Profit function
    p = 1; % price of fish
    c = 0.0001; % cost of fishing
    out = p * h - c ./ h;
  end

  % f is the Beverton-Holt function
  function x2 = f(x1,pars)
    x1 = max(0,x1);
    x2 = max(0,pars(1)*x1/(1+pars(2)*x1));
  end 

  profit_handle = @profit;
  fhandle = @f; 

  SDP_Mat = determine_SDP_matrix(@f, pars, n_vec, HVec, dev);
  [V1, D] = find_dp_optim(SDP_Mat, n_vec, HVec, OptTime, 0, @profit, delta); 
  [xt_h,xt,x_ph] = ForwardSimulate(K/2,pars,D,dev,n_vec,HVec);
  draw_plots(HVec, n_vec, V1, D, xt_h, xt, x_ph, OptTime);
% end




 












