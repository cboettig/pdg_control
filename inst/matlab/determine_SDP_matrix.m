function SDP_Mat = determine_SDP_matrix(fhandle, pars, n_vec, HVec, dev)
% Creates transition matrices corr to alternate actions    
% author: Carl Boettiger, modified from Michael Bode
% date: Nov 8, 2011
%
% Inputs are:
%  fhandle  - A function handle for the state equation, f(x,pars)
%  pars     - the parameters for the state equation f(x,pars)
%  n_vec    - the grid for the population dynamics
%  HVec     - the grid for harvest levels
%  dev      - sd of noise term
%
% Output
%  SDP_Mat  - the stochastic transition matrix
%
% Description: 
%   Create the stochastic transition matrices 
%   for transitions from any state (on grid, n_vec), to any 
%   other state, for each value in havest grid, HVec

  L_H = length(HVec);   % number of havest states 
  S = length(n_vec);    % number of states
  SDP_Mat = zeros(S, S, L_H); % initialize transition matrix
  % Cycle through all the harvest options -- VECTORIZE_ME
  for q = 1:L_H 
      h = HVec(q); % Harvest option being considered in this round of the loop
      for i = 1:S                  % Cycle through state-space -- VECTORIZE ME 
          x1 = n_vec(i);                % Pop is in state i, with abundance x1
          x2_exp = fhandle(x1 - h,pars);  % expected next abundance
          if x2_exp == 0; 
            SDP_Mat(i, 1, q) = 1; 
          else 
              % relative probability of a transition to that state
              PropChange = n_vec ./ x2_exp; 
              % Since the noise is lognormal. (an approximation of integral) 
              Prob = lognpdf(PropChange, log(1) - dev^2/2, dev);
              % Store the normalised transition prob
              SDP_Mat(i,:,q) = Prob ./ sum(Prob); 
            end
      end
  end
end


