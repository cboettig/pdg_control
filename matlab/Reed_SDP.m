function [D,SDP_Mat,xt_h] = Reed_SDP(DISCRETISE,SIM)
% REED_SDP implements a numerical version of the SDP described in:
%       Reed, W.J., 1979. Optimal Escapement Levels in Stochastic and
%       Deterministic Harvesting Models. Journal of Environmental Economics
%       and Management. 6: 350-363.
%
%       Updated: November 1 2011
% 
%   [D,SDP_Mat,xt_h,xt,x_ph] = REED_SDP(DISCRETISE,SIM)
%   
%   Inputs are:
%           DISCRETISE  - the resolution at which the state space (i.e., the 
%                         number of fish in the population) is discretised.
%           SIM         - a binary stating whether a stoch forward simulation
%                         of the population should be run and outputted. 
%
%   Outputs are: 
%           D           - The #States x Tmax state-dependent optimal decision
%                         matrix
%           SDP_Mat     - The markovian state transition matrix
%           xt_h        - The population level through time
% 
%  if nargin < 2
    DISCRETISE = 100; 
    SIM = 1;
%  end
  pars = [2,4]; % Beverton-Holt parameters
  K = (pars(1)-1)/pars(2); % K is the equilib w/o stochasticity & harvest
  delta = 0.1;  % The economic discount rate
  dev = 0.4;    % Standard deviation of the noise process
  OptTime = 25; % Stopping time 

  % Create the grid 
  n_vec = linspace(0,2*K,DISCRETISE); % Vector of state-spaces
  HVec = n_vec; 
 

  function out = profit(h) 
    % Define a Profit function
    p = 1; % price of fish
    c = 0.0001; % cost of fishing
    out = p*h - c/h;
  end

  % f is the Beverton-Holt function
  function x2 = f(x1,pars)
    x1 = max(0,x1);
    x2 = max(0,pars(1)*x1/(1+pars(2)*x1));
  end 

  fhandle = @f; 
  SDP_Mat = determine_SDP_matrix(fhandle, pars, n_vec, HVec, dev);
  [V, D] = find_dp_optim(SDP_Mat, n_vec, HVec, OptTime, 0, profit, delta); 
  [xt_h,xt,x_ph] = ForwardSimulate(K/2,pars,D,dev,n,H);
  draw_plots(Hvec, n_vec, V1, D, xt_h, xt, x_ph);
end



function SDP_Mat = determine_SDP_matrix(fhandle, pars, n_vec, HVec, dev)
% CREATE TRANSITION MATRICES CORRESPONDING TO ALTERNATE ACTIONS               
  L_H = length(HVec);   % number of havest states 
  S = length(n_vec);    % number of states
  SDP_Mat = zeros(S,S,L_H); % initialize transition matrix
  % Cycle through all the harvest options -- VECTORIZE_ME
  for q = 1:L_H 
      h = HVec(q); % Harvest option being considered in this round of the loop
      for i = 1:S                  % Cycle through state-space -- VECTORIZE ME 
          x1 = n_vec(i);                % Pop is in state i, with abundance x1
          x2_exp = fhandle(x1-h,pars);  % expected next abundance
          if x2_exp == 0; 
            SDP_Mat(i,1,q) = 1; 
          else 
              % relative probability of a transition to that state
              PropChange = n_vec./x2_exp; 
              % lognormal due to multiplicative Gaussian noise
              Prob = lognpdf(PropChange,0,dev);
              % Store the normalised transition prob
              SDP_Mat(i,:,q) = Prob./sum(Prob); 
            end
      end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IDENTIFY THE DYNAMIC OPTIMUM USING BACKWARD ITERATION                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [V, D] = find_dp_optim(SDP_Mat, n_vec, HVec, OptTime, xT, profit, delta) 
  L_H = length(HVec);   % number of havest states 
  S = length(n_vec);    % number of states
  V = ones(S,1) * xT;    % Scrap value 
  V1 = zeros(S,S);        % initial allocation
  D = zeros(S, OptTime);  % initial allocation
  for Timestep = 1:OptTime
      for i = 1:L_H % Try out all the potential harvest rates  -- VECTORIZE ME
        % Consider the ongoing reward
        V1(:,i) = SDP_Mat(:,:,i)*V + profit(min(HVec(i),n_vec')) * exp(-delta*(OptTime-Timestep)); 
      end
      [V,D(:,OptTime-Timestep+1)] = max(V1,[],2); % Choose the optimal harvest
  end
end

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
      zt = dev*randn+1; % Realised stochasticity this year
      x_ph(t) = xt_h(t)-ht; % The population size after harvesting takes place
      xt_h(t+1) = zt*f(x_ph(t), pars); % Population dynamics with harvest
      xt(t+1) = zt*f(xt_h(t), pars); % What would have happened without harvest
  end
end




function draw_plots(Hvec, n_vec, V1, D, xt_h, xt, x_ph)
  figure(1), clf % PLOTTING FUNCTIONS
  subplot(2,2,1), hold on, set(gca,'fontsize',14)
  P=pcolor(HVec,n_vec,V1); set(P,'edgecolor','none'), colorbar, ylim([0 75])
  title('Reward in 1st timestep','fontsize',14)
  xlabel('Alternative harvest rates','fontsize',14)
  ylabel('Population','fontsize',14), axis tight
  subplot(2,2,2), cla, hold on, xlim([1 OptTime]), set(gca,'fontsize',14)
  ZZ = HVec(D); ZZ(ZZ==0)=nan;
  P=pcolor([OptTime:-1:1],n_vec,ZZ(:,end:-1:1)); set(P,'edgecolor','none'), colorbar, ylim([0 75])

  title('Optimal decision (no scrap value)','fontsize',14)
  xlabel('Management timeline (years)','fontsize',14)
  ylabel('Population','fontsize',14), axis tight
  subplot(2,2,[3,4]), hold on, set(gca,'fontsize',14)

  ReedThreshold = n(sum(D(:,1)==1));
  plot([0 OptTime],[ReedThreshold ReedThreshold],'r:')
  MixX = [1:length(xt); 1:length(xt)]; MixX = MixX(1:end-1);
  Mix = [xt;xt_h]; Mix = Mix(1:end-1); plot(MixX,Mix,'g--','linewidth',2)
  plot([1:OptTime],xt_h,'linewidth',2)
  plot([1:OptTime],[x_ph xt_h(end)],'m--','linewidth',2), ylim([0 1.1*max([xt_h,xt])])
  xlabel('Management timeline (years)','fontsize',16)
  ylabel('Population','fontsize',16), axis tight
  L=legend('Reed''s S','N(t) without harvest','N(t)','Post-harvest population',2); set(L,'fontsize',12)
end






