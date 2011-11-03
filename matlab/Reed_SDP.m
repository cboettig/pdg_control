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
%           SIM         - a binary stating whether a stochastic forward simulation
%                         of the population should be run and outputted. 
%
%   % Outputs 
%           D           - The #States x Tmax state-dependent optimal decision matrix
%           SDP_Mat     - The markovian state transition matrix
%           xt_h        - The population level through time
% 


  if nargin < 2; 
    DISCRETISE = 100; 
    SIM = 1; 
    disp('Running defaults'); 
  end

  A = 2; 
  B = 4;        % Beverton-Holt parameters
  K = (A-1)/B;  % K is the equilib in the absence of stochasticity and harvest
  delta = 0.1;  % The economic ddiscount rate
  dev = 0.4;    % Standard deviation of the noise process
  n_vec = linspace(0,2*K,DISCRETISE);       % Vector of state-spaces
  S = length(n_vec);                        % number of states
  d = n_vec(2)-n_vec(1);                    % grid size 
  HVec = n_vec; 
  L_H = length(HVec); 
  SDP_Mat = zeros(S,S,L_H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE TRANSITION MATRICES CORRESPONDING TO ALTERNATE ACTIONS               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Cycle through all the harvest options -- VECTORIZE_ME
  for q = 1:L_H 
      h = HVec(q); % Harvest option being considered in this round of the loop
      for i = 1:S                  % Cycle through state-space -- VECTORIZE ME 
          x1 = n_vec(i);                % Pop is in state i, with abundance x1
          x2_exp = SubBevHolt(x1,A,B,h);             % expected next abundance
          if x2_exp == 0; 
            SDP_Mat(i,1,q) = 1; 
          else % If all transitions are to negative abundances
              PropChange = n_vec./x2_exp; 
% The associated probability of each transition according to the Lognormal PDF
              Prob = lognpdf(PropChange,0,dev);
% Store the normalised transition probability in the transition matrix
              SDP_Mat(i,:,q) = Prob./sum(Prob); 
            end
      end
  end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IDENTIFY THE DYNAMIC OPTIMUM USING BACKWARD ITERATION                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  V = 0.*n_vec'; % No scrap value
  OptTime = 25;
  V1 <- zeros(S,S)        % initial allocation
  D <- zeros(S, OptTime)  % initial allocation
  for Timestep = 1:OptTime
      for i = 1:L_H % Try out all the potential harvest rates  -- VECTORIZE ME
        % Consider the ongoing reward
        V1(:,i) = SDP_Mat(:,:,i)*V + min(HVec(i),n_vec') * exp(-delta*(OptTime-Timestep)); 
      end
      [V,D(:,OptTime-Timestep+1)] = max(V1,[],2); % Choose the optimal harvest rate
  end






  if SIM == 1
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
      [xt_h,xt,x_ph] = ForwardSimulate(K/2,A,B,D,dev,n_vec,HVec);
  end

  function x2 = SubBevHolt(x1,A,B,h)
  x1_minus_h = max(0,x1-h);
  x2 = max(0,A*(x1_minus_h)/(1+B*(x1_minus_h)));

  function [xt_h,xt,x_ph] = ForwardSimulate(x0,A,B,D,dev,n,H)
    OptTime = size(D,2); 
  %% FIXME these should be initialized as vectors!
    xt_h = x0;  
    xt = x0; 
    x_ph = x0;% Initial conditions
  for t = 1:OptTime-1
      [dummy,St] = min(abs(n-xt_h(end))); % Current state
      %xt_h[end] shouldn't really be xt_h[t] ?  doesn't matter if we didn't initialize
      ht = H(D(St,t+1)); % The corresponding optimal harvest
      zt = dev*randn+1; % Realised stochasticity this year
      x_ph(t) = xt_h(t)-ht; % The population size after harvesting takes place
      xt_h(t+1) = zt*max(0,A*x_ph(t)/(1+B*x_ph(t))); % Population dynamics with harvest
      xt(t+1) = zt*max(0,A*xt_h(t)/(1+B*xt_h(t))); % What would have happened without harvest
  end


  ReedThreshold = n(sum(D(:,1)==1));
  plot([0 OptTime],[ReedThreshold ReedThreshold],'r:')
  MixX = [1:length(xt); 1:length(xt)]; MixX = MixX(1:end-1);
  Mix = [xt;xt_h]; Mix = Mix(1:end-1); plot(MixX,Mix,'g--','linewidth',2)
  plot([1:OptTime],xt_h,'linewidth',2)
  plot([1:OptTime],[x_ph xt_h(end)],'m--','linewidth',2), ylim([0 1.1*max([xt_h,xt])])
  xlabel('Management timeline (years)','fontsize',16)
  ylabel('Population','fontsize',16), axis tight
  L=legend('Reed''s S','N(t) without harvest','N(t)','Post-harvest population',2); set(L,'fontsize',12)
