function [D,SDP_Mat,xt_h] = Reed_SDP_with_changecosts(DISCRETISE,SIM,ChangeCost,OptTime)
% REED_SDP_with_changecosts implements a numerical version of a variant of the SDP described in:
%       Reed, W.J., 1979. Optimal Escapement Levels in Stochastic and
%       Deterministic Harvesting Models. Journal of Environmental Economics
%       and Management. 6: 350-363.
%
%       Reed's original formulation has been changed to include a cost
%       associated with increasing or decreasing the harvest.
%
%       Updated: December 7 2011
% 
%   [D,SDP_Mat,xt_h,xt,x_ph] = REED_SDP_with_changecosts(DISCRETISE,SIM,ChangeCost,OptTime)
%   
%   Inputs are:
%           DISCRETISE  - the resolution at which the state space (i.e., the 
%                         number of fish in the population) is discretised.
%           SIM         - a binary stating whether a stochastic forward simulation
%                         of the population should be run and outputted. 
%           ChangeCost  - the cost of changing the harvest rate. Value is
%                         in the units of harvest, and is per-unit-harvest-change. 
%           OptTime     - length of the project period 
%
%   % Outputs 
%           D           - The #States x Tmax state-dependent optimal decision matrix
%           SDP_Mat     - The markovian state transition matrix
%           xt_h        - The population level through time
% 

if nargin < 3; DISCRETISE = 30; SIM = 1; ChangeCost = 0.05; OptTime = 25; disp('Running defaults'); end
A = 1.5; B = 4; % Beverton-Holt parameters
K = (A-1)/B; % K is the equilibrium population in the absence of stochasticity and harvest
delta = 0.05; % The economic discount rate
dev = 0.1; % Standard deviation of the noise process
n_vec = linspace(0,1.75*K,DISCRETISE); % Vector of population sizes

% Details of all the states. Col 1 = current abundance. Col 2 = previous% harvest rate
States = [repmat(n_vec',DISCRETISE,1) reshape(repmat(n_vec,DISCRETISE,1),DISCRETISE^2,1)]; 
S = length(States); % Number of states

% CREATE TRANSITION MATRICES CORRESPONDING TO ALTERNATE ACTIONS
SDP_Mat = zeros(S,S,DISCRETISE);
for q = 1:DISCRETISE % Cycle through all the harvest options
    h = n_vec(q); % Harvest option being considered in this round of the loop
    for i = 1:DISCRETISE % Pop is in state i
        x1 = States(i,1); % Pop has abundance x1
        x2_exp = SubBevHolt(x1,A,B,h); % What the next abundance is expected to be
        if x2_exp == 0; SDP_Mat(i+[0:DISCRETISE:DISCRETISE^2-DISCRETISE],(q-1)*DISCRETISE+1,q) = 1; else % If all transitions are to negative abundances
            PropChange = n_vec./x2_exp; % The proportional change required for a given state transition
            Prob = lognpdf(PropChange,0,dev); % The associated probability of each transition according to the Lognormal PDF
            SDP_Mat(i,(q-1)*DISCRETISE + [1:DISCRETISE],q) = Prob./sum(Prob); % Store the normalised transition probability in the transition matrix
            SDP_Mat(i+[0:DISCRETISE:DISCRETISE^2-DISCRETISE],(q-1)*DISCRETISE + [1:DISCRETISE],q) = ...
                repmat(Prob./sum(Prob),DISCRETISE,1); % Store the normalised transition probability in the transition matrix
        end
    end
end

V = zeros(S,1); % No scrap value
% IDENTIFY THE DYNAMIC OPTIMUM USING BACKWARD ITERATION 
for Timestep = 1:OptTime
    for i = 1:DISCRETISE % Try out all the potential harvest rates
        V1(:,i) = SDP_Mat(:,:,i)*V ... % The value of expected future dynamics
                  %% shouldn't this be the min of HVec and n_vec?  I guess n_vec = HVec, but still
                  + min(n_vec(i),States(:,1))*exp(-delta*(OptTime-Timestep)) ... % The value of the current harvest
                  - ChangeCost*abs(n_vec(i)-States(:,2)); % The cost of changing the harvest strategy
    end
    [V,D(:,OptTime-Timestep+1)] = max(V1,[],2); % Choose the optimal harvest rate
end

if SIM > 0
    if SIM < 5
        figure(1), clf % PLOTTING FUNCTIONS
        subplot(3,2,1), cla, hold on, xlim([1 OptTime]), set(gca,'fontsize',10), ZZ = n_vec(D); ZZ(ZZ==0)=nan;
        P=pcolor([OptTime:-1:1],n_vec,ZZ((0*DISCRETISE)+[1:DISCRETISE],end:-1:1)); set(P,'edgecolor','none'), colorbar
        title('h_{t-1} = 0','fontsize',10), axis tight
        
        subplot(3,2,2), cla, hold on, xlim([1 OptTime]), set(gca,'fontsize',10), ZZ = n_vec(D); ZZ(ZZ==0)=nan;
        P=pcolor([OptTime:-1:1],n_vec,ZZ((2*DISCRETISE)+[1:DISCRETISE],end:-1:1)); set(P,'edgecolor','none'), colorbar, axis tight
        title(['h_{t-1} = ' num2str(n_vec(2),2)],'fontsize',10)
        subplot(3,2,3), cla, hold on, xlim([1 OptTime]), set(gca,'fontsize',10), ZZ = n_vec(D); ZZ(ZZ==0)=nan;
        P=pcolor([OptTime:-1:1],n_vec,ZZ((6*DISCRETISE)+[1:DISCRETISE],end:-1:1)); set(P,'edgecolor','none'), colorbar, axis tight
        title(['h_{t-1} = ' num2str(n_vec(6),2)],'fontsize',10)
        ylabel('Population','fontsize',10), axis tight
        subplot(3,2,4), cla, hold on, xlim([1 OptTime]), set(gca,'fontsize',10), ZZ = n_vec(D); ZZ(ZZ==0)=nan;
        P=pcolor([OptTime:-1:1],n_vec,ZZ((10*DISCRETISE)+[1:DISCRETISE],end:-1:1)); set(P,'edgecolor','none'), colorbar, axis tight
        title(['h_{t-1} = ' num2str(n_vec(10),2)],'fontsize',10)
        xlabel('Management timeline (years)','fontsize',10)
        
        subplot(3,2,[5 6]), hold on, set(gca,'fontsize',14)
        [xt_h,xt,x_ph] = ForwardSimulate(K/2,A,B,D,dev,n_vec,n_vec);
    else
    	figure(2), subplot(2,2,SIM-4), cla, hold on, xlim([1 OptTime]), set(gca,'fontsize',14)
        ZZ = n_vec(D); ZZ(ZZ==0)=nan;
        P=pcolor([OptTime:-1:1],n_vec,ZZ((0*DISCRETISE)+[1:DISCRETISE],end:-1:1)); set(P,'edgecolor','none'), colorbar
        title(['Penalty = ' num2str(ChangeCost,2) 'x magnitude'],'fontsize',14)
        xlabel('Management timeline (years)','fontsize',14)
        ylabel('Population','fontsize',14), axis tight, figure(1)
    end
end
save ResultsOfChangingCost
function x2 = SubBevHolt(x1,A,B,h)
x1_minus_h = max(0,x1-h);
x2 = max(0,A*(x1_minus_h)/(1+B*(x1_minus_h)));

function [xt_h,xt,x_ph] = ForwardSimulate(x0,A,B,D,dev,n,H)
OptTime = size(D,2); xt_h = x0; xt = x0; x_ph = x0;% Initial conditions
for t = 1:OptTime-1
    [~,St] = min(abs(n-xt_h(end))); % Current state
    ht = H(D(St,t+1)); % The corresponding optimal harvest
    zt = dev*randn+1; % Realised stochasticity this year
    zt = dev*lognrnd(0,dev)+1; % Realised stochasticity this year
    x_ph(t) = xt_h(t)-ht; % The population size after harvesting takes place
    xt_h(t+1) = zt*max(0,A*x_ph(t)/(1+B*x_ph(t))); % Population dynamics with harvest
    xt(t+1) = zt*max(0,A*xt_h(t)/(1+B*xt_h(t))); % What would have happened without harvest
end
MixX = [1:length(xt); 1:length(xt)]; MixX = MixX(1:end-1);
Mix = [xt;xt_h]; Mix = Mix(1:end-1); plot(MixX,Mix,'g--','linewidth',2)
plot([1:OptTime],xt_h,'linewidth',2)
plot([1:OptTime],[x_ph xt_h(end)],'m--','linewidth',2), ylim([0 1.25*max([xt_h,xt])])
xlabel('Management timeline (years)','fontsize',16)
ylabel('Population','fontsize',16)
L=legend('N(t) without harvest','N(t)','Post-harvest population',2); set(L,'fontsize',12)
