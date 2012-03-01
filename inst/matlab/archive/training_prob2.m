% training_prob2.M


%% Some parameters for the Allee model
r = 2;
alpha = 20;
K = 100;
beta = 1;
pars = [r, alpha, K, beta];



function [n] = F(N, t, E_h, pars)
% Define the discrete-time dynamics with Allee effect & harvesting
% Args:
%   N - the starting population size. Can be applied over a vector.  
%   t - Scalar number of timesteps to iterate the dynamics. 
%   E_h - harvesting effort
%   pars - a vector of other parameters for the allee functional response
%          depends on choice of f, but generally gives (1) rate, (2) allee 
%          parameter, (3) carrying capacity, (4) fishing efficiency.  
% Returns:
%   n - a length(N) by t matrix of population size at time t.  


   %% Allee function, Beverton-Holt style  
%  f = @(x, E_h, pars) max(pars(1)*x.^pars(2)./(1+x.^pars(2)/pars(3)) - x*E_h*pars(4),0);

  %% Allee function, Ricker style
  f = @(x, E_h, pars) max(x.*exp(pars(1)*(1-x/pars(3)).*((x-pars(2))/pars(3))) - x*E_h*pars(4),0);

  %% Loop forward in time  
  n = zeros(length(N), t);
  n(:,1) = N;
  for i=1:t-1 
    n(:,i+1) = f(n(:,i), E_h, pars); 
  end
end


% show/check the dynamics
plot(1:T,F(25, T, 0, pars));


function [out] = Pi(E_h, E_s, t, pars)
% Cost/profit function on time t
% Args:
%   E_h - Harvesting effort
%   E_s - Sampling effort
%   t - time at which to evaluate cost 
%   pars - parameters for the Allee model
% Returns:
%   vector of profits at time 1 to t
 N_true = 50;
 P_0 = @(x) lognpdf(x, mu = log(N_true), sd=1/E_s);
 N = 0:500; % Integrate (sum) over a range of N values
 pop = F(N, t, E_h, pars);
 out = E_h * N .*  P_0(N) * pop;
 out(1) = out(1) - sample_cost*E_s;
end

% test function call
Pi(.01,.5,4,pars) 



% Integrated costs
function [out] = target(E_h, E_s, pars)
  R = .1; % economic discount rate
  T = 20; % Finite time horizon problem
  discount = exp(-R*1:T);
  out = -Pi(E_h, E_s, 1:T, pars)*discount';
end

target(.01, .01, pars)

%% This is the same target function I was passing to the optim routine (Nelder-Mead).  fminsearch function will do this in Matlab, but it might be worth exploring doing something more intelligent and enforcing the boundary condition at T.  Add and edit at will!







