function [V1, D] = find_dp_optim(SDP_Mat, n_vec, HVec, OptTime, xT, profit, delta) 
% Identify the dynamic optimum using backward iteration  
  L_H = length(HVec);   % number of havest states 
  S = length(n_vec);    % number of states
  V = zeros(S,1);       %  No scrap value
  V(n_vec >= xT) = 10;   % a non-zero reward for having >= xT fish at the end
  V1 = zeros(S,S);        % initial allocation
  D = zeros(S, OptTime);  % initial allocation
  for Timestep = 1:OptTime
      for i = 1:L_H % Try out all the potential harvest rates  -- VECTORIZE ME
        % Consider the ongoing reward
        sale = profit(min(HVec(i), n_vec'));
        sale = min(HVec(i), n_vec');
        V1(:,i) = SDP_Mat(:,:,i) * V + sale * exp(- delta * (OptTime - Timestep)); 
      end
      [V, D(:,OptTime-Timestep + 1)] = max(V1, [], 2); % Choose the optimal harvest
  end
end


