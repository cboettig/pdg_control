  delta = 0.1;  % The economic discount rate
  dev = 0.4;    % Standard deviation of the noise process
  OptTime = 25; % Stopping time 
  DISCRETISE = 100; % Gridsize

  pars = [2,4]; % Beverton-Holt parameters
  K = (pars(1)-1)/pars(2); % K is unharvested equilib 
 
   
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

  % Create the grid 
  n_vec = linspace(0,2*K,DISCRETISE); % Vector of state-spaces
  HVec = n_vec; 

  % Create the transition matrix
  SDP_Mat = determine_SDP_matrix(@f, pars, n_vec, HVec, dev);

  % Find the optimum solution
  [V1, D] = find_dp_optim(SDP_Mat, n_vec, HVec, OptTime, 0, @profit, delta); 

  % Simulate an example population
  [xt_h,xt,x_ph] = ForwardSimulate(K/2,pars,D,dev,n_vec,HVec);

  % plot the result
  draw_plots(HVec, n_vec, V1, D, xt_h, xt, x_ph, OptTime);




 












