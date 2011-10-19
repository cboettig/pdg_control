%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   stochastic_growth.m:  A Matlab program to solve a simple stochastic 
%   growth model via either value function or policy function iteration,
%   which is also know as Howard's Policy Improvement Algorithm).
%
%   Characteristics of the problem are:
%           (1) Infinite horizon
%           (2) Continuous state space and discrete time
%           (3) Income shock takes on two discrete states
%
%   Assumptions are:
%        (1) U(c)=c^(1-g)/(1-g); (2) discount rate =10%;  (3)F(K)=(1/a)*K^a
% 
% ARE 254 Dynamic Optimization Fall 2009
%
% Parts of this code are based on code writen by George Hall, July 2001 
% Yale University. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stochastic_growth

clear all
close all
format short g
disp('STOCHASTIC GROWTH MODEL');
disp('');

% Specifying which approach you want to use
valueiteration=0;    
if valueiteration==1
    disp('Value Function Iteration')
else
    disp('Policy Function Iteration')
end
 
%  Economic parameter values
a  = 0.20;               % production parameter f(K)=A*K^a          
b  = 0.10;               % discount rate
b  =1/(1+b);             % discount factor
g  = 1.50;               % g>=1, U(c) is CRRA Utility function

% Transition probability matrix
P=[.8 .2; .2 .8];    % Pij=Prob(A(t+1)=j|A(t)=i) where j is high, i is low               
A_high = 1.250;      % high value for technology
A_low  = 0.50;       % low value for technology

% Discretizing the continuous state variable, Capital.
N=1501; % N is the grid size, greater N the longer the solution will take

%%%%%%%%%
% Typically, you would solve for the steady-state level and use that to 
% determine the range. I know for this problem that the steady-state is 
% around K=1. 
%%%%%%%%%
MaxK=1.5;  
MinK=0.1;
K=linspace(MinK,MaxK,N);        % Uniform grid, which is typical default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tabulating conumption and utility for all combinations of K, K(t+1) for
% the grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For all values of K(t) and K(t+1), calculate what consumption would be
[Kt,Kt1]=meshgrid(K,K);
consumeh=A_high*(1/a)*Kt.^a-Kt1;
consumel=A_low*(1/a)*Kt.^a-Kt1;

% Replace all negative consumption with NaN.
% You will not necessarily need this depending on your problem. E.g, is
% u(c)=ln(c)
consumeh(consumeh<=0)=NaN;
consumel(consumel<=0)=NaN;

% Calculate the utility at all levels of consumption
utilityh=(consumeh.^(1-g))/(1-g);
utilityl=(consumel.^(1-g))/(1-g);

% Replace NaN in utility to -Inf to ensure that they are not chosen
% You will not necessarily need this depending on your problem. E.g, is
% u(c)=ln(c)
utilityh(isnan(utilityh))=-inf;  
utilityl(isnan(utilityl))=-inf;

% Initializing iteration variables
iter=0;
v=repmat(0,N,2);  % Takes a vector of 0's of size N and duplicates it twice
% wtf? why not zeros(N,2)
decis=repmat(0,N,2);
VError=10; 
Tol=1e-7;
    
if valueiteration==1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Iteration on value function%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while VError > Tol
        % Calculate at the value at high
        [tvh,tdecish]=max(utilityh+b*repmat(v*P(1,:)',1,N));
        % Calculate at the value at low
        [tvl,tdecisl]=max(utilityl+b*repmat(v*P(2,:)',1,N));
        % combine in one vector
        tv=[tvh' tvl'];
        tdecis=[tdecish' tdecisl'];
        % Check the value function
        VError=max(abs((tv-v)./v));
        v=tv; % Update the value function
        decis=tdecis;% Keeps track of the index of decision
        iter=iter+1;
        if iter>500, break; end
    end
    
else
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Solve for fixed point using policy iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tv=zeros(N,2); 
    
    while VError > Tol;
        
       % Calculate at the value at high
        [tvh,tdecish]=max(utilityh+b*repmat(v*P(1,:)',1,N));
        % Calculate at the value at low
        [tvl,tdecisl]=max(utilityl+b*repmat(v*P(2,:)',1,N));
        
        tdecis=[tdecish' tdecisl'];
        
        % Solving out for the consumption that corresponds to the value
        % function using the index from the above calculation of the value
        % function. x<-argmax(U(c)+bPV)
        r1 = zeros(N,1); r2=r1;
        for i=1:N
            r1(i) = utilityh(tdecish(i),i);
            r2(i) = utilityl(tdecisl(i),i);
        end
        
        % Finding the sequence of the transition probabilities
        g1=sparse(N,N); g2=g1;
        for i=1:N
            g1(i,tdecish(i))=1;
            g2(i,tdecisl(i))=1;
        end
        
        Q=[P(1,1)*g1 P(1,2)*g1; P(2,1)*g2 P(2,2)*g2];
        
        % Newton step v<-inv(I-bQ)*Utility(policy)
        % Records for both high and low in tv
        tv(:) = (speye(2*N) - b*Q)\[ r1; r2 ];
        
        % Check the value function
        VError=max(max(abs((tv-v)./tv)));
        
        % update the value function v
        v=tv;
        % Updates the index for K', which is consumption
        decis=tdecis;
        iter = iter+1;
        if iter>500, break; end
    end
end

disp(['Fixed point solved in ' num2str(iter) ' iterations']);

% Decision rule
DEC=K(decis); % pulls out the K's that correspond to the optimal K from the max

%
%   form transition matrix
%   trans is the transition matrix from state at t (row)
%   to the state at t+1 (column) 
%   The eigenvector associated with the unit eigenvalue
%   of trans' is  the stationary distribution. 
% 
 g1=sparse(N,N); g2=g1;
 for i=1:N
  g1(i,tdecish(i))=1;
  g2(i,tdecisl(i))=1;
end
Q=[P(1,1)*g1 P(1,2)*g1; P(2,1)*g2 P(2,2)*g2];
Q= Q';
probst = (1/(2*N))*ones(2*N,1); % initial guess  
test = 1;
while test > 10^(-8);
   probst1 = Q*probst; % Finding the stationary distribution using fixed
                       % point iteration (pie=Q*pie)
   test=max(abs(probst1-probst));
   probst = probst1;
end;

%   vectorize the decision rule to be conformable with probst
%   to calculate mean level of capital
meanK=probst'*DEC(:);

%
%   calculate stationary distribution of capital 
%
lambda=zeros(N,2);
lambda(:)=probst;
probk=sum(lambda');     
probk=probk';

%
%   print out results
%
disp('PARAMETER VALUES');
disp([' a=' num2str(a) ',  b=' num2str(b)]); 
disp('RESULTS ');
disp('');
disp(['Mean of K =' num2str(meanK)]);

%%%%%%%%%%Plotting the results%%%%%%%%%%%%%%%%%%%%%
figure
plot(K,v(:,1),'-',K,v(:,2),':');
%title({'Stochastic growth model','Value Function'});
title('Stochastic growth model Value Function');
ylabel('Value Function')
xlabel('Capital ')
legend('High','Low','Location','SouthEast')

figure
plot(K,DEC(:,1),'.-',K,DEC(:,2),':',K,K,'-');
%title({'Stochastic growth model','Policy Function'});
title('Stochastic growth model, Policy Function');
ylabel('K in next period')
xlabel('K in this period')
legend('High','Low','45^o line','Location','SouthEast')
axis([ MinK MaxK MinK MaxK ]);

figure
plot(K,probk);
title('DISTRIBUTION OF CAPITAL');
xlabel('Capital');
ylabel('Probability');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    simulate a number of consumption paths over n periods 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kgrid = linspace(MinK,MaxK,N)';  % capital grid  
kmark = 10;
k = kgrid(kmark,1);        % initial level of assets
n = 100;                   % number of periods to simulate
s0 = 1;                    % initial state, high=1, low=0
states   = zeros(n-1,2);
controls = zeros(n-1,2);
runs=100;

for J=1:runs
    [chain, dummy] = markov(P,n,s0);  % Code returns the chain and states
    for i = 1:n-1;
        if chain(i) == 1;           % If we realize a high income
            kprime = DEC(kmark,1);   % given the level of K calculate optimal K'
            cons   = A_high*(1/a)*k^a - kprime;  % calculate consumption
            kmark = tdecis(kmark,1); % at the level of consumption calculate new
            % level of K in next period
        elseif chain(i) == 2;       % If we realize a low income
            kprime = DEC(kmark,2);
            cons   = A_high*(1/a)*k^a - kprime;
            kmark = tdecis(kmark,2);
        else
            disp('something is wrong with chain');
        end;
        states(i,:) =   [ k chain(i) ];
        controls(i,:) = [ cons kprime ];
        k = kprime;               % update the level of K for next run
    end;
    Control{J}=controls;
end

figure
subplot(211)
plot((1:n-1)',Control{1}(:,1),'r',(1:n-1)',Control{2}(:,1),'b',...
     (1:n-1)',Control{3}(:,1),'k',(1:n-1)',Control{4}(:,1),'g')
legend('Run 1', 'Run 2', 'Run 3', 'Run 4','Orientation','Horizontal',...
       'Location','SouthEast')
title('Stochastic growth model: simulated consumption path');
ylabel('Consumption')
xlabel('Stage')

subplot(212)
Con=0;
% Calculating the average consumption in each stage across the differnt
% runs
for t=1:n-1
    for J=1:runs
        Con=Con+Control{J}(t,1);
    end
    avg(t)=(1/runs)*Con;
    Con=0;
end
plot(1:n-1,avg)
title(['Stochastic growth model: Avg.is over ' num2str(runs) ' runs']);
ylabel('Consumption')
xlabel('Stage')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chain,state]=markov(T,n,s0,V)
%function [chain,state]=markov(T,n,s0,V);
%  chain generates a simulation from a Markov chain of dimension
%  the size of T
%
%  T is transition matrix
%  n is number of periods to simulate
%  s0 is initial state
%  V is the quantity corresponding to each state
%  state is a matrix recording the number of the realized state at time t
%
%
[r c]=size(T);
if nargin == 1, V=[1:r]; s0=1; n=100; end;
if nargin == 2, V=[1:r]; s0=1; end;
if nargin == 3, V=[1:r]; end;
%
if r ~= c;
    disp('error using markov function');
    disp('transition matrix must be square');
    return;
end;
%
for k=1:r;
    if sum(T(k,:)) ~= 1;
        disp('error using markov function')
        disp(['row ',num2str(k),' does not sum to one']);
        disp(' it sums to :');
        disp([ sum(T(k,:)) ]);
        disp(['normalizing row ',num2str(k),'']);
        T(k,:)=T(k,:)/sum(T(k,:));
    end;
end;
[v1 v2]=size(V);
if v1 ~= 1 || v2 ~=r
    disp('error using markov function');
    disp(['state value vector V must be 1 x ',num2str(r),''])
    if v2 == 1 && v2 == r;
        disp('transposing state valuation vector');
        V=V';
    else
        return;
    end;
end
if s0 < 1 || s0 > r;
    disp(['initial state ',num2str(s0),' is out of range']);
    disp(['initial state defaulting to 1']);
    s0=1;
end;
%
%T
%rand('uniform');
X=rand(n-1,1);
s=zeros(r,1);
s(s0)=1;
cum=T*triu(ones(size(T)));
%
for k=1:length(X);
    state(:,k)=s;
    ppi=[0 s'*cum];
    s=((X(k)<=ppi(2:r+1)).*(X(k)>ppi(1:r)))';
end;
chain=V*state;

end
