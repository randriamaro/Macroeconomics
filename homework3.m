% ----- save output to file ----- %
% ----- Tiana Randriamaro Homework#3 ----- %
% ----- Question 4: Discrete State Growth Model ----- %
clear all
diary homework3.mout
diary off
delete homework3.mout
diary homework3.mout
more off

% ----- display date and time of computation ----- %
format short g
date
time0 = clock;

fprintf('\na) Policy function vector through value function iteration with continuous capital stock\n\n')
% Parameter values %
A = 1;
alpha = .36;
beta = .9;

% Steady state %
kss = (alpha*beta)^(1/(1-alpha));
k1bar = .95*kss;
k2bar = 1.05*kss;

% Initialization of endogeneous variables %
v = zeros(2,1);                     % v(i): value function at k(i)
tv = zeros(2,1);                    % tv(i): mapped value function at k(i)
kprime = zeros(2,1);                % kprime(i): next-period capital given current state k(i)

k = linspace(k1bar,k2bar,2)';       % grid for capital

crit = 1.0;                         % init criterion for stopping iterations
tol = 1e-6;                         % tolerance 
it = 0;                             % number of iterations

% consumption and current-period utility %
c = A*repmat(k.^(alpha),1,2)-repmat(k',2,1);    % cons(i,j): consumption given state k(i) and decision k(j)
isFeasible = (c>=0);                              % indicator of feasibility
c(~isFeasible) = 0.0;                           % gives -Inf utility
util = log(c);                              % util(i,j): current utility at state k(i) and decision k(j)

% iterate over the value function %
while crit > tol
    vVec = repmat(v',2,1);                      % vectorize next-period value
    [tv,ind] = max(util+beta*vVec,[],2);        % apply Bellman operator 
    kprime = k(ind);                            % update decision 
    crit = max(abs(tv-v));                      % criteria for tolerance
    v = tv;                                     % update value function
    it = it+1;                                  % update number of iterations
end

% policy function %
kprime

fprintf('\nb) Continuous state value function evaluated at k1bar, k2bar using closed form solution\n\n')
% v(k) = E + F log(k)
E = 1/(1-beta)*(log(A*(1-alpha*beta))+alpha*beta/(1-alpha*beta)*log(A*alpha*beta));
F = alpha/(1-alpha*beta);                       

disp(E+F*log(k));

fprintf('\nc) Plot of solutions from a and b\n\n')

figure1 = figure('Name','Value Function');
plot1 = plot(k,[v,E+F*log(k)]);
set(plot1(1),'Marker','o','DisplayName','Discrete');
set(plot1(2),'Marker','.','Color',[1 0 0],'DisplayName','Continuous');
xlabel('Capital');
ylabel('Value of Capital');
title('Value Function');
legend('show');
set(legend,'Location','NorthWest');

figure2 = figure('Name','Policy Function');
plot2 = plot(k,[kprime,alpha*beta*A*k.^(alpha)]);
set(plot2(1),'Marker','o','DisplayName','Discrete');
set(plot2(2),'Marker','.','Color',[1 0 0],'DisplayName','Continuous');
xlabel('Capital');
ylabel('Next-period Capital');
title('Policy Function');
legend('show');
set(legend,'Location','NorthWest');


fprintf('\nd) Plot of discrete-state and continuous state solution\n\n')
% parameter values
L = 0.6;                         % factor lower bound
U = 1.4;                         % factor upper bound
n = 100;                         % number of gridpoints

k = linspace(L*kss,U*kss,n)';

% solving for valued function and policy function through iteration
% initialization of endogeneous variables
v = zeros(n,1);                     % v(i): value function at k(i)
tv = zeros(n,1);                    % tv(i): mapped value function at k(i)
kprime = zeros(n,1);                % kprime(i): next-period capital given current state k(i)

% consumption and current-period utility
c = A*repmat(k.^(alpha),1,n)-repmat(k',n,1);    % cons(i,j): consumption given state k(i) and decision k(j)
isFeasible = (c>=0);                              % indicator of feasibility
c(~isFeasible) = 0.0;                           % gives -Inf utility
util = log(c);                                  % util(i,j): current utility at state k(i) and decision k(j)

% iterate over the value function
while crit > tol
    vVec = repmat(v',n,1);                      % vectorize next-period value
    [tv,ind] = max(util+beta*vVec,[],2);        % apply Bellman operator 
    kprime = k(ind);                            % update decision 
    crit = max(abs(tv-v));                      % criteria for tolerance
    v = tv;                                     % update value function
    it = it+1;                                  % updare number of iterations
end
kprime;
v;

% continuous state
E = 1/(1-beta)*(log(A*(1-alpha*beta))+alpha*beta/(1-alpha*beta)*log(A*alpha*beta));
F = alpha/(1-alpha*beta);  

% graph
figure3 = figure('Name','Value function n100');
plot3 = plot(k,[v,E+F*log(k)]);
set(plot3(1),'Marker','o','DisplayName','Discrete');
set(plot3(2),'Marker','.','Color',[1 0 0],'DisplayName','Continuous');
xlabel('Capital');
ylabel('Value of Capital');
title('Value function n100');
legend('show');
set(legend,'Location','NorthWest');

figure4 = figure('Name','Policy function n100');
plot4 = plot(k,[kprime,alpha*beta*A*k.^(alpha)]);
set(plot4(1),'Marker','o','DisplayName','Discrete');
set(plot4(2),'Marker','.','Color',[1 0 0],'DisplayName','Continuous');
xlabel('Capital');
ylabel('Next-period Capital');
title('Policy function n100');
legend('show');
set(legend,'Location','NorthWest');

% ----- turn off diary ----- %
comptime = etime(clock, time0)
diary off