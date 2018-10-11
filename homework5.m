
% ----- save output to file ----- %
clear all
diary hw5.mout
diary off
delete hw5.mout
diary hw5.mout
more off

% ----- display date and time of computation ----- %
format short g
date
time0 = clock;

% ---- Question #3 ----- %
% parameters
beta = .9;
alpha = .36;

% shock
z1 = .9;
z2 = 1.1;
z = [z1;z2];

% transition matrix
p11 = .95;
p22 = .9;
P = [p11 .05 
    .1 p22];
pz1 = [p11 .05];
pz2 = [.1 p22];

% steady state
kss = (alpha*beta)^(1/(1-alpha));

% capital stock
k1bar = .95*kss;
k2bar = 1.05*kss;
%k = [k1bar k2bar]';  % grid for capital
kD = linspace(k1bar,k2bar,2)';
[r c]= size(kD);

% Initialization of endogeneous variables %
v = zeros(r,2);                     % v(i): value function at k(i),z(i)
tv = zeros(r,2);                    % tv(i): mapped value function at k(i),z(i)
%kprime = zeros(r,2);                % kprime(i): next-period capital given current state k(i)

crit = 1.0;                         % init criterion for stopping iterations
tol = 1e-6;                         % tolerance 
it = 0;                             % number of iterations

while crit > tol
    for i=1:r
        k = kD(i);
        c = [z(1)*k^alpha-kD z(2)*k^alpha-kD]; 
        [tv(i,:) I] = max(log(c)+ beta*v*P');
        g(i,:) = [kD(I(1)) kD(I(2))];
    end
    crit = max(max(abs(tv-v)));
    v = tv;

end

fprintf('\nb) Decision rule for policy function \n\n')
disp([g(:,1);g(:,2)])

fprintf('\nc) Value function and policy function with n=1000 \n\n')
kD_1000 = linspace(k1bar,k2bar,1000)';
[r c]= size(kD_1000);

% Initialization of endogeneous variables %
v = zeros(r,2);                     % v(i): value function at k(i),z(i)
tv = zeros(r,2);                    % tv(i): mapped value function at k(i),z(i)

crit = 1.0;                         % init criterion for stopping iterations
tol = 1e-6;                         % tolerance 
it = 0;                             % number of iterations

while crit > tol
    for i=1:r
        k = kD_1000(i);
        c = [z(1)*k^alpha-kD_1000 z(2)*k^alpha-kD_1000]; 
        [tv(i,:) I] = max(log(c)+ beta*v*P');
        g(i,:) = [kD_1000(I(1)) kD_1000(I(2))];
    end
    crit = max(max(abs(tv-v)));
    v = tv;

end
kprime = [g(:,1);g(:,2)];
v;

% graph
figure1 = figure('Name','Value function n1000');
plot1 = plot(kD_1000,v);
set(plot1(1),'DisplayName','z1');
set(plot1(2),'Color',[1 0 0],'DisplayName','z2');
xlabel('Capital');
ylabel('Value of Capital');
title('Value function n1000');
legend('show');
set(legend,'Location','NorthWest');

figure2 = figure('Name','Policy function n1000');
plot2 = plot(kD_1000,g);
set(plot2(1),'DisplayName','z1');
set(plot2(2),'Color',[1 0 0],'DisplayName','z2');
xlabel('Capital');
ylabel('Next-period Capital');
title('Policy function n1000');
legend('show');
set(legend,'Location','NorthWest');

% ----- turn off diary ----- %
comptime = etime(clock, time0)
diary off