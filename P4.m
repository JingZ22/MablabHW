%% Effect of restart in accelerated Nesterov
% Consider one instance of problem 2

% form a quadratic optimization problem
n=100;
B=randn(n,n);
D=diag(log(10+exp(randn(n,1))));
Q=B*B'+D;
b=10*randn(n,1);
%errs = zeros(n,1);
% Optimal salution by first order optimality condition
xstar=-Q\b; 
% Accelerated gradient method---------------------------------------------------
% set some fixed values
x0 = zeros(n,1); % initial point
L= max(eig(Q)); % Lipschitz constant: max eigenvalue of Q
Sigma = min(eig(Q)); % min eigenvalue of Q
K=L/Sigma; % condition number



% Accelerated gradient method with no restart
x=x0;
xr=x0;
a=0;
ar=(1/2)*(1+sqrt(4*a^2+1));
t= (a-1)/ar;
errANs1=zeros(n,1);
for r=1:1000
    erAN=norm(xr-xstar)/norm(x0-xstar);
    errANs1(r)=log(erAN);
    arr=(1/2)*(1+sqrt(4*ar^2+1));
    tr=(ar-1)/arr;
    %ar=(1/2)*(1+sqrt(4*a^2+1));
    yrr=(1+tr)*xr-tr*x;
    xrr=yrr-(1/L)*(Q*(yrr)+b);
    a=ar;
    ar=arr;
    x=xr;
    xr=xrr;
end
epsilonAN1 = log(norm(xrr-xstar)/norm(x0-xstar));
subplot(2,2,1);
plot(errANs1);
title('epsilon of accelarate Nesterov with no restart');
ylabel('log-relative accuracy');
xlabel('number of iterates');

% Accelerated gradient method with first T restart
T=floor(sqrt(K)); % restart at every T iterations
% set some initial values
xr=x0;
a=0;
ar=(1/2)*(1+sqrt(4*a^2+1));
t= (a-1)/ar;
errANs2=zeros(n,1);
for r=1:1000
    if (mod(r-1,T)==0)
        a=0;
        ar=(1/2)*(1+sqrt(4*a^2+1));
        t= (a-1)/ar;
    end   
    erAN=norm(xr-xstar)/norm(x0-xstar);
    errANs2(r)=log(erAN);
    arr=(1/2)*(1+sqrt(4*ar^2+1));
    tr=(ar-1)/arr;
    %ar=(1/2)*(1+sqrt(4*a^2+1));
    yrr=(1+tr)*xr-tr*x;
    xrr=yrr-(1/L)*(Q*(yrr)+b);
    a=ar;
    ar=arr;
    x=xr;
    xr=xrr;
end
epsilonAN2 = log(norm(xrr-xstar)/norm(x0-xstar));
subplot(2,2,2);
plot(errANs2);
title('epsilon of accelarate Nesterov with first T restart');
ylabel('log-relative accuracy');
xlabel('number of iterates');




% Accelerated gradient method with second T restart
T=floor(5*sqrt(K)); % restart at every T iterations
% set some initial values
xr=x0;
a=0;
ar=(1/2)*(1+sqrt(4*a^2+1));
t= (a-1)/ar;
errANs3=zeros(n,1);
for r=1:1000
    if (mod(r-1,T)==0)
        a=0;
        ar=(1/2)*(1+sqrt(4*a^2+1));
        t= (a-1)/ar;
    end   
    erAN=norm(xr-xstar)/norm(x0-xstar);
    errANs3(r)=log(erAN);
    arr=(1/2)*(1+sqrt(4*ar^2+1));
    tr=(ar-1)/arr;
    %ar=(1/2)*(1+sqrt(4*a^2+1));
    yrr=(1+tr)*xr-tr*x;
    xrr=yrr-(1/L)*(Q*(yrr)+b);
    a=ar;
    ar=arr;
    x=xr;
    xr=xrr;
end
epsilonAN3 = log(norm(xrr-xstar)/norm(x0-xstar));
subplot(2,2,3);
plot(errANs3);
title('epsilon of accelarate Nesterov with second T restart');
ylabel('log-relative accuracy');
xlabel('number of iterates');



% Accelerated gradient method with third T restart
T=floor(20*sqrt(K)); % restart at every T iterations
% set some initial values
xr=x0;
a=0;
ar=(1/2)*(1+sqrt(4*a^2+1));
t= (a-1)/ar;
errANs4=zeros(n,1);
for r=1:1000
    if (mod(r-1,T)==0)
        a=0;
        ar=(1/2)*(1+sqrt(4*a^2+1));
        t= (a-1)/ar;
    end   
    erAN=norm(xr-xstar)/norm(x0-xstar);
    errANs4(r)=log(erAN);
    arr=(1/2)*(1+sqrt(4*ar^2+1));
    tr=(ar-1)/arr;
    %ar=(1/2)*(1+sqrt(4*a^2+1));
    yrr=(1+tr)*xr-tr*x;
    xrr=yrr-(1/L)*(Q*(yrr)+b);
    a=ar;
    ar=arr;
    x=xr;
    xr=xrr;
end
epsilonAN4 = log(norm(xrr-xstar)/norm(x0-xstar));
subplot(2,2,4);
plot(errANs4);
title('epsilon of accelarate Nesterov with third T restart');
ylabel('log-relative accuracy');
xlabel('number of iterates');











