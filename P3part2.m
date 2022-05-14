%% Effect of Lipschitz constant estimation
% Consider one instance of Problem 2
% Change Lipchitz constant gradually and see how the convergence of the
% given algorithm goes

% set the initial values of the problem
n=100;
B=randn(n,n);
D=diag(log(10+exp(randn(n,1))));
Q=B*B'+D;
b=10*randn(n,1);
xstar=-Q\b;  % Optimal salution by first order optimality condition

% classical gradient descent algorithm
% set some fixed values
x0 = zeros(n,1); % initial point
L= max(eig(Q)); % Lipschitz constant: max eigenvalue of Q
Sigma = min(eig(Q)); % min eigenvalue of Q
K=L/Sigma; % condition number

% plot log-relative accuracy
% Constant---------------------------------------------------
alphaC = 1/L; % set an intial step-size
x=x0;
errCs = zeros(n,1);
for r=1:1000
    xCr = x-alphaC*(Q*x+b);
    x = xCr;
    erC=norm(xCr-xstar)/norm(x0-xstar);
    errCs(r)=log(erC);
end
epsilonC = log(norm(xCr-xstar)/norm(x0-xstar)); % log-relative accuracy of Constant method
figure(1);
plot(errCs);
title('epsilon of constant method versus iterations');
ylabel('log-relative accuracy');
xlabel('number of iterates');

% Accelerated gradient method---------------------------------------------------
%r=0;
x=x0;
xr=x0;
a=0;
ar=(1/2)*(1+sqrt(4*a^2+1));
t= (a-1)/ar;
errANs=zeros(n,1);
for r=1:1000
    erAN=norm(xr-xstar)/norm(x0-xstar);
    errANs(r)=log(erAN);
    arr=(1/2)*(1+sqrt(4*ar^2+1));
    tr=(ar-1)/arr;
    yrr=(1+tr)*xr-tr*x;
    xrr=yrr-(1/L)*(Q*(yrr)+b);
    a=ar;
    ar=arr;
    x=xr;
    xr=xrr;
end
epsilonAN = log(norm(xrr-xstar)/norm(x0-xstar));
figure(2);
plot(errANs);
title('epsilon of accelerate gradient method versus iterations');
ylabel('log-relative accuracy');
xlabel('number of iterates');

%  change the step-size gradually
% Constant---------------------------------------------------
AccuracyC = zeros(22,1); % generate an zero vector to store log-relative accuracy of the algorithm

for i=1:22
    x=x0;
    alphaCnew = (0.1*i)/L; % set an intial step-size
    for r=1:1000
        xCrnew = x-alphaCnew*(Q*x+b);
        x = xCrnew;
    end
    epsilonCnew = log(norm(xCrnew-xstar)/norm(x0-xstar)); % log-relative accuracy of Constant method
    AccuracyC(i) = epsilonCnew;
end
figure(3)
plot(AccuracyC);
title('epsilon of constant method versus step-size');
ylabel('log-relative accuracy');
xlabel('step-size');

% Accelerated gradient method---------------------------------------------------
%r=0;
AccuracyAN=zeros(22,1);
for i=1:22
    x=x0;
    xr=x0;
    a=0;
    ar=(1/2)*(1+sqrt(4*a^2+1));
    t= (a-1)/ar;
    for r=1:1000
        arr=(1/2)*(1+sqrt(4*ar^2+1));
        tr=(ar-1)/arr;
        yrr=(1+tr)*xr-tr*x;
        xrr=yrr-((0.1*i/L))*(Q*(yrr)+b);
        a=ar;
        ar=arr;
        x=xr;
        xr=xrr;
    end
    epsilonANnew = log(norm(xrr-xstar)/norm(x0-xstar));
    AccuracyAN(i)= epsilonANnew;
end
figure(4)
plot(AccuracyAN);
title('epsilon of accelerate gradient method versus step-size');
ylabel('log-relative accuracy');
xlabel('step-size');