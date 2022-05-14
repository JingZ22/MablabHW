%clc;clear;
%% form a quadratic optimization problem and make five realizaitons for the problem 2

Epsilon = zeros(6,5); % set a zero matrix to store epsilon of six algorithms over five realizations
Kappa = zeros(5,1); % set a zero vector to store conditon numbers of five realizations
for k=1:5
    
    % form a quadratic optimization problem
    n=100;
    B=randn(n,n);
    D=diag(log(10+exp(randn(n,1))));
    Q=B*B'+D;
    b=10*randn(n,1);
    % Optimal salution by first order optimality condition
    xstar=-Q\b;
    % classical gradient descent algorithm
    % set some fixed values
    x0 = zeros(n,1); % initial point
    L= max(eig(Q)); % Lipschitz constant: max eigenvalue of Q
    Sigma = min(eig(Q)); % min eigenvalue of Q
    K=L/Sigma; % condition number
    Kappa(k,1)=K;
    
    % Exact line search--------------------------------
    x=x0;
    for r=1:1000
        alphaLr=(Q*x+b)'*(Q*x+b)/((Q*x+b)'*Q*(Q*x+b));
        xLr= x-alphaLr*(Q*x+b);
        x=xLr;
    end
    epsilonL = log(norm(xLr-xstar)/norm(x0-xstar)); % log-relative accuracy of Diminishing method
    Epsilon(1,k)=epsilonL;
    
    % Armijo rule----------------------------------------
    beta=0.5; % set an initial parameter beta
    sigma=0.5; % set an initial parameter sigma
    alphaA=1; % set an intial step-size
    x=x0;
    flag = 1;
    for r=1:1000
        i=1;
        while(flag)      
            if ((1/2)*x'*Q*x+b'*x)-((1/2)*(x-alphaA*beta^i*(Q*x+b))'*Q*(x-alphaA*beta^i*(Q*x+b))+b'*(x-alphaA*beta^i*(Q*x+b)))>=sigma*alphaA*beta^i*(Q*x+b)'*(Q*x+b);
                break;
            end
            i = i+1;
        end
        xAr = x-alphaA*beta^i*(Q*x+b);
        x=xAr;
    end
    epsilonA = log(norm(xAr-xstar)/norm(x0-xstar)); % log-relative accuracy of Armijo method
     Epsilon(2,k)=epsilonA;
    
     % Diminishing---------------------------------------------------
    alphaD = 0.1; % set an intial step-size
    x=x0;
    for r=1:1000
        alphaDr = alphaD/r;
        xDr = x-alphaDr*(Q*x+b);
        x = xDr;
    end
    epsilonD = log(norm(xDr-xstar)/norm(x0-xstar)); % log-relative accuracy of Diminishing method
     Epsilon(3,k)=epsilonD;
    
     % Constant---------------------------------------------------
    alphaC = 1/L; % set an intial step-size
    x=x0;
    for r=1:1000
        xCr = x-alphaC*(Q*x+b);
        x = xCr;
    end
    epsilonC = log(norm(xCr-xstar)/norm(x0-xstar)); % log-relative accuracy of Constant method
     Epsilon(4,k)=epsilonC;
    
    % Accelerated gradient method---------------------------------------------------
    x=x0;
    xr=x0;
    a=0;
    ar=(1/2)*(1+sqrt(4*a^2+1));
    t= (a-1)/ar;
    for r=1:1000
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
    Epsilon(5,k)=epsilonAN;
     
    % Newton's method---------------------------------------------------
    x=x0;
    for r=1:1000
        xNr=x-Q\(Q*x+b);
        x=xNr;
    end
    epsilonN = log(norm(xNr-xstar)/norm(x0-xstar)); % log-relative accuracy of Newton's method
     Epsilon(6,k)=epsilonN;
end
% compute the average epsilon value of five realizations
aveEpsilon = mean(Epsilon,2)















