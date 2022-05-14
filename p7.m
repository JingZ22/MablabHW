%% Gradient Projection
% set some initial values
n=100;
B=randn(100,100);
D=diag(exp(randn(n,1)));
Q=B*B'+D;
b=10*randn(n,1);
x0=zeros(n,1); % intial point
L= max(eig(Q)); % Lipschitz constant: max eigenvalue of Q
iters=2000;% number of iterations
% projection to a feasible point
% constant step-size
alphaC = 1/L; 

logThetoListC=zeros(iters,1);
x=x0;
yCr= zeros(n,1);
for r=1:iters
    xCr = x-alphaC*(Q*x+b);
    for i = 1:n
        if xCr(i) > 1;
            xCr(i) = 1;
        elseif xCr(i) <0;
            xCr(i) = 0;
        end
    end
    yCr = xCr-(Q*xCr+b);
    for j = 1:n
        if yCr(j) > 1;
            yCr(j) = 1;
        elseif yCr(j) <0;
            yCr(j) = 0;
        end
    end
    logThetaC= log((norm(xCr-yCr))^2);
    x=xCr;
    logThetoListC(r)=logThetaC;

end
figure(1);
plot(logThetoListC);
title('log theta of constant step-size versus iterations');
ylabel('log-theta error');
xlabel('iterations');

% diminishing step-size
logThetaListD=zeros(iters,1);
yDr=zeros(n,1);
x=x0;
for r=1:iters
    alphaDr = 5*log(r)/(r*L);
    xDr = x-alphaDr*(Q*x+b);
   
    for i = 1:n
        if xDr(i) > 1;
            xDr(i) = 1;
        elseif xDr(i) <0;
            xDr(i) = 0;
        end
    end
    yDr = xDr-(Q*xDr+b);
    for j = 1:n
        if yDr(j) > 1;
            yDr(j) = 1;
        elseif yDr(j) <0;
            yDr(j) = 0;
        end
    end
    logThetaD= log((norm(xDr-yDr))^2);
    logThetaListD(r)=logThetaD;
    x=xDr;
end

figure(2);
plot(logThetaListD);
title('log theta of diminishing step-size versus iterations');
ylabel('log-theta error');
xlabel('iterations');









