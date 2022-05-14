% Jing Zheng 5-05-2017
% Triangular Arbitrage findinng
clear;clc;
% The input data
CC=load('12-31-2009.csv');
C=CC';
n=size(C,1);
lgC=zeros(n,n);
for i=1:n
    for j=1:n
        lgC(i,j)=log(C(i,j));
    end
end
% Create and solve the model
cvx_begin 
% decision variables
variable S(3,1);
variable T(3,1);
variables X1(4,1) X2(4,4) X3(4,1);
variables Y1(3,1) Y2(3,3) Y3(3,1);
variables Z1(2,1) Z2(2,2) Z3(2,1);

% objective function
sub1=lgC(1,2:5)*X1+lgC(2:5,1)'*X3+lgC(2,2:5)*(X2(1,:)')+lgC(3,2:5)*(X2(2,:)')+lgC(4,2:5)*(X2(3,:)')+lgC(5,2:5)*(X2(4,:)');
sub2=lgC(2,3:5)*Y1+lgC(3:5,2)'*Y3+lgC(3,3:5)*(Y2(1,:)')+lgC(4,3:5)*(Y2(2,:)')+lgC(5,3:5)*(Y2(3,:)');
sub3=lgC(3,4:5)*Z1+lgC(4:5,3)'*Z3+lgC(4,4:5)*(Z2(1,:)')+lgC(5,4:5)*(Z2(2,:)');
obj=sub1+sub2+sub3;

maximize (obj)

% constraints
sum(S)==1;
S(1,1)==sum(X1);
S(2,1)==sum(Y1);
S(3,1)==sum(Z1);
T(1,1)==sum(X3);
T(2,1)==sum(Y3);
T(3,1)==sum(Z3);
for i=1:4
    X1(i,1)==sum(X2(i,:));
    X3(i,1)==sum(X2(:,i));
    for j=1:4
        0<=X2(i,j)<=1;
    end
    0<=X1(i,1)<=1;
     0<=X3(i,1)<=1;
    X2(i,i)==0;
    
end
for i=1:3
    Y1(i,:)==sum(Y2(i,:));
    Y3(i,:)==sum(Y2(:,i));
    for j=1:3
        0<=Y2(i,j)<=1;
    end
    0<=Y1(i,1)<=1;
     0<=Y3(i,1)<=1;
    Y2(i,i)==0;
    0<=S(i,1)<=1;
    0<=T(i,1)<=1;
end
for i=1:2
    Z1(i,:)==sum(Z2(i,:));
    Z3(i,:)==sum(Z2(:,i));
    for j=1:2
        0<=Z2(i,j)<=1;
    end
    0<=Z1(i,1)<=1;
     0<=Z3(i,1)<=1;
    Z2(i,i)==0;
    
end
cvx_end

% output the result
if cvx_optval >=0
   disp(['The Triangular Arbitrage exists.'])
else
   disp(['The Triangular Arbitrage does not exists'])
end
