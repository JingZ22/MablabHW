%% Implement of ADMM
A=load('A.txt'); % Adjacency matrix of the graph
figure(1);
imagesc(A);
title('Adjacency matrix A');
% implement ADMM
n=100;
lam=5;
I = ones(n,1);
A_=lam*I*I'-A; 
%X=eye(n); % set the initial value 
Z=eye(n);
mu=ones(n,n);%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho=0.5;

for r=1:1000
    Xr=Z-rho^(-1)*mu-rho^(-1)*A_;
    [V D]=eig(Xr);
    for i=1:n
        if D(i,i)<0
            D(i,i)=0;
        end
    end
    Xr=V*D*inv(V);
    Zr=Xr+rho^(-1)*mu;
    for j=1:n
        Zr(j,j)=1;
    end
    mur=mu+rho*(Xr-Zr);
    mu=mur;
    Z=Zr;
end
% imagesc(Xr);
% rank one approximation
[B C]=eig(Xr);
eigva=diag(C);
[alpha,in]=max(eigva);
x=B(:,in);
X_=alpha*x*x';
% partition
% get the index of two partitions
Index=zeros(n,1);
ii=1;
jj=1;
for k=1:n
    if x(k,1)>0
        Index(ii,1)=k;
        ii=ii+1;
    else
        Index(n-jj+1,1)=k;
        jj=jj+1;
    end
end
nA=zeros(n,n);

for t1=1:n
    i=Index(t1,1);
    nA(:,t1)=X_(:,i);
end
nnA=zeros(n,n);
for t2=1:n
    j=Index(t2,1);
    nnA(t2,:)=nA(j,:);
end
figure(2);
imagesc(nnA);
    





