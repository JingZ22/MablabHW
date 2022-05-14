%% Implement PGM

%downlaod sampling data
X = load('X.txt'); % matrix of samples x
Y = load('Y.txt'); % vector of label y

%set some intial values
n = 2000; % number of simples
d = 1000; % nunber of features
n_tr = 200; % number of samples for training
n_test = 500; % number of samples for testing

X_tr = X(1:n_tr, 1:d); % samples for training
Y_tr = Y(1:n_tr, 1);
X_test = X(n-n_test+1:n, 1:d); % samples for test
Y_test = Y(n-n_test+1:n, 1);

alpha = 0.001; % initial step-size
lambda = 0.85; % parameter
wr=zeros(d,1);
for r=1:1000
    
    % compute the gradient of f1(wr)
    ff1=zeros(n_tr,d);
    for j=1:n_tr
        ff1(j,:) = X_tr(j,:)*(exp(X_tr(j,:)*wr)/(1+exp(X_tr(j,:)*wr))-Y_tr(j));
        %aaaa=exp(X_tr(j,:)*wr)/(1+exp(X_tr(j,:)*wr))-Y_tr(j);
    end
    gf11=mean(ff1);
    gf1=gf11';
    % soft-thresholding function
    
    Ar=wr-alpha*gf1;
    for j=1:d
        if Ar(j)>lambda*alpha
            wr(j)=Ar(j)-lambda*alpha;
        elseif Ar(j)<-lambda*alpha
            wr(j)=Ar(j)+lambda*alpha;
        else
            wr(j)=0;
        end
    end
end





