%% NOTICE FOR READERS
% This code is for paper
% "Data-driven Quantized Control of Partially
% Unknown Linear Systems with Noises"
% 
% Writen by Xingchen Li
% lixc21@mails.tsinghua.edu.cn
% Last modification at 2021-12-05
%
% You should have YALMIP with MOSEK


%% Data Initialization
clear
% System matrix
A=[-0.192, -0.936, -0.814;
      -0.918, +0.729, -0.724;
      -0.412, -0.135, -0.516];
B=[-0.554; 0.735; 0.528];
% dimension of state
n=size(B,1);
m=size(B,2);
% Caculate the eigenvalues of A
disp('The eigenvalues of A is:')
disp(eig(A)')
% Check for controllability
disp('The rank of Qk is:')
disp(rank([B A*B A*A*B]))

%%
c=0;
% while 1
% Data length
T=20;
% Second norm bound of noise: ||w||^2_2<=w_max
w_max=0.05;
% Generate data
[Xm,Xp,U,W]=gen_data(A,B,T,w_max);


%% Noise bound matrix
phi11=w_max*T*eye(n);
phi12=zeros(n,T);
phi22=-eye(T);


%% Theorem 14 in "From noisy data ... via a matrix S-lemma"
Y=sdpvar(n,n); % Matrix Y
X=sdpvar(1,n); % Matrix X
a=sdpvar(); % Scalars alpha
b=sdpvar(); % Scalars beta
d=sdpvar(); % Scalars delta^2

% Matrix defination
M=[Y-d*(B*B')-b*eye(n), zeros(n), B*X, zeros(n,1);
     zeros(n), zeros(n), Y, zeros(n,1);
     X'*B', Y', Y, X';
      zeros(1,n), zeros(1,n), X, 1];
N=[eye(n), Xp-B*U; 
    zeros(n,n), -Xm;
    zeros(n+1,n+T)]* ...
    [phi11 phi12; phi12' phi22]* ...
	[eye(n), Xp-B*U; 
    zeros(n,n), -Xm;
    zeros(n+1,n+T)]';

% SDP constraints
constraints=[Y>=0, a>=0, b>=0, [Y X'; X 1]>=0, M-a*N>=0];

% Options to YALMIP
options=sdpsettings('verbose',0,'solver','MOSEK');
% Optimize!
solution=optimize(constraints,-d,options);
% Get feasible solution
K=value(X)/value(Y);

disp(K)
disp(max(abs(eig(A+B*K))))
disp(value(d))


%% Function for generating data
function [Xm,Xp,U,W]=gen_data(A,B,T,w_max)
    % dimension of state
    n=size(B,1);
    m=size(B,2);
    % Data initialization for generation
    X=[randn(n,1) zeros(n,T)];
    U=zeros(m,T);
    W=zeros(n,T);
    
    % Uniform distribution noise in high-dimensional sphere
    w_num=0;
    while w_num<T
        % Uniform distribution noise in every dimension
        w=sqrt(w_max)*(2*rand(n,1)-1);
        % If in high-dimensional sphere
        if norm(w)^2<=w_max
            w_num=w_num+1;
            W(:,w_num)=w;
        end
    end
    
    % Get noisy data by recursion
    for i=1:T
        % Generate noisy data by gaussian
        U(:,i)=randn(m,1);
        % System recursion
        X(:,i+1)=A*X(:,i)+B*U(:,i)+W(:,i);
    end
    % Translate data for formality
    Xm=X(:,1:end-1); % X_minus
    Xp=X(:,2:end); % X_plus
end



