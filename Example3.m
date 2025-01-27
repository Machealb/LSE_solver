% verify the basis matrix form recurrences of Q_A and Q_L
clear, clc;
directory = pwd;
path(directory, path)
addpath(genpath('..'))
rng(2024);


%%------------- test matrices---------------
C = mmread('lp_bnl2.mtx');
[~,n] = size(C);
A = 1*get_l(n,1);

% AL = mmread('cage9.mtx');
% m = 2500;
% A = AL(1:m, :);
% C = AL(m+1:end, :);
% [~,n] = size(A);

type   = '2';
[x, b, d] = gen_prob1(A, C, type);
nx = norm(x);


%%---------gLSQR ------------------------------------------------
type1 = 'semi';
k0  = 50;
tol = 0;
reorth = 1;

k1 = 100;
k2 = 100;
k3 = 100;
[x11, X12, res11, res12] = KIDS2(A, C, b, d, 1e-10, 1e-8, k3);
[x21, X22, res21, res22] = KIDS2(A, C, b, d, 1e-8, 1e-6, k3);

kk1 = size(X12,2);
kk2 = size(X22,2);
% kk  = min(kk1,kk2);
X11 = kron(x11,ones(1,kk1));
X1  = X11+ X12;
kkk = size(X22,2);
X21 = kron(x21,ones(1,kk2));
X2 = X21 + X22;

er1 = zeros(kk1, 1); 
er2 = zeros(kk2, 1); 
for i=1:kk1
    er1(i) = norm(x-X1(:,i)) / nx;
end
for i=1:kk2
    er2(i) = norm(x-X2(:,i)) / nx;
end

x01 = lse_solv(A, C, b, d, '1');
x02 = lse_solv(A, C, b, d, '2');
er01 = norm(x-x01) / nx;
er02 = norm(x-x02) / nx;


%%-----------plot--------------------------------------
lw = 1.5;

% figure; 
% semilogy(l,res,'->','Color','[0.8500 0.3250 0.0980]','MarkerIndices',1:5:k0,...
%     'MarkerSize',5,'MarkerFaceColor','[0.8500 0.3250 0.0980]','LineWidth',1.5);
% hold on;
% semilogy(l,bnd+1e-3*tol*ones(k,1),'-o','Color','[0 0.4470 0.7410]','MarkerIndices',1:5:k0,...
%     'MarkerSize',6,'LineWidth',1.5);
% legend('residual norm','estimation', 'Fontsize',14);
% xlabel('Iteration','Fontsize',15);
% ylabel('Relative  residual','Fontsize',15);
% grid on;
% set(gca, 'GridAlpha', 0.2);
% set(gca, 'MinorGridAlpha', 0.01);

figure; 
semilogy(1:kk1, er1,'-o','Color','[0.9290 0.6940 0.1250]','MarkerIndices',1:5:kk1,...
    'MarkerSize',5,'MarkerFaceColor', '[0.9290 0.6940 0.1250]','LineWidth',1.5);
hold on;
semilogy(1:kk2, er2,'-d','Color','[0.4940 0.1840 0.5560]','MarkerIndices',1:5:k2,...
    'MarkerSize',5,'MarkerFaceColor','[0.4940 0.1840 0.5560]','LineWidth',1.5);
hold on;
semilogy(1:kk1, er01*ones(kk1,1),'k-','LineWidth',2.0);
hold on;
semilogy(1:kk1, er02*ones(kk1,1),'r--','LineWidth',2.0);
xlabel('Iteration','Fontsize',15);
ylabel('Relative  error','Fontsize',15);
legend('KIDS-II, $\tau_{1}=10^{-10}$','KIDS-II, $\tau_{1}=10^{-8}$','NS','DE','interpreter','latex','fontsize',18);
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);

figure; 
plot(1:n, x,'-','Color','b','LineWidth',1.5);
hold on;
plot(1:n,X2(:,end),'--','Color','m','LineWidth',1);
legend('True solution','Computed solution','fontsize',14);
% xlabel('Iteration','Fontsize',15);
% ylabel('Relative error','Fontsize',15);

