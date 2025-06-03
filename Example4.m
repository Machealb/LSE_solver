% test for scaling
clear, clc;
directory = pwd;
path(directory, path)
addpath(genpath('..'))
rng(2024);


%%------------- test matrices---------------
n= 8000;
A_type   = '3';
C_type   = '1';
x_type   = '1';
[A, C, x, b, d, x1, x2] = gen_prob3(n);
nx = norm(x);


%%---------gLSQR ------------------------------------------------
type1 = 'semi';
k0  = 50;
tol = 1e-12;
reorth = 1;

k1 = 80;
k2 = 80;
k3 = 80;
tic;
[X11, X12, res11, res12] = KIDS1(A, C, b, d, tol, tol, k1, k2, type1);
toc;
time1 = toc-tic;
tic;
[x21, X22, res21, res22] = KIDS2(A, C, b, d, 0.01*tol, 0.1*tol, k3);
toc;
time2 = toc-tic;

kk1 = size(X11,2);
kk2 = size(X12,2);
kk  = min(kk1,kk1);
X1  = X11(:,1:kk) + X12(:,1:kk);
kkk = size(X22,2);
X21 = kron(x21,ones(1,k3));
X2 = X21 + X22;

er1 = zeros(kk, 1); 
er2 = zeros(kkk, 1); 
for i=1:kk
    er1(i) = norm(x-X1(:,i)) / nx;
end
for i=1:kkk
    er2(i) = norm(x-X2(:,i)) / nx;
end

% tic;
% x01 = lse_solv(A, C, b, d, '1');
% toc;
% tic;
% x02 = lse_solv(A, C, b, d, '2');
% toc;
% er01 = norm(x-x01) / nx;
% er02 = norm(x-x02) / nx;


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
semilogy(1:kk, er1,'-o','Color','[0.3010 0.7450 0.9330]','MarkerIndices',1:9:kk,...
    'MarkerSize',5,'MarkerFaceColor', '[0.3010 0.7450 0.9330]','LineWidth',1.5);
hold on;
semilogy(1:kkk, er2,'-d','Color','[0.4660 0.6740 0.1880]','MarkerIndices',1:9:kkk,...
    'MarkerSize',5,'MarkerFaceColor','[0.4660 0.6740 0.1880]','LineWidth',1.5);
hold on;
semilogy(1:kk, er01*ones(kk,1),'k-','LineWidth',2.0);
hold on;
semilogy(1:kk, er02*ones(kk,1),'r--','LineWidth',2.0);
xlabel('Iteration','Fontsize',15);
ylabel('Relative  error','Fontsize',15);
legend('KIDS-I','KIDS-II','NS','DE','interpreter','latex','fontsize',18);
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

