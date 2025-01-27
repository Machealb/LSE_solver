% verify the basis matrix form recurrences of Q_A and Q_L
clear, clc;
directory = pwd;
path(directory, path)
addpath(genpath('..'))
rng(2024);

%%------------- test matrices---------------
L = mmread('fxm2-6.mtx');
[~,n] = size(L);
A = 1*get_l(n,1);

% M = load('TF15.mat');
% A = M.Problem.A;
% [m,n] = size(A);
% L = 1*get_l(n,2);

% M = load('ch.mat');
% A = M.Problem.A;
% [m,n] = size(A);
% L = 1*get_l(n,1);

% M = load('t2d_q4.mat');
% AL = M.Problem.A;
% [m1,n] = size(AL);
% A = AL(1:8000, :);
% L = AL(8001:end, :);

type   = '1';
[x, x1, x2, b, d] = gen_prob2(A, L, type);
nx = norm(x1);


%%---------gLSQR ------------------------------------------------
type1 = 'semi';
k0  = 100;
tol = 0;
reorth = 1;
[X, res, bnd] = gLSQR(L, A, d, 1, k0, tol, reorth, type1);

k = size(res,1);
er = zeros(k, 1);  % relative error w.r.t. x
for i=1:k
    er(i) = norm(x1-X(:,i)) / nx;
end


%%-----------plot--------------------------------------
lw = 1.5; l = 1:k;

figure; 
semilogy(l,res,'->','Color','[0.8500 0.3250 0.0980]','MarkerIndices',1:5:k0,...
    'MarkerSize',5,'MarkerFaceColor','[0.8500 0.3250 0.0980]','LineWidth',1.5);
hold on;
semilogy(l,bnd,'-o','Color','[0 0.4470 0.7410]','MarkerIndices',1:5:k0,...
    'MarkerSize',6,'LineWidth',1.5);
legend('residual norm','estimation', 'Fontsize',14);
xlabel('Iteration','Fontsize',15);
ylabel('Relative  residual','Fontsize',15);
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);

figure; 
semilogy(l, er,'-^','Color','[0.4940 0.1840 0.5560]','MarkerIndices',1:5:k0,...
    'MarkerSize',5,'MarkerFaceColor','[0.4940 0.1840 0.5560]','LineWidth',1.5);
xlabel('Iteration','Fontsize',15);
ylabel('Relative  error','Fontsize',15);
% legend('$\mathtt{tol}=10^{-12}$','$\mathtt{tol}=10^{-10}$','$\mathtt{tol}=10^{-8}$','interpreter','latex','fontsize',18);
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);

figure; 
plot(1:n, x1,'-','Color','b','LineWidth',1.5);
hold on;
plot(1:n, X(:,end),'--','Color','m','LineWidth',1);
legend('True solution','Computed solution','fontsize',14);
% xlabel('Iteration','Fontsize',15);
% ylabel('Relative error','Fontsize',15);

