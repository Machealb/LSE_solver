% test for function handle matrix
clear, clc;
directory = pwd;
path(directory, path)
addpath(genpath('..'))
rng(2024);


%%------------- test matrices---------------
n= 10000;
t = 500;
x_type   = '4';
[A, C, x, b, d, x1, x2] = gen_prob4(n, t, x_type);
nx = norm(x);
nx1 = norm(x1);
nx2 = norm(x2);


%%---------gLSQR ------------------------------------------------
type1 = 'semi';
k0  = 50;
tol = 1e-8;
reorth = 0;

k1 = 20;
k2 = 20;
k3 = 20;
tic;
[X11, X12, res11, res12] = KIDS1(A, C, b, d, tol, tol, k1, k2, type1);
toc;
tic;
[x21, X22, res21, res22] = KIDS2(A, C, b, d, 1*tol, 1*tol, k3);
toc;

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



%%-----------plot--------------------------------------
lw = 1.5;


figure; 
semilogy(1:kk, er1,'-o','Color','[0.3010 0.7450 0.9330]','MarkerIndices',1:1:kk,...
    'MarkerSize',5,'MarkerFaceColor', '[0.3010 0.7450 0.9330]','LineWidth',1.5);
hold on;
semilogy(1:kkk, 5*er2,'-d','Color','[0.4660 0.6740 0.1880]','MarkerIndices',1:1:kkk,...
    'MarkerSize',5,'MarkerFaceColor','[0.4660 0.6740 0.1880]','LineWidth',1.5);
xlabel('Iteration','Fontsize',15);
ylabel('Relative  error','Fontsize',15);
legend('KIDS-I','KIDS-II','interpreter','latex','fontsize',18);
grid on;
set(gca, 'GridAlpha', 0.2);
set(gca, 'MinorGridAlpha', 0.01);

figure; 
plot(1:n, x,'-','Color','b','LineWidth',1.5);
hold on;
plot(1:n,X1(:,end),'--','Color','m','LineWidth',1);
legend('True solution','Computed solution','fontsize',14);

