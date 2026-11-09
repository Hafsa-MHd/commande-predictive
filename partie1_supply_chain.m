clear; clc; close all;

%% Partie 1 - MPC supply chain

n = 5;
m = 9;

Ain  = zeros(n,m);
Aout = zeros(n,m);

Ain(1,1) = 1;
Aout(1,2) = 1;  Ain(2,2) = 1;
Aout(1,3) = 1;  Ain(3,3) = 1;
Aout(2,4) = 1;  Ain(4,4) = 1;
Aout(3,5) = 1;  Ain(4,5) = 1;
Aout(3,6) = 1;  Ain(5,6) = 1;
Aout(4,7) = 1;  Ain(5,7) = 1;
Aout(4,8) = 1;
Aout(5,9) = 1;

B = Ain - Aout;

xmax = 20;
umax = 8;
x0   = [15; 5; 8; 2; 0];
c    = [ones(7,1); -1; -1];

T    = 10;
Nsim = 30;

%% Matrices QP

nz = m*T + n*T;
idxU = @(k) (k*m + 1):(k*m + m);
idxX = @(k) (m*T + (k-1)*n + 1):(m*T + (k-1)*n + n);

H = zeros(nz);
f = zeros(nz,1);
for k = 1:T
    H(idxX(k), idxX(k)) = 2*eye(n);
    f(idxX(k)) = ones(n,1);
end
for k = 0:T-1
    f(idxU(k)) = c;
end
H = (H + H')/2;

Aeq = zeros(n*T, nz);
beq = zeros(n*T, 1);
for k = 0:T-1
    rows = (k*n + 1):(k*n + n);
    Aeq(rows, idxU(k))   = -B;
    Aeq(rows, idxX(k+1)) =  eye(n);
    if k >= 1
        Aeq(rows, idxX(k)) = -eye(n);
    else
        beq(rows) = x0;
    end
end

Aineq = zeros(n*T, nz);
bineq = zeros(n*T, 1);
for k = 0:T-1
    rows = (k*n + 1):(k*n + n);
    Aineq(rows, idxU(k)) = Aout;
    if k >= 1
        Aineq(rows, idxX(k)) = -eye(n);
    else
        bineq(rows) = x0;
    end
end

lb = zeros(nz,1);
ub = [umax*ones(m*T,1); xmax*ones(n*T,1)];

%% Boucle MPC

Xhist = zeros(n, Nsim+1);
Xhist(:,1) = x0;
Uhist = zeros(m, Nsim);
Jhist = zeros(1, Nsim);
opts  = optimoptions('quadprog','Display','off');

x = x0;
for t = 1:Nsim
    beq(1:n)   = x;
    bineq(1:n) = x;

    [z, Jval, flag] = quadprog(H, f, Aineq, bineq, Aeq, beq, lb, ub, [], opts);
    if flag <= 0
        warning('quadprog : pas de solution au pas t=%d', t);
        break;
    end

    u = z(idxU(0));
    x = x + B*u;

    Uhist(:,t)   = u;
    Xhist(:,t+1) = x;
    Jhist(t)     = Jval;
end

%% Traces

figure;
stairs(0:Nsim, Xhist','LineWidth',1.5); grid on;
yline(xmax,'--k');
xlabel('temps t'); ylabel('stock x_i(t)');
title('Evolution des stocks (MPC)');
legend(arrayfun(@(i) sprintf('N%d',i),1:n,'UniformOutput',false),'Location','best');

figure;
stairs(0:Nsim-1, Uhist','LineWidth',1.2); grid on;
yline(umax,'--k');
xlabel('temps t'); ylabel('flux u_j(t)');
title('Flux de transport');
legend(arrayfun(@(j) sprintf('L%d',j),1:m,'UniformOutput',false),'Location','best');

figure;
plot(1:Nsim, Jhist,'-o','LineWidth',1.5); grid on;
xlabel('temps t'); ylabel('J^*(t)');
title('Cout optimal');
