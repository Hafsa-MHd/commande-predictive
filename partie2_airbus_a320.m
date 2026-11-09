clear; clc; close all;

%% Partie 2 - MPC Airbus A320

A = [ 0.9996  0.0383  0.0131 -0.0322;
     -0.0056  0.9647  0.7446  0.0001;
      0.0020 -0.0097  0.9543  0;
      0.0001 -0.0005  0.0978  1.0000];

B = [ 0.0002  0.1002;
     -0.1230  0.0183;
     -0.2266  0.0586;
     -0.0114  0.0029];

C = [1  0  0  0;
     0 -1  0  7.74];

n  = size(A,1);
m  = size(B,2);
p  = size(C,1);
nz = p;
na = n + nz + m;

Ainc = [ A,            zeros(n,nz),  B;
         C,            eye(nz),      zeros(nz,m);
         zeros(m,n),   zeros(m,nz),  eye(m) ];

Binc = [ B;
         zeros(nz,m);
         eye(m) ];

Cy = [C,           zeros(p,nz),  zeros(p,m)];
Cz = [zeros(nz,n), eye(nz),      zeros(nz,m)];
Ex = [eye(n),      zeros(n,nz),  zeros(n,m)];
Eu = [zeros(m,n),  zeros(m,nz),  eye(m)];

T  = 20;
Qy = diag([10, 10]);
Qz = diag([1, 1]);
R  = diag([1, 1]);
P  = 10*Qy;

xmax = [12; 12; 1; pi/3];
xmin = -xmax;
umax = [pi/3; 3];
umin = -umax;

Nsim = 120;
Ts   = 0.1;

%% Matrices de prediction

Phi   = zeros(na*T, na);
Gamma = zeros(na*T, m*T);
Apow  = eye(na);
for i = 1:T
    Apow = Ainc*Apow;
    Phi((i-1)*na+1:i*na, :) = Apow;
    for j = 1:i
        Gamma((i-1)*na+1:i*na, (j-1)*m+1:j*m) = (Ainc^(i-j))*Binc;
    end
end

CyBig = kron(eye(T), Cy);
CzBig = kron(eye(T), Cz);
ExBig = kron(eye(T), Ex);
EuBig = kron(eye(T), Eu);

WyCell = repmat({Qy}, 1, T);
WyCell{T} = P;
Wy = blkdiag(WyCell{:});
Wz = kron(eye(T), Qz);
Wu = kron(eye(T), R);

My = CyBig*Gamma;
Mz = CzBig*Gamma;

Hqp = 2*(My'*Wy*My + Mz'*Wz*Mz + Wu);
Hqp = (Hqp + Hqp')/2;

Gx = ExBig*Gamma;
Gu = EuBig*Gamma;
nrx = size(Gx,1);
nru = size(Gu,1);
sCol = [ -ones(2*nrx,1); zeros(2*nru,1) ];
Aineq = [ [Gx; -Gx; Gu; -Gu], sCol ];
XmaxRep = repmat(xmax, T, 1);
XminRep = repmat(xmin, T, 1);
UmaxRep = repmat(umax, T, 1);
UminRep = repmat(umin, T, 1);

rhoQ = 1e4;
rhoL = 1e3;
Hqp = blkdiag(Hqp, 2*rhoQ);
lb  = [-inf(m*T,1); 0];

%% Scenario

yref = zeros(p, Nsim+T+1);
yref(1, :) = 5;
yref(2, 60:end) = 2;

dwind = zeros(n, Nsim+1);
dwind(1, 40:end) = 2.0;
dwind(2, 40:end) = -1.0;

%% Boucle MPC

x      = zeros(n,1);
z      = zeros(nz,1);
u_prev = zeros(m,1);

Xhist = zeros(n, Nsim+1);
Yhist = zeros(p, Nsim);
Uhist = zeros(m, Nsim);
Rhist = zeros(p, Nsim);
opts  = optimoptions('quadprog','Display','off');

for k = 1:Nsim
    xt0 = [x; z; u_prev];
    Rref = reshape(yref(:, k:k+T-1), [], 1);

    cy = CyBig*Phi*xt0;
    cz = CzBig*Phi*xt0;
    fqp = [ 2*( My'*Wy*(cy - Rref) + Mz'*Wz*cz ) ; rhoL ];

    bineq = [ XmaxRep - ExBig*Phi*xt0;
             -XminRep + ExBig*Phi*xt0;
              UmaxRep - EuBig*Phi*xt0;
             -UminRep + EuBig*Phi*xt0 ];

    [w, ~, flag] = quadprog(Hqp, fqp, Aineq, bineq, [], [], lb, [], [], opts);
    if flag <= 0
        warning('quadprog infaisable au pas k=%d', k);
        w = zeros(m*T+1,1);
    end
    dU = w(1:m*T);

    du = dU(1:m);
    u  = u_prev + du;
    u  = min(max(u, umin), umax);

    y = C*x;
    z = z + (y - yref(:,k));

    d = dwind(:,k);
    x = A*(x - d) + B*u + d;

    Yhist(:,k)   = y;
    Uhist(:,k)   = u;
    Rhist(:,k)   = yref(:,k);
    Xhist(:,k+1) = x;
    u_prev       = u;
end

%% Traces

tt = (0:Nsim-1)*Ts;

figure;
subplot(2,1,1);
plot(tt, Yhist(1,:),'b','LineWidth',1.5); hold on;
stairs(tt, Rhist(1,:),'--r','LineWidth',1.2); grid on;
ylabel('y_1'); legend('y_1','y_{ref,1}'); title('Poursuite de reference');
subplot(2,1,2);
plot(tt, Yhist(2,:),'b','LineWidth',1.5); hold on;
stairs(tt, Rhist(2,:),'--r','LineWidth',1.2); grid on;
ylabel('y_2'); xlabel('temps (s)'); legend('y_2','y_{ref,2}');

figure;
stairs(tt, Uhist','LineWidth',1.5); grid on; hold on;
yline(umax(1),'--k'); yline(umin(1),'--k');
yline(umax(2),'--m'); yline(umin(2),'--m');
xlabel('temps (s)'); ylabel('commandes');
legend('\delta','p','Location','best');
title('Commandes');

figure;
plot(tt, Xhist(:,1:Nsim)','LineWidth',1.3); grid on;
xlabel('temps (s)'); ylabel('etats');
legend('v_x','v_y','w','\theta','Location','best');
title('Etats');
