% -----------------------------------------------------------------------------
% File: main.m
% Author: Antoine Leeman (aleeman@ethz.ch)
% Date: 15th May 2023
% License: MIT
% Reference:
%{
@inproceedings{leeman2023a,
  title = {Predictive safety filter using system level synthesis},
  year = {2023},
  booktitle = {Proceedings of the 5th Annual Learning for Dynamics and Control Conference, PMLR},
  volume = {211},
  author={Leeman, Antoine P. and K{\"o}hler, Johannes and Benanni, Samir and Zeilinger, Melanie N.},
  pages = {1180-1192},
  doi = {10.3929/ethz-b-000615512}
}
%}
% Link: https://arxiv.org/abs/2212.02111
% -----------------------------------------------------------------------------
%%
clear all;
close all;
yalmip('clear');
clc;
%% Initialisation - System
m = Integrator();  % create an instance of the Integrator class
nx = m.nx;         % get the number of states
nu = m.nu;         % get the number of inputs
Bw = m.Bw;          % get the disturbance matrix from the Integrator class
W = Bw*Polyhedron([eye(nx);-eye(nx)],ones(2*nx,1)); % create a polyhedron representing the disturbance set
A = m.A;
B = m.B;

X = m.x_max*Polyhedron([eye(nx);-eye(nx)],ones(2*nx,1)); % create a polyhedron representing the state constraints
U = m.u_max*Polyhedron([eye(nu);-eye(nu)],ones(2*nu,1)); % create a polyhedron representing the input constraints
N = m.N;

Kf = -dlqr(A,B,m.Q_cost,m.R_cost); % compute the feedback gain matrix
A_cl = A+B*Kf; % compute the closed-loop system matrix

% compute mRPI
E = W; % initialize E to the disturbance set
k=1;
while k <= 25 % repeat the following loop up to 25 times
    k = k + 1;
    E_prev = E;
    E = E + A_cl^k * W; % compute the set E using the closed-loop system matrix and the disturbance set
    E = E.minVRep(); % compute the minimal vertex representation of E
    if eq(E, E_prev) % break the loop if E has converged
        break
    end
end

mRPI = E;
X_E = minus(X,E); % compute the set difference between the state constraints and E
U_KE = minus(U,Kf*E); % compute the set difference between the input constraints and Kf*E

% compute max PI
system = LTISystem('A',A_cl); % create an instance of the LTISystem class with the closed-loop system matrix
Xp = Polyhedron('A',[X_E.H(:,1:nx); U_KE.H(:,1:m.nu)*Kf], 'b', [X_E.H(:,nx+1); U_KE.H(:,1+nu);]);
system.x.with('setConstraint');
system.x.setConstraint = Xp; % set the state constraints to Xp
max_PI = system.invariantSet(); % compute the maximal positively invariant set
terminal_set_MPSC = max_PI; % set the terminal set to the maximal positively invariant set

system = ULTISystem('A',A_cl,'D',W); % create an instance of the ULTISystem class with the closed-loop system matrix and the disturbance set
Xp = Polyhedron('A',[X.H(:,1:nx); U.H(:,1:nu)*Kf], 'b', [X.H(:,nx+1); U.H(:,1+nu);]);
system.x.with('setConstraint');
system.x.setConstraint = Xp; % set the state constraints to Xp
system.d.max = [Bw(1,1); Bw(2,2)]; % set the maximum disturbance values
system.d.min = -[Bw(1,1); Bw(2,2)]; % set the minimum disturbance values
max_RPI = system.invariantSet(); % compute the maximal robust positively invariant

system = ULTISystem('A',A,'B',B,'D',W); % create an instance of the ULTISystem class with the system matrix and the disturbance set
Xp = Polyhedron('A',X.H(:,1:nx), 'b', X.H(:,nx+1)); % create a polyhedron representing the state constraints
Up = Polyhedron('A',U.H(:,1:nu), 'b',  U.H(:,1+nu)); % create a polyhedron representing the input constraints
system.x.with('setConstraint');
system.x.setConstraint = Xp; % set the state constraints to Xp
system.u.with('setConstraint');
system.u.setConstraint = Up; % set the input constraints to Up
system.d.max = [Bw(1,1); Bw(2,2)]; % set the maximum disturbance values
system.d.min = -[Bw(1,1); Bw(2,2)]; % set the minimum disturbance values
max_RCI = system.invariantSet(); % compute the maximal robust control invariant set

%% Initialisation - System Level Model Predictive Safety Controller 
N=5;
% Define decision variables
Z = sdpvar(nx, N + 1, 'full'); % State trajectory variables
V = sdpvar(nu, N, 'full');     % Input trajectory variables
X0 = sdpvar(nx, 1, 'full');    % Initial state variable
U_L = sdpvar(nu, 1, 'full');   % learned-input bound variable

Phi_x = sdpvar( (N + 1) * nx, (N + 1) * nx, 'full');
Phi_u = sdpvar( (N + 1) * nu, (N + 1) * nx, 'full');

objective =(V(:,1) - U_L)^2;

% Construct the sigma matrix
sigma_seq = kron(eye(N), m.Bw);
Sigma_mat = blkdiag(eye(nx),sigma_seq);

% % Define the objective function
% objective = Z(:,N+1)'*m.Q_cost*Z(:,N+1);
% for k=1:N
%     objective = objective + Z(:,k)'*m.Q_cost*Z(:,k) + V(:,k)'*m.R_cost*V(:,k);
% end

%objective = objective + norm([kron(eye(N+1),m.Q_cost)* Phi_x;kron(eye(N+1),m.Q_cost)* Phi_u],'fro')^2;
%objective = objective + norm([Phi_x;Phi_u],'fro')^2;


% Initialise the constraints
% constraints = X0 == [0;0];
constraints = [];
% Add structure constraints for Phi_x and Phi_u
for k = 1 : N
    constraints = [constraints, Phi_x( (k - 1)*nx + 1: k*nx, k*nx + 1: end) == zeros(nx, (N + 1 - k)*nx)];
end

for k = 1 : N
    constraints = [constraints, Phi_u( (k - 1)*nu + 1: k*nu, k*nx + 1: end) == zeros(nu, (N + 1 - k)*nx)];
end

% Define block downshift operator
Z_block = kron(diag(ones(1,N),-1), eye(nx));
ZA_block = Z_block*blkdiag(kron(eye(N), A), zeros(nx, nx));
ZB_block = Z_block*blkdiag(kron(eye(N), B), zeros(nx, nu));
Id = eye((N + 1)*nx);

% Add System Level Parametrisation constraint
constraints = [constraints, (Id - ZA_block)*Phi_x - ZB_block*Phi_u == Sigma_mat];

% Add initial state constraint
constraints = [ constraints, Z(:,1)==X0 ];

% Add state dynamics constraints
for k=1:N
    constraints = [ constraints, Z(:,k+1)==A*Z(:,k)+B*V(:,k)];
end
% 
% state constraints
Fx = m.F_x;
bx = m.b_x;
nFx = length(bx);
for ii = 1:N
    for jj = 1: nFx
        f = Fx(jj,:); b = bx(jj);
        LHS = f*Z(:,ii);
        for kk = 1:ii-1
            LHS = LHS + norm(f*Phi_x((ii-1)*nx+1:ii*nx,kk*nx+1:(kk+1)*nx), 2);
        end
        constraints = [constraints, LHS <= b];
    end
end

% % terminal constraint   
% X_f = max_RPI.H;
% Ft= X_f(:,1:nx);
% bt = X_f(:,nx+1);
% nFt = length(bt);
% for jj = 1:nFt
%     f = Ft(jj,:); b = bt(jj);
%     LHS = f*Z(:,N+1);
%     for kk = 1:N
%         LHS = LHS + norm(f*Phi_x(N*nx+1:(N+1)*nx,kk*nx+1:(kk+1)*nx),2);
%     end
%     constraints = [constraints, LHS <= b];
% end

% add input constraint
Fu = m.F_u;
bu = m.b_u;
nFu = length(bu);
for ii = 1:N
    for jj = 1: nFu
        f = Fu(jj,:); b = bu(jj);
        LHS = f*V(:,ii);
        for kk = 1:ii-1
            LHS = LHS + norm(f*Phi_u((ii-1)*nu+1:ii*nu,kk*nx+1:(kk+1)*nx),2);
        end
        constraints = [constraints, LHS <= b];
    end
end

%options = sdpsettings('verbose',1);
%optimize(constraints,objective,options);

options = sdpsettings('verbose',1,'solver','sedumi');
sol_SL_MPSF = optimizer(constraints,objective,options,[X0;U_L],V(1));


Z = value(Z);
V = value(V);

plot(Z(1,:), Z(2,:))
%%

ii=3;
kk=1; a = value(Phi_u((ii-1)*nu+1:ii*nu,kk*nx+1:(kk+1)*nx))*inv(value(Phi_x((ii-1)*nx+1:ii*nx,kk*nx+1:(kk+1)*nx)))
kk=2; b=  value(Phi_u((ii-1)*nu+1:ii*nu,kk*nx+1:(kk+1)*nx))*inv(value(Phi_x((ii-1)*nx+1:ii*nx,kk*nx+1:(kk+1)*nx)))
a-b


ii=9;
kk=7; a = value(Phi_u((ii-1)*nu+1:ii*nu,kk*nx+1:(kk+1)*nx))*inv(value(Phi_x((ii-1)*nx+1:ii*nx,kk*nx+1:(kk+1)*nx)));
kk=2; b=  value(Phi_u((ii-1)*nu+1:ii*nu,kk*nx+1:(kk+1)*nx))*inv(value(Phi_x((ii-1)*nx+1:ii*nx,kk*nx+1:(kk+1)*nx)));

a-b


