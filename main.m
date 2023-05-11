% -----------------------------------------------------------------------------
% File: sls_safety_filter_leeman2023.m
% Author: Antoine Leeman (aleeman@ethz.ch)
% Date: 15th May 2023
% License: MIT
% Reference:
% A.P. Leeman, J. Köhler, S. Benanni, M.N. Zeilinger, "Predictive Safety Filter
% Using System Level Synthesis", 2023.
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
Bw = m.E;          % get the disturbance matrix from the Integrator class
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

%% Model Predictive Safety Controller 
% Define optimisation variables
Z = sdpvar(nx, N + 1, 'full'); % State trajectory variables
V = sdpvar(nu, N, 'full');     % Input trajectory variables
X0 = sdpvar(nx, 1, 'full');    % Initial state variable
U_L = sdpvar(nu, 1, 'full');   % learned-input bound variable

% Define the objective function
objective =(V(:,1) - U_L)^2;

% Initialise the constraints
constraints=[];

% Initial state constraint
constraints = [constraints, mRPI.H(:,1:2)*(X0-Z(:,1)) <= mRPI.H(:,3)];

% Loop over the horizon length
for k=1:N
    constraints = [ constraints, U_KE.H(:,1)*V(:,k) <= U_KE.H(:,2)];% Input constraint
    constraints = [ constraints, X_E.H(:,1:2)*Z(:,k) <= X_E.H(:,3)];% State constraint
    constraints = [ constraints, Z(:,k+1)==A*Z(:,k)+B*V(:,k)]; % System dynamics constraint
end

% Terminal state constraint
constraints = [constraints, terminal_set_MPSC.H(:,1:2)*(Z(:,N+1)) <= terminal_set_MPSC.H(:,3)];

% Set optimization options
options = sdpsettings('verbose',1,'solver','mosek');

% Define the optimizer
sol = optimizer(constraints,objective,options,[X0;U_L],V(1));

yalmip('clear');
clear Z V X0 U_L objective constraints;
%% Initialisation - System Level Model Predictive Safety Controller 

% Define decision variables
Z = sdpvar(nx, N + 1, 'full'); % State trajectory variables
V = sdpvar(nu, N, 'full');     % Input trajectory variables
X0 = sdpvar(nx, 1, 'full');    % Initial state variable
U_L = sdpvar(nu, 1, 'full');   % learned-input bound variable

Phi_x = sdpvar( (N + 1) * nx, (N + 1) * nx, 'full');
Phi_u = sdpvar( (N + 1) * nu, (N + 1) * nx, 'full');

% Construct the sigma matrix
sigma_seq = kron(eye(N), m.Bw);
Sigma_mat = blkdiag(eye(nx),sigma_seq);

% Define the objective function
objective =(V(:,1) - U_L)^2;

% Initialise the constraints
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

% state constraints
Fx = m.F_x;
bx = m.b_x;
nFx = length(bx);
for ii = 1:N
    for jj = 1: nFx
        f = Fx(jj,:); b = bx(jj);
        LHS = f*Z(:,ii);
        for kk = 1:ii-1
            LHS = LHS + norm(f*Phi_x((ii-1)*nx+1:ii*nx,kk*nx+1:(kk+1)*nx), 1);
        end
        constraints = [constraints, LHS <= b];
    end
end

% terminal constraint   
X_f = max_RPI.H;
Ft= X_f(:,1:nx);
bt = X_f(:,nx+1);
nFt = length(bt);
for jj = 1:nFt
    f = Ft(jj,:); b = bt(jj);
    LHS = f*Z(:,N+1);
    for kk = 1:N
        LHS = LHS + norm(f*Phi_x(N*nx+1:(N+1)*nx,kk*nx+1:(kk+1)*nx),1);
    end
    constraints = [constraints, LHS <= b];
end

% add input constraint
Fu = m.F_u;
bu = m.b_u;
nFu = length(bu);
for ii = 1:N
    for jj = 1: nFu
        f = Fu(jj,:); b = bu(jj);
        LHS = f*V(:,ii);
        for kk = 1:ii-1
            LHS = LHS + norm(f*Phi_u((ii-1)*nu+1:ii*nu,kk*nx+1:(kk+1)*nx),1);
        end
        constraints = [constraints, LHS <= b];
    end
end

options = sdpsettings('verbose',1,'solver','mosek');
sol_SL_MPSF = optimizer(constraints,objective,options,[X0;U_L],V(1));

yalmip('clear');
clear Z V Phi_x Phi_u X0 U_L objective constraints;
%% Initialisation - Explicit safety filter based on System Level Synthesis %%

Z = sdpvar(nx,N+1,'full');
V = sdpvar(nu,N,'full');
alpha = sdpvar(1,1,'full');

Phi_x = sdpvar( (N + 1) * nx, (N + 1) * nx, 'full');
Phi_u = sdpvar( (N + 1) * nu, (N + 1) * nx, 'full');
sigma_seq = kron(eye(N), m.Bw);
Sigma_mat = blkdiag(alpha*eye(nx),sigma_seq);

objective = -alpha;

constraints = [];
constraints = [constraints, alpha >= 0, alpha <=5];

for k = 1 : N
    constraints = [constraints, Phi_x( (k - 1)*nx + 1: k*nx, k*nx + 1: end) == zeros(nx, (N + 1 - k)*nx)];
end

for k = 1 : N
    constraints = [constraints, Phi_u( (k - 1)*nu + 1: k*nu, k*nx + 1: end) == zeros(nu, (N + 1 - k)*nx)];
end

% block downshift operator
Z_block = kron(diag(ones(1,N),-1), eye(nx));
ZA_block = Z_block*blkdiag(kron(eye(N), A), zeros(nx, nx));
ZB_block = Z_block*blkdiag(kron(eye(N), B), zeros(nx, nu));
Id = eye((N + 1)*nx);

% add affine constraint
constraints = [constraints, (Id - ZA_block)*Phi_x - ZB_block*Phi_u == Sigma_mat];
for k=1:N
    constraints = [ constraints, Z(:,k+1)==A*Z(:,k)+B*V(:,k)];
end

% state constraints
Fx = m.F_x;
bx = m.b_x;
nFx = length(bx);
for ii = 1:N
    for jj = 1: nFx
        f = Fx(jj,:); b = bx(jj);
        %LHS = f*Phi_x((ii-1)*nx+1:ii*nx,1:nx); % remove
        LHS = f*Z(:,ii) + norm(f*Phi_x((ii-1)*nx+1:ii*nx,1:nx),1);
        for kk = 1:ii-1
            LHS = LHS + norm(f*Phi_x((ii-1)*nx+1:ii*nx,kk*nx+1:(kk+1)*nx), 1);% indices?
        end
        constraints = [constraints, LHS <= b];
    end
end

% terminal constraint
Ft = eye(nx);
nFt = nx;
for jj = 1: nFt
    f = Ft(jj,:);
    LHS = norm(f*(Z(:,N+1) - Z(:,1)),1);
    for kk = 1:N
        LHS = LHS + norm(f*Phi_x(N*nx+1:(N+1)*nx,kk*nx+1:(kk+1)*nx),1);
    end
    constraints = [constraints, LHS <= alpha];
end

%add input constraint
Fu = m.F_u;
bu = m.b_u;
nFu = length(bu);
for ii = 1:N
    for jj = 1: nFu
        f = Fu(jj,:); b = bu(jj);
        LHS = f*V(:,ii) + norm(f*Phi_u((ii-1)*nu+1:ii*nu,1:nx),1);
        for kk = 1:ii-1
            LHS = LHS + norm(f*Phi_u((ii-1)*nu+1:ii*nu,kk*nx+1:(kk+1)*nx),1);
        end
        constraints = [constraints, LHS <= b];
    end
end

options = sdpsettings('verbose',1,'solver','mosek');
optimize(constraints,objective,options);

alpha_star = value(alpha);
clear Z V Phi_x Phi_u objective constraints alpha;

%% RESULTS %%
%% Performance explicit safe set
S_e_alpha_star = alpha_star*Polyhedron([eye(nx);-eye(nx)],ones(2*nx,1));

BS = minus(S_e_alpha_star, W);
UL= [-m.u_max, m.u_max];

grid_density = 100;
x1_range = linspace(-5,5,grid_density);
x2_range = linspace(-5,5,grid_density);

timming_explicit = [];

for i = 1:length(x1_range)
    for j = 1:length(x2_range)
        x0 = [x1_range(i); x2_range(j)];
        u = 2*m.u_max*rand(1)-m.u_max;
        tic
        BS.contains(A*x0 + B*u);
        U.contains(u);
        time =toc;
        timming_explicit = [timming_explicit;time];
    end
end

mean(timming_explicit)
std(timming_explicit)
%% Performance MPSF vs. SL-MPSF

max_v0_MPSF = nan(length(x1_range), length(x2_range));
max_v0_SL_MPSF = nan(length(x1_range), length(x2_range));

timming_MPSF = [];
timming_SL_MPSF  = [];

RoA_SL_MPSF = [];
RoA_MPSF = [];
for i = 1:length(x1_range)
    for j = 1:length(x2_range)
        x0 = [x1_range(i); x2_range(j)];
        tic
        [v_star, err] = sol([x0;m.u_max]);
        time = toc;
            timming_MPSF = [timming_MPSF;time];
        if err==0
            [v_star2, err] = sol([x0;-m.u_max]);
            RoA_MPSF = [RoA_MPSF,[x0;max(norm(v_star-m.u_max)^2,norm(v_star2+m.u_max)^2)]];
            max_v0_MPSF(i,j) = full(max(norm(v_star-m.u_max)^2,norm(v_star2+m.u_max)^2));
        end
        tic
        [v_star, err] = sol_SL_MPSF([x0;m.u_max]);
        time = toc;
            timming_SL_MPSF = [timming_SL_MPSF;time];
        if err==0
            [v_star2, err] = sol_SL_MPSF([x0;-m.u_max]);
            RoA_SL_MPSF = [RoA_SL_MPSF,[x0;max(norm(v_star-m.u_max)^2,norm(v_star2+m.u_max)^2)]];
            max_v0_SL_MPSF(i,j) = full(max(norm(v_star-m.u_max)^2,norm(v_star2+m.u_max)^2));
        end
    end
end

grid_density^2
mean(timming_SL_MPSF)
mean(timming_MPSF)

std(timming_SL_MPSF)
std(timming_MPSF)

disp(table([
mean(timming_explicit); 
mean(timming_SL_MPSF); mean(timming_MPSF)], [std(timming_explicit); std(timming_SL_MPSF); std(timming_MPSF)], 'VariableNames', {'mean', 'std'}, 'RowNames', {'Explicit','SL_MPSF', 'MPSF'}))


RoA_MPSF = RoA_MPSF'; 
RoA_SL_MPSF = RoA_SL_MPSF';
[X, Y] = meshgrid(x1_range, x2_range);
%ave('data.mat');
%% figure
clear all;
close all;
clc;
load('data.mat')
%%
f = figure(1);
clf;
% colors =     [0.0504    0.0298    0.5280
%     0.4176    0.0006    0.6584
%     0.6928    0.1651    0.5645
%     0.8814    0.3925    0.3832
%     0.9883    0.6523    0.2114
%     0.9400    0.9752    0.1313
%     ];
% colors =[    0.0504    0.0298    0.5280
%     0.3656    0.0030    0.6499
%     0.6107    0.0902    0.6200
%     0.7964    0.2780    0.4713
%     0.9283    0.4730    0.3261
%     0.9936    0.7018    0.1845
%     0.9400    0.9752    0.1313
%     ];

colors =[    0.2670    0.0049    0.3294
    0.2673    0.2259    0.5132
    0.1906    0.4071    0.5561
    0.1282    0.5651    0.5509
    0.2080    0.7187    0.4729
    0.5604    0.8413    0.2661
    0.9932    0.9062    0.1439
    ];

color_MPSF = colors(1,:);
color_Se = colors(2,:);
color_SL_MPSF = colors(3,:);
color_max_RPI = colors(4,:);
color_X = colors(5,:);
color_max_RCI = colors(6,:);


%[X, Y] = meshgrid(x1_range, x2_range);
subplot(1,3,1);
kav = convhull(RoA_MPSF(:,1:2));
hold on
[M,c] = contourf(X,Y, sqrt(max_v0_MPSF)');
num_colors = 256;
white_color = [1, 1, 1];
gray_color = [0.5, 0.5, 0.5];
custom_gray_colormap = interp1([0, 1], [gray_color; white_color], linspace(0, 1, num_colors));
colormap(custom_gray_colormap);
c.LineWidth = 0.1;
plot(RoA_MPSF(kav,1),RoA_MPSF(kav,2),'color',color_MPSF,'LineStyle','-.','Linewidth',3);
axis equal
xlim([-1.1*m.x_max 1.1*m.x_max]);
ylim([-1.1*m.x_max 1.1*m.x_max]);
rectangle('Position',[-m.x_max -m.x_max 2*m.x_max 2*m.x_max],'EdgeColor',	color_X,'LineStyle','--','Linewidth',2);
hline = line(NaN,NaN,'Color',	"#D95319",'LineStyle','--','Linewidth',2);
set(gca, 'XTick', [], 'YTick', []);
axis off;

s = subplot(1,3,2);
kav = convhull(RoA_SL_MPSF(:,1:2));
hold on
[M,c] =contourf(X, Y, sqrt(max_v0_SL_MPSF)');
colormap(custom_gray_colormap);
c.LineWidth = 0.2;
plot(RoA_SL_MPSF(kav,1),RoA_SL_MPSF(kav,2),'color',color_SL_MPSF,'LineStyle','-','Linewidth',3);
axis equal
xlim([-1.1*m.x_max 1.1*m.x_max]);
ylim([-1.1*m.x_max 1.1*m.x_max]);
rectangle('Position',[-m.x_max -m.x_max 2*m.x_max 2*m.x_max],'EdgeColor',	color_X,'LineStyle','--','Linewidth',2);
hline = line(NaN,NaN,'Color',"#D95319",'LineStyle','--','Linewidth',2);
set(gca,'FontSize',8)
set(gca, 'XTick', [], 'YTick', []);
axis off;

h = axes(f,'visible','off'); 
c = colorbar(h,'south');
colormap(c,custom_gray_colormap)
c.Label.Interpreter = 'latex';
c.Label.String = 'max $|u_\mathcal{L} - \mathbf{v}_0|$';
c.Ticks = [0 0.5 1];
c.TickLabels = {"0", "3", '6'};
c.Label.FontSize = 12;
c.Position = c.Position*diag([1,1,1,0.5]);

s = subplot(1,3,3);
hold on;
kav = convhull(max_RCI.V);
x = max_RCI.V(kav, 1);
y = max_RCI.V(kav, 2);
transparency_color_max_RCI = 0.65;
fill(x, y, color_max_RCI, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_max_RCI);

alpha = alpha_star;
transparency_color_Se = 1;
x = [-alpha, alpha, alpha, -alpha, -alpha];
y = [-alpha, -alpha, alpha, alpha, -alpha];
fill(x, y, color_Se, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_Se);

% MAX RPI
kav = convhull(max_RPI.V);
x = max_RPI.V(kav, 1);
y = max_RPI.V(kav, 2);
transparency_color_max_RPI = 0.8;
fill(x, y, color_max_RPI, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_max_RPI);


kav = convhull(RoA_SL_MPSF(:,1:2));
plot(RoA_SL_MPSF(kav,1),RoA_SL_MPSF(kav,2),'color',color_SL_MPSF,'LineStyle','-','Linewidth',3);

kav = convhull(RoA_MPSF(:,1:2));
hold on
plot(RoA_MPSF(kav,1),RoA_MPSF(kav,2),'color',color_MPSF,'LineStyle','-.','Linewidth',3);
rectangle('Position',[-m.x_max -m.x_max 2*m.x_max 2*m.x_max],'EdgeColor',	color_X,'LineStyle','--','Linewidth',2);
hline = line(NaN,NaN,'Color',	"#D95319",'LineStyle','-.','Linewidth',2);

axis equal;
xlim([-1.1*m.x_max 1.1*m.x_max]);
ylim([-1.1*m.x_max 1.1*m.x_max]);

set(gca, 'XTick', [], 'YTick', []);
axis off;

dummy_axes = axes('Position', [0.92 0.3 0.1 0.1], 'visible', 'off');
hold(dummy_axes, 'on');

dummy_line1 = line(dummy_axes, NaN, NaN, 'Color', colors(1,:), 'LineStyle', '-.', 'LineWidth', 2);
dummy_line2 = fill(nan, nan, color_Se, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_Se);
dummy_line3 = line(dummy_axes, NaN, NaN, 'Color', colors(3,:), 'LineStyle', '-', 'LineWidth', 2);
dummy_line4 = fill(nan, nan, color_max_RPI, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_max_RPI);
dummy_line5 = line(dummy_axes, NaN, NaN, 'Color', color_X, 'LineStyle', '--', 'LineWidth', 2);
dummy_line6 = fill(nan, nan, color_max_RCI, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_max_RCI);

l = legend(dummy_axes, [dummy_line1, dummy_line2, dummy_line3, dummy_line4, dummy_line5, dummy_line6], ...
    {'$\mathcal{S}_\mathrm{MPSF}$','$\mathcal{S}_\mathrm{e}$', '$\mathcal{S}_\mathrm{SL-MPSF}$',  '$\Omega_\mathrm{max}$','$\mathcal{X}$', '$\Xi_\mathrm{max}$'}, ...
    'Interpreter', 'latex', 'FontSize', 12, 'Location', 'east', 'Box', 'off','Position',[0.8997 0.3366 0.1216 0.3819]);

dummy_axes = axes('Position', [0.05 0.325 0.1 0.35], 'visible', 'off');
hold(dummy_axes, 'on');

dummy_colorbar = colorbar(dummy_axes, 'Location', 'west');
colormap(dummy_colorbar, custom_gray_colormap);
dummy_colorbar.Label.Interpreter = 'latex';
dummy_colorbar.Label.String = 'max $|u_\mathcal{L} - \mathbf{v}_0|$';
dummy_colorbar.Ticks = [0 0.5 1];
dummy_colorbar.TickLabels = {"0", "3", '6'};
dummy_colorbar.Label.FontSize = 12;

set(gca,'FontSize',8)
set(gcf,'units','centimeters','Position', [0 0 25 15]);

exportgraphics(gcf,strcat('fig1.pdf'),'ContentType','vector');

