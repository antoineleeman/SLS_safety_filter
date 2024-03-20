clear all;
close all;
clc;
load('data_paper.mat')

f = figure(1);
clf;

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
color_max_PI = colors(4,:);
color_X = colors(5,:);
color_max_RCI = colors(6,:);


B_inf = Polyhedron('lb', -1*ones(nx*(N+1),1), 'ub', ones(nx*(N+1),1));
c = magma(N);
transparency_color_Se = 0.7;
Linewidth_x = 3;
%%
figure(1);
clf;
c = magma(N);
hold on;
rectangle('Position',[-m.x_max -m.x_max 2*m.x_max 2*m.x_max],'EdgeColor',	color_X,'LineStyle','--','Linewidth',Linewidth_x);
kav = convhull(RoA_SL_MPSF(:,1:2));
plot(RoA_SL_MPSF(kav,1),RoA_SL_MPSF(kav,2),'color',color_SL_MPSF,'LineStyle','-','Linewidth',3);

kav = convhull(max_RPI.V);
x = max_RPI.V(kav, 1);
y = max_RPI.V(kav, 2);
transparency_color_max_RPI = 0.3;
fill(x, y, color_max_RPI, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_max_RPI);


for kk=1:N    
    DRSkk = (m.Bw*Phi_x_example( (kk - 1)*nx + 1: kk*nx,:)*B_inf);
    DRSkk.computeVRep;
    Vk = DRSkk.V;
    if length(Vk) <=2
        continue
    end
    kav = convhull(Vk);    
    patch(Vk(kav,1) +Z_example(1,kk+1) ,Vk(kav,2)+Z_example(2,kk+1),c(kk,:),'FaceAlpha',0.5);
end

p = plot(Z_example(1,:), Z_example(2,:),'--ko','Linewidth',2 );

set(gca,'FontSize',8)
set(gcf,'units','centimeters','Position', [0 0 17 17]);

set(gca, 'XTick', [], 'YTick', []);
axis off;
axis equal;

h=gcf;
%h.Color = 'none';  % Sets the background color to none
exportgraphics(gcf,strcat('poster_fig1.pdf'),'ContentType','vector');
%saveas(gcf, 'poster_fig1.png');
exportgraphics(h, 'poster_fig1.png', 'BackgroundColor', 'none','Resolution',600)

%%
figure(2);
clf;
hold on;
rectangle('Position',[-m.x_max -m.x_max 2*m.x_max 2*m.x_max],'EdgeColor',color_X,'LineStyle','--','Linewidth',Linewidth_x);

x = [-alpha_star, alpha_star, alpha_star, -alpha_star, -alpha_star];
y = [-alpha_star, -alpha_star, alpha_star, alpha_star, -alpha_star];
fill(x, y, color_Se, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_Se);

for kk=2:N    
    DRSkk = (Phi_x_explicit( (kk - 1)*nx + 1: kk*nx,:)*B_inf);
    DRSkk.computeVRep;
    Vk = DRSkk.V;
    if length(Vk) <=2
        continue
    end
    kav = convhull(Vk);    
    patch(Vk(kav,1) ,Vk(kav,2),c(kk,:),'FaceAlpha',0.5); % assume Z=0
end

Z_cl_vec = Phi_x_explicit*reshape([ones(2,1), ones(2,N)],[nx*(N+1),1]);
Z_cl = reshape(Z_cl_vec,[nx, N+1]);

plot(Z_cl(1,:),Z_cl(2,:),'-k','Linewidth',2);

set(gcf,'units','centimeters','Position', [0 0 17 17]);

set(gca, 'XTick', [], 'YTick', []);
axis off;
axis equal;
h=gcf;
exportgraphics(gcf,strcat('poster_fig2.pdf'),'ContentType','vector');
%saveas(gcf, 'poster_fig1.png');
exportgraphics(h, 'poster_fig2.png', 'BackgroundColor', 'none','Resolution',600);

%%

figure(3);
clf;
% hold on;
% kav = convhull(max_RCI.V);
% x = max_RCI.V(kav, 1);
% y = max_RCI.V(kav, 2);
% transparency_color_max_RCI = 0.65;
% fill(x, y, color_max_RCI, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_max_RCI);

alpha = alpha_star;
x = [-alpha, alpha, alpha, -alpha, -alpha];
y = [-alpha, -alpha, alpha, alpha, -alpha];
fill(x, y, color_Se, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_Se);
hold on
% MAX RPI
% kav = convhull(max_RPI.V);
% x = max_RPI.V(kav, 1);
% y = max_RPI.V(kav, 2);
% transparency_color_max_RPI = 0.8;
% fill(x, y, color_max_RPI, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_max_RPI);

kav = convhull(RoA_SL_MPSF(:,1:2));
plot(RoA_SL_MPSF(kav,1),RoA_SL_MPSF(kav,2),'color',color_SL_MPSF,'LineStyle','-','Linewidth',3);

kav = convhull(RoA_MPSF(:,1:2));
hold on
plot(RoA_MPSF(kav,1),RoA_MPSF(kav,2),'color',color_MPSF,'LineStyle','-.','Linewidth',3);
rectangle('Position',[-m.x_max -m.x_max 2*m.x_max 2*m.x_max],'EdgeColor',	color_X,'LineStyle','--','Linewidth',Linewidth_x);
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
%dummy_line4 = fill(nan, nan, color_max_RPI, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_max_RPI);
dummy_line5 = line(dummy_axes, NaN, NaN, 'Color', color_X, 'LineStyle', '--', 'LineWidth', 2);
%dummy_line6 = fill(nan, nan, color_max_RCI, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_max_RCI);

l = legend(dummy_axes, [dummy_line1, dummy_line2, dummy_line3, dummy_line5], ...
    {'[1]','Offline', 'Online','$\mathcal{X}$'}, ...
    'Interpreter', 'latex', 'FontSize', 12, 'Location', 'east', 'Box', 'off','Position',[0.7472 0.3337 0.1286 0.3819]);

dummy_axes = axes('Position', [0.05 0.325 0.1 0.35], 'visible', 'off');
hold(dummy_axes, 'on');

set(gca,'FontSize',8)
set(gcf,'units','centimeters','Position', [0 0 20 12]);
h=gcf;
exportgraphics(gcf,strcat('poster_fig3.pdf'),'ContentType','vector');
exportgraphics(h, 'poster_fig3.png', 'BackgroundColor', 'none' ,'Resolution',600);

%%

f = figure(4);
clf;
subplot(1,2,1);
title('State-of-the-art [1]','FontName','arial')
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
rectangle('Position',[-m.x_max -m.x_max 2*m.x_max 2*m.x_max],'EdgeColor',	color_X,'LineStyle','--','Linewidth',Linewidth_x);
hline = line(NaN,NaN,'Color',	"#D95319",'LineStyle','--','Linewidth',2);
set(gca, 'XTick', [], 'YTick', []);
set(gca,'FontSize',12)
axis off;

s = subplot(1,2,2);
title('Online method','FontName','arial','fontsize',18)
kav = convhull(RoA_SL_MPSF(:,1:2));
hold on
[M,c] =contourf(X, Y, sqrt(max_v0_SL_MPSF)');
colormap(custom_gray_colormap);
c.LineWidth = 0.2;
plot(RoA_SL_MPSF(kav,1),RoA_SL_MPSF(kav,2),'color',color_SL_MPSF,'LineStyle','-','Linewidth',3);
axis equal
xlim([-1.1*m.x_max 1.1*m.x_max]);
ylim([-1.1*m.x_max 1.1*m.x_max]);
rectangle('Position',[-m.x_max -m.x_max 2*m.x_max 2*m.x_max],'EdgeColor',	color_X,'LineStyle','--','Linewidth',Linewidth_x);
hline = line(NaN,NaN,'Color',"#D95319",'LineStyle','--','Linewidth',2);
set(gca,'FontSize',12)
set(gca, 'XTick', [], 'YTick', []);
axis off;

h = axes(f,'visible','off'); 
c = colorbar(h,'south');
colormap(c,custom_gray_colormap)
c.Label.Interpreter = 'latex';
c.Label.String = 'max. filter intervention';
c.Ticks = [0 0.5 1];
c.TickLabels = {"0", "3", '6'};
c.FontSize = 10;
c.Label.FontSize = 16;
c.Position = c.Position*diag([1,1,1,0.5]);
c.FontName = 'Arial';


% dummy_axes = axes('Position', [0.92 0.3 0.1 0.1], 'visible', 'off');
% hold(dummy_axes, 'on');

% dummy_line1 = line(dummy_axes, NaN, NaN, 'Color', colors(1,:), 'LineStyle', '-.', 'LineWidth', 2);
% dummy_line2 = fill(nan, nan, color_Se, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_Se);
% dummy_line3 = line(dummy_axes, NaN, NaN, 'Color', colors(3,:), 'LineStyle', '-', 'LineWidth', 2);
% dummy_line4 = fill(nan, nan, color_max_RPI, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_max_RPI);
% dummy_line5 = line(dummy_axes, NaN, NaN, 'Color', color_X, 'LineStyle', '--', 'LineWidth', 2);
% dummy_line6 = fill(nan, nan, color_max_RCI, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_max_RCI);
% 
% l = legend(dummy_axes, [dummy_line1, dummy_line2, dummy_line3, dummy_line4, dummy_line5, dummy_line6], ...
%     {'$\mathcal{S}_\mathrm{MPSF}$','$\mathcal{S}_\mathrm{e}$', '$\mathcal{S}_\mathrm{SL-MPSF}$',  '$\Omega_\mathrm{max}$','$\mathcal{X}$', '$\Xi_\mathrm{max}$'}, ...
%     'Interpreter', 'latex', 'FontSize', 12, 'Location', 'east', 'Box', 'off','Position',[0.8997 0.3366 0.1216 0.3819]);
% 
% dummy_axes = axes('Position', [0.05 0.325 0.1 0.35], 'visible', 'off');
% hold(dummy_axes, 'on');

% dummy_colorbar = colorbar(dummy_axes, 'Location', 'west');
% colormap(dummy_colorbar, custom_gray_colormap);
% dummy_colorbar.Label.Interpreter = 'latex';
% dummy_colorbar.Label.String = 'max intervention';
% dummy_colorbar.Ticks = [0 0.5 1];
% dummy_colorbar.TickLabels = {"0", "3", '6'};
% dummy_colorbar.Label.FontSize = 12;

set(gca,'FontSize',8)
set(gcf,'units','centimeters','Position', [0 0 18 13]);
h=gcf;
exportgraphics(gcf,strcat('poster_fig4.pdf'),'ContentType','vector');
exportgraphics(h, 'poster_fig4.png', 'BackgroundColor', 'none','Resolution',600);

% %%
% figure(5);
% clf;
% c = magma(N);
% hold on;
% rectangle('Position',[-m.x_max -m.x_max 2*m.x_max 2*m.x_max],'EdgeColor',	color_X,'LineStyle','--','Linewidth',Linewidth_x);
% kav = convhull(RoA_MPSF(:,1:2));
% plot(RoA_MPSF(kav,1),RoA_MPSF(kav,2),'color',color_MPSF,'LineStyle','-','Linewidth',3);
% 
% kav = convhull(max_PI.V);
% x = max_PI.V(kav, 1);
% y = max_PI.V(kav, 2);
% transparency_color_max_PI = 0.3;
% fill(x, y, color_max_PI, 'EdgeColor', 'none', 'FaceAlpha', transparency_color_max_PI);
% 
% kav_E = convhull(E.V);
% x_E = E.V(kav_E, 1);
% y_E = E.V(kav_E, 2);
% 
% 
% for kk=1:N    
% %     DRSkk = (m.Bw*Phi_x_example( (kk - 1)*nx + 1: kk*nx,:)*B_inf);
% %     DRSkk.computeVRep;
% %     Vk = DRSkk.V;
% %     if length(Vk) <=2
% %         continue
% %     end
% %     kav = convhull(Vk);    
%     patch(x_E +Z_example_MPC(1,kk+1) ,y_E+Z_example_MPC(2,kk+1),c(kk,:),'FaceAlpha',0.5);
% end
% 
% p = plot(Z_example_MPC(1,:), Z_example_MPC(2,:),'--ko','Linewidth',2 );
% 
% set(gca,'FontSize',8)
% set(gcf,'units','centimeters','Position', [0 0 17 17]);
% 
% set(gca, 'XTick', [], 'YTick', []);
% axis off;
% axis equal;
% 
% h=gcf;
% %h.Color = 'none';  % Sets the background color to none
% exportgraphics(gcf,strcat('poster_fig5.pdf'),'ContentType','vector');
% %saveas(gcf, 'poster_fig1.png');
% exportgraphics(h, 'poster_fig5.png', 'BackgroundColor', 'none','Resolution',600)

%% table:


disp(table([mean(timming_explicit); mean(timming_SL_MPSF); mean(timming_MPSF)],...
    [std(timming_explicit); std(timming_SL_MPSF); std(timming_MPSF)],...
    [S_e_alpha_star.volume/ max_RCI.volume; Polyhedron(RoA_SL_MPSF(:,1:2)).volume/max_RCI.volume; Polyhedron(RoA_MPSF(:,1:2)).volume/max_RCI.volume], 'VariableNames', {'mean', 'std', 'volume ratio'}, 'RowNames', {'Explicit','SL_MPSF', 'MPSF'}))




