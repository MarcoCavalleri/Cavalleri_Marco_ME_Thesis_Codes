%% VISUALIZATION E-TYPE

% This is the same code for every steady state, it will
% change only the datadir, filenames and the value of chi and simulation
% parameters 

clear; clc; close all;

% For chemical potential visualization (to avoid numerical errors)
const_chem = false;

% For visualize norm of chemical potentials
norm_chem = false;

%% Load data
D_2 = 1;        % Diffusion constant for c_2
mu0_2 = -2;
chi = 1;
k1 = 0.5;
k2 = 0.1;
s_2 = 0.001;
s_3 = 0.001;
s_4 = 0.001;
mu_Y2 = -1;

% --- File paths ---
datadir = 'C:\Users\user\Desktop\UniPd_Thesis\Pb0\fig1\results_fig1_ss1\'; 

filenames = {'c_1_fig1','c_2_fig1','c_3_fig1','c_s_fig1'};

% --- Load DATA ---
for i = 1:4
    opts=detectImportOptions(fullfile(datadir, [filenames{i}]), 'FileType', 'text', 'CommentStyle','#');
    T=readtable(fullfile(datadir, [filenames{i}]),opts);
    X{i} = T{:,1};
    Y{i} = T{:,2};
    C{i} = T{:,3};
    M{i} = T{:,4};    
    D{i} = T{:,5};
end

x = X{2};
y = Y{2};

N = sqrt(length(x));
assert(N==round(N),'Not a square grid!');

xvals = unique(x);  % sorted
yvals = unique(y);  % sorted
[Xg, Yg] = meshgrid(linspace(min(xvals),max(xvals),N), linspace(min(yvals),max(yvals),N));

color1 = [250 221  95]/255; % #CADB37
color2 = [228 156 149]/255; % #E49C95
color3 = [109 114 209]/255; % #92ABD3

customMap = interp1( ...
    [0 0.5 1], ...
    [color3; color2; color1], ...
    linspace(0,1,256) ...
);

...

dx_i = 1 - 15/30;
dx_f = 1 - dx_i;
dy_i = 1 - 30/30;
dy_f = 1 - dy_i;

% Reshape fields into 2D matrices
c_1g_m = reshape(C{1},N,N);
c_2g_m = reshape(C{2},N,N);
c_3g_m = reshape(C{3},N,N);
c_sg_m = reshape(C{3},N,N);

c_1g = zeros(N,N);
c_1g(1:dx_i*N,1:dy_i*N) = c_1g_m(dx_f*N+1:N,dy_f*N+1:N);
c_1g(1:dx_i*N,dy_i*N+1:N) = c_1g_m(dx_f*N+1:N,1:dy_f*N);
c_1g(dx_i*N+1:N,1:dy_i*N) = c_1g_m(1:dx_f*N,dy_f*N+1:N);
c_1g(dx_i*N+1:N,dy_i*N+1:N) = c_1g_m(1:dx_f*N,1:dy_f*N);

c_2g = zeros(N,N);
c_2g(1:dx_i*N,1:dy_i*N) = c_2g_m(dx_f*N+1:N,dy_f*N+1:N);
c_2g(1:dx_i*N,dy_i*N+1:N) = c_2g_m(dx_f*N+1:N,1:dy_f*N);
c_2g(dx_i*N+1:N,1:dy_i*N) = c_2g_m(1:dx_f*N,dy_f*N+1:N);
c_2g(dx_i*N+1:N,dy_i*N+1:N) = c_2g_m(1:dx_f*N,1:dy_f*N);

c_3g = zeros(N,N);
c_3g(1:dx_i*N,1:dy_i*N) = c_3g_m(dx_f*N+1:N,dy_f*N+1:N);
c_3g(1:dx_i*N,dy_i*N+1:N) = c_3g_m(dx_f*N+1:N,1:dy_f*N);
c_3g(dx_i*N+1:N,1:dy_i*N) = c_3g_m(1:dx_f*N,dy_f*N+1:N);
c_3g(dx_i*N+1:N,dy_i*N+1:N) = c_3g_m(1:dx_f*N,1:dy_f*N);

c_sg = zeros(N,N);
c_sg(1:dx_i*N,1:dy_i*N) = c_sg_m(dx_f*N+1:N,dy_f*N+1:N);
c_sg(1:dx_i*N,dy_i*N+1:N) = c_sg_m(dx_f*N+1:N,1:dy_f*N);
c_sg(dx_i*N+1:N,1:dy_i*N) = c_sg_m(1:dx_f*N,dy_f*N+1:N);
c_sg(dx_i*N+1:N,dy_i*N+1:N) = c_sg_m(1:dx_f*N,1:dy_f*N);

[dc_1g, dc_2g, dc_3g, dc_sg, m_1g, m_2g, m_3g, m_sg] = compute_derivatives(c_1g, c_2g, c_3g, c_sg, norm_chem);

sigma = 100;
gaussianKernel = exp(-(Xg.^2 + Yg.^2)/(2*sigma^2));
gaussianKernel = gaussianKernel / sum(gaussianKernel(:));  

dc_1g = conv2(dc_1g, gaussianKernel, 'same');
dc_2g = conv2(dc_2g, gaussianKernel, 'same');
dc_3g = conv2(dc_3g, gaussianKernel, 'same');
dc_sg = conv2(dc_sg, gaussianKernel, 'same');

if const_chem
    m_1g = conv2(m_1g, gaussianKernel, 'same');
    m_2g = conv2(m_2g, gaussianKernel, 'same');
    m_3g = conv2(m_3g, gaussianKernel, 'same');
    m_sg = conv2(m_sg, gaussianKernel, 'same');
end

% Concentrations
figure('Visible', 'on');
subplot(2,2,1); 
surf(Xg, Yg, c_1g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$c_1$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'a)', 'Position', [min(xvals)-9, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(2,2,2);
surf(Xg, Yg, c_2g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16);
title('$c_2$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'b)', 'Position', [min(xvals)-9, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(2,2,3); 
surf(Xg, Yg, c_3g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$c_3$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'c)', 'Position', [min(xvals)-9, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(2,2,4); 
surf(Xg, Yg, c_sg,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$c_s$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'd)', 'Position', [min(xvals)-9, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

% Derivatives
figure('Visible', 'off');
subplot(2,2,1); 
surf(Xg, Yg, dc_1g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$\frac{dc_1}{dt}$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'a)', 'Position', [min(xvals)-9, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(2,2,2);
surf(Xg, Yg, dc_2g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16);
title('$\frac{dc_2}{dt}$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'b)', 'Position', [min(xvals)-9, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(2,2,3); 
surf(Xg, Yg, dc_3g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$\frac{dc_3}{dt}$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'c)', 'Position', [min(xvals)-9, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(2,2,4); 
surf(Xg, Yg, dc_sg,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$\frac{dc_s}{dt}$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'd)', 'Position', [min(xvals)-9, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

% Chemical potentials
figure('Visible', 'off');
subplot(2,2,1); 
surf(Xg, Yg, m_1g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$\mu_1$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'a)', 'Position', [min(xvals)-9, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(2,2,2);
surf(Xg, Yg, m_2g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16);
title('$\mu_2$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'b)', 'Position', [min(xvals)-9, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(2,2,3); 
surf(Xg, Yg, m_3g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$\mu_3$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'c)', 'Position', [min(xvals)-9, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(2,2,4); 
surf(Xg, Yg, m_sg,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$\mu_s$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'd)', 'Position', [min(xvals)-9, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

%% Computation chemical potentials and derivatives
function [dc_1g, dc_2g, dc_3g, dc_sg, m_1g, m_2g, m_3g, m_sg] = compute_derivatives(c_1g, c_2g, c_3g, c_sg, norm_chem)
    
    % Parameters
    D_1 = 1; D_2 = 1; D_3 = 1; D_s = 1;
    mu0_1 = 0; mu0_2 = -2; mu0_3 = 0; mu0_s = 0;
    k1 = 0.5; k2 = 0.1;
    s = 1e-3; s_1 = s; s_2 = s; s_3 = s; s_4 = s; s_5 = s;
    mu_Y1 = 1; mu_Y2 = -1;
    l = 1; % chi
    dx = 0.33; % grid size
    
    % Finite differences under PBCs
    laplace = @(f) (circshift(f, [1 0]) + circshift(f, [-1 0]) + ...
                   circshift(f, [0 1]) + circshift(f, [0 -1]) - 4*f)/dx^2;
    
    gradient_x = @(f) (0.5 * (circshift(f, [-1 0]) - circshift(f, [1 0])))/dx;
    gradient_y = @(f) (0.5 * (circshift(f, [0 -1]) - circshift(f, [0 1])))/dx;
    
    dot_grad = @(fx, fy, gx, gy) fx .* gx + fy .* gy;
    
    % Chemical potentials
    lap_c1 = laplace(c_1g);
    lap_c2 = laplace(c_2g);
    lap_c3 = laplace(c_3g);
    lap_cs = laplace(c_sg);
    
    m_1g = mu0_1 + log(c_1g) + l * c_2g - k1 * lap_c1 - k2 * lap_c2;
    m_2g = mu0_2 + log(c_2g) + l * (c_1g + c_3g + c_sg) - k1 * lap_c2 - k2 * laplace(c_1g + c_3g + c_sg);
    m_3g = mu0_3 + log(c_3g) + l * c_2g - k1 * lap_c3 - k2 * lap_c2;
    m_sg = mu0_s + log(c_sg) + l * c_2g - k1 * lap_cs - k2 * lap_c2;
    
    % Reaction currents
    j_1 = s_1 * (exp(mu_Y1) - exp(m_1g));
    j_2 = s_2 * (exp(m_1g) - exp(m_2g));
    j_3 = s_3 * (exp(m_2g) - exp(mu_Y2));
    j_4 = s_4 * (exp(m_2g) - exp(m_3g));
    j_5 = s_5 * (exp(m_3g) - exp(m_1g));
    
    chem_term_1 = j_1 - j_2 + j_5;
    chem_term_2 = j_2 - j_3 - j_4;
    chem_term_3 = j_4 - j_5;
    
    % Current fluxes
    [m1_x, m1_y] = deal(gradient_x(m_1g), gradient_y(m_1g));
    [m2_x, m2_y] = deal(gradient_x(m_2g), gradient_y(m_2g));
    [m3_x, m3_y] = deal(gradient_x(m_3g), gradient_y(m_3g));
    [ms_x, ms_y] = deal(gradient_x(m_sg), gradient_y(m_sg));
    
    [c1_x, c1_y] = deal(gradient_x(c_1g), gradient_y(c_1g));
    [c2_x, c2_y] = deal(gradient_x(c_2g), gradient_y(c_2g));
    [c3_x, c3_y] = deal(gradient_x(c_3g), gradient_y(c_3g));
    [cs_x, cs_y] = deal(gradient_x(c_sg), gradient_y(c_sg));
    
    diff_term_1 = D_1 * (dot_grad(c1_x, c1_y, m1_x, m1_y) + c_1g .* laplace(m_1g));
    diff_term_2 = D_2 * (dot_grad(c2_x, c2_y, m2_x, m2_y) + c_2g .* laplace(m_2g));
    diff_term_3 = D_3 * (dot_grad(c3_x, c3_y, m3_x, m3_y) + c_2g .* laplace(m_3g));
    diff_term_s = D_s * (dot_grad(cs_x, cs_y, ms_x, ms_y) + c_sg .* laplace(m_sg));
    
    % Temporal derivative
    dc_1g = diff_term_1 + chem_term_1;
    dc_2g = diff_term_2 + chem_term_2;
    dc_3g = diff_term_3 + chem_term_3;
    dc_sg = diff_term_s;
    
    if norm_chem
        m_1g = sqrt(gradient_x(m_1g).^2 + gradient_y(m_1g).^2);
        m_2g = sqrt(gradient_x(m_2g).^2 + gradient_y(m_2g).^2);
        m_3g = sqrt(gradient_x(m_3g).^2 + gradient_y(m_3g).^2);
        m_sg = sqrt(gradient_x(m_sg).^2 + gradient_y(m_sg).^2);
    end
end
