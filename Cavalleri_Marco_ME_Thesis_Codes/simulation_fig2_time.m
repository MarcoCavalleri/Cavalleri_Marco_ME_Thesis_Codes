%% VISUALIZATION BRUSSELLATOR 

% This is the same code for every steady state, it will
% change only the datadir, filenames and the value of chi and simulation
% parameters 

clear; clc; close all;

% For chemical potential visualization (to avoid numerical errors)
const_chem = false;

% For visualize norm of chemical potentials
norm_chem = false;

%% Load data
datadir = 'C:\Users\user\Desktop\UniPd_Thesis\Pb0\fig2\results_fig2'; 

filenames = {'c_1_fig2','c_2_fig2','c_s_fig2'};

for i = 1:3
    opts = detectImportOptions(fullfile(datadir, [filenames{i} '.txt']), 'FileType', 'text', 'CommentStyle','#');
    T = readtable(fullfile(datadir, [filenames{i} '.txt']), opts);
    X{i} = T{:,1};
    Y{i} = T{:,2};
    C{i} = T{:,3};
    M{i} = T{:,4};
    D{i} = T{:,5};
end

% Same x and y, generate grid
x = X{2};
y = Y{2};

N = sqrt(length(x));
assert(N==round(N),'Not a square grid!');

xvals = unique(x);
yvals = unique(y);
[Xg, Yg] = meshgrid(linspace(min(xvals),max(xvals),N), linspace(min(yvals),max(yvals),N));

color1 = [250 221  95]/255;
color2 = [228 156 149]/255;
color3 = [109 114 209]/255;

customMap = interp1( ...
    [0 0.5 1], ...
    [color3; color2; color1], ...
    linspace(0,1,256) ...
);

...

% For rearrangement
dx_i = 1 - 30/30;
dx_f = 1 - dx_i;
dy_i = 1 - 30/30;
dy_f = 1 - dy_i;

% Concentrations
c_1g_m = reshape(C{1},N,N);
c_2g_m = reshape(C{2},N,N);
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

c_sg = zeros(N,N);
c_sg(1:dx_i*N,1:dy_i*N) = c_sg_m(dx_f*N+1:N,dy_f*N+1:N);
c_sg(1:dx_i*N,dy_i*N+1:N) = c_sg_m(dx_f*N+1:N,1:dy_f*N);
c_sg(dx_i*N+1:N,1:dy_i*N) = c_sg_m(1:dx_f*N,dy_f*N+1:N);
c_sg(dx_i*N+1:N,dy_i*N+1:N) = c_sg_m(1:dx_f*N,1:dy_f*N);

% Compute derivative and chemical potentials
[dc_1g, dc_2g, dc_sg, m_1g, m_2g, m_sg] = compute_derivatives(c_1g, c_2g, c_sg, norm_chem);

sigma = 10;
gaussianKernel = exp(-(Xg.^2 + Yg.^2)/(2*sigma^2));
gaussianKernel = gaussianKernel / sum(gaussianKernel(:));  

dc_1g = conv2(dc_1g, gaussianKernel, 'same');
dc_2g = conv2(dc_2g, gaussianKernel, 'same');
dc_sg = conv2(dc_sg, gaussianKernel, 'same');

if const_chem
    m_1g = conv2(m_1g, gaussianKernel, 'same');
    m_2g = conv2(m_2g, gaussianKernel, 'same');
    m_sg = conv2(m_sg, gaussianKernel, 'same');
end

%% Plots

% Concentration
figure('Visible', 'on');

subplot(1,3,1); 
surf(Xg, Yg, c_1g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$c_1$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'a)', 'Position', [min(xvals)-8, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(1,3,2);
surf(Xg, Yg, c_2g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16);
title('$c_2$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'b)', 'Position', [min(xvals)-8, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(1,3,3); 
surf(Xg, Yg, c_sg,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$c_s$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'c)', 'Position', [min(xvals)-8, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

 
% Derivatives
figure('Visible', 'on');

subplot(1,3,1); 
surf(Xg, Yg, dc_1g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$\frac{dc_1}{dt}$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'a)', 'Position', [min(xvals)-8, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(1,3,2);
surf(Xg, Yg, dc_2g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16);
title('$\frac{dc_2}{dt}$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'b)', 'Position', [min(xvals)-8, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(1,3,3); 
surf(Xg, Yg, dc_sg,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex', 'FontSize', 16); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$\frac{dc_s}{dt}$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
text('String', 'c)', 'Position', [min(xvals)-8, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

% Chemical potentials
figure('Visible', 'on');
subplot(1,3,2);
surf(Xg, Yg, m_2g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex'); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16);
title('$\mu_2$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
text('String', 'b)', 'Position', [min(xvals)-8, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(1,3,1); 
surf(Xg, Yg, m_1g,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex'); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$\mu_1$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
text('String', 'a)', 'Position', [min(xvals)-8, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

subplot(1,3,3); 
surf(Xg, Yg, m_sg,'EdgeColor','none'); 
colormap(customMap); colorbar('TickLabelInterpreter', 'latex'); view(2)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16); ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 16); 
title('$\mu_s$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
text('String', 'c)', 'Position', [min(xvals)-8, max(yvals)+2, 1], 'FontSize', 18, 'Interpreter', 'latex');
axis([0 30 0 30]);

%% Computation chemical potentials and derivatives
function [dc_1g, dc_2g, dc_sg, m_1g, m_2g, m_sg] = compute_derivatives(c_1g, c_2g, c_sg, norm_chem)
    
    % Parameters
    D_1 = 0.5; D_2 = 10; D_s = 1;
    mu0_1 = 0; mu0_2 = 0; mu0_s = 0;
    k1 = 0.1; k2 = 0.05;
    s = 1; s_1 = s; s_2 = s; s_3 = s; s_4 = s;
    mu_Y1 = log(2); mu_Y2 = log(5); mu_Y3 = log(0.01); mu_Y4 = log(0.1);
    l = 0; % chi
    dx = 0.25; % mesh size
    
    % Finite differences under PBCs
    laplace = @(f) (circshift(f, [1 0]) + circshift(f, [-1 0]) + ...
                   circshift(f, [0 1]) + circshift(f, [0 -1]) - 4*f)/dx^2;
    
    gradient_x = @(f) (0.5 * (circshift(f, [-1 0]) - circshift(f, [1 0])))/dx;
    gradient_y = @(f) (0.5 * (circshift(f, [0 -1]) - circshift(f, [0 1])))/dx;
    
    dot_grad = @(fx, fy, gx, gy) fx .* gx + fy .* gy;
    
    % Chemical potentials
    lap_c1 = laplace(c_1g);
    lap_c2 = laplace(c_2g);
    lap_cs = laplace(c_sg);
    
    m_1g = mu0_1 + log(c_1g) + l * c_2g - k1 * lap_c1 - k2 * lap_c2;
    m_2g = mu0_2 + log(c_2g) + l * (c_1g + c_sg) - k1 * lap_c2 - k2 * laplace(c_1g + c_sg);
    m_sg = mu0_s + log(c_sg) + l * c_2g - k1 * lap_cs - k2 * lap_c2;
    
    % Reaction currents
    j_1 = s_1 * (exp(mu_Y1) - exp(m_1g));
    j_2 = s_2 * (exp(2*m_1g + m_2g) - exp(3*m_1g));
    j_3 = s_3 * (exp(m_1g + mu_Y2) - exp(m_2g + mu_Y4));
    j_4 = s_4 * (exp(m_1g) - exp(mu_Y3));
    
    chem_term_1 = j_1 + j_2 - j_3 - j_4;
    chem_term_2 = -j_2 + j_3;
    
    % Current fluxes
    [m1_x, m1_y] = deal(gradient_x(m_1g), gradient_y(m_1g));
    [m2_x, m2_y] = deal(gradient_x(m_2g), gradient_y(m_2g));
    [ms_x, ms_y] = deal(gradient_x(m_sg), gradient_y(m_sg));
    
    [c1_x, c1_y] = deal(gradient_x(c_1g), gradient_y(c_1g));
    [c2_x, c2_y] = deal(gradient_x(c_2g), gradient_y(c_2g));
    [cs_x, cs_y] = deal(gradient_x(c_sg), gradient_y(c_sg));
    
    diff_term_1 = D_1 * (dot_grad(c1_x, c1_y, m1_x, m1_y) + c_1g .* laplace(m_1g));
    diff_term_2 = D_2 * (dot_grad(c2_x, c2_y, m2_x, m2_y) + c_2g .* laplace(m_2g));
    diff_term_s = D_s * (dot_grad(cs_x, cs_y, ms_x, ms_y) + c_sg .* laplace(m_sg));
    
    % Temporal derivative
    dc_1g = diff_term_1 + chem_term_1;
    dc_2g = diff_term_2 + chem_term_2;
    dc_sg = diff_term_s;
    
    if norm_chem
        m_1g = sqrt(gradient_x(m_1g).^2 + gradient_y(m_1g).^2);
        m_2g = sqrt(gradient_x(m_2g).^2 + gradient_y(m_2g).^2);
        m_sg = sqrt(gradient_x(m_sg).^2 + gradient_y(m_sg).^2);
    end
end
