%% ROUGH SOLUTIONS FOR BRUSSELLATOR

clear; clc; close all;

%% Parameters
cs    = 1; % non-reactive
muY1  = log(2);
muY2  = log(5);
muY3  = log(0.01);
muY4  = log(0.1);

% Grid
x_range = linspace(0.01,2,600);
y_range = linspace(0.01,6,600);
[X, Y] = meshgrid(x_range, y_range);

color1 = [255 221  95]/255;
color2 = [228 156 149]/255;
color3 = [109 114 209]/255;

%% l = chi
l_vals = linspace(0, 10, 2000);
tol = 5e-2; % Tolerance of solutions

%% Colormap
n = length(l_vals);
cmap = zeros(n, 3);
% First segment color1 --> color2
n1 = round(n/20);
n2 = round(n/5);
for k = 1:3
    cmap(1:n1,k) = linspace(color3(k), color2(k), n1);
end
% Second segment color2 --> color3
for k = 1:3
    cmap(n1+1:n2,k) = linspace(color2(k), color1(k), n2-n1);
end

for k = 1:3
    cmap(n2+1:n,k) = linspace(color1(k), color1(k), n-n2);
end

figure; hold on;

%% Searching for points that satisfies the homogeneous RD equation under the tolerance
for ind = 1:length(l_vals)
    l = l_vals(ind);

    % Chemical potentials
    mu1 = log(X) + l * Y;
    mu2 = log(Y) + l * (X + cs);

    % Currents
    j1 = exp(muY1) - exp(mu1);
    j2 = exp(2*mu1 + mu2) - exp(3*mu1);
    j3 = exp(mu1 + muY2) - exp(mu2 + muY4);
    j4 = exp(mu1) - exp(muY3);

    eq1 = j1 + j2 - j3 - j4;
    eq2 = -j2 + j3;

    % Rough satisfaction
    mask = (abs(eq1) < tol) & (abs(eq2) < tol);

    % Plots
    [yy, xx] = find(mask);
    if ~isempty(xx)
        scatter(x_range(xx), y_range(yy), 12, cmap(ind,:), 'filled', ...
            'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0);
    end
end

xlabel('$c_1$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$c_2$', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
colormap(cmap);
cbar = colorbar('TickLabelInterpreter', 'latex', 'Ticks', linspace(0, 1, 6), ...
    'TickLabels', arrayfun(@(x) sprintf('$%d$', x), linspace(0, 10, 6), 'UniformOutput', false));
cbar.Limits = [0 1];
ylabel(cbar, '$\chi$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
axis([0 1.2 0 6]);
hold off;