%% INTEGRAL OF NON-REACTIVE SPECIES BRUSSELLATOR 

% This is the same code for every system and every steady state, it will
% change only the fileName

clear; clc; close all;

fileName = 'C:\Users\user\Desktop\UniPd_Thesis\Pb0\fig2\results_fig2\c_s_fig2';

% Skipping lines starting with '#'
opts = detectImportOptions(fileName, 'FileType', 'text', 'CommentStyle','#');
T = readtable(fileName, opts);

c_s0 = 1; % Homogeneous nr

x = T{:,1};
y = T{:,2};
c_s = T{:,3};
c_s = c_s - c_s0*ones(size(c_s));

xy = [x, y];

x_unique = unique(x);
y_unique = unique(y);

% Check if the grid is structured
if numel(x_unique) * numel(y_unique) == numel(c_s)
    c_s_grid = reshape(c_s, [numel(x_unique), numel(y_unique)]);
    
    % Mesh size
    dx = mean(diff(x_unique));
    dy = mean(diff(y_unique));

    % Compute the integral using the trapezoidal rule
    integral_c_s = trapz(y_unique, trapz(x_unique, c_s_grid, 1)) * dx * dy;

else
    % Delaunay for irregular data (convex hull)
    tri = delaunay(x, y);
    integral_c_s = 0;
    for k = 1:size(tri,1)
        % Vertices of the triangle
        idx = tri(k,:);
        xv = x(idx);
        yv = y(idx);
        fvals = c_s(idx);

        % Area of the triangle
        area = abs((xv(2)-xv(1))*(yv(3)-yv(1)) - (xv(3)-xv(1))*(yv(2)-yv(1))) / 2;

        % Average value of c_s over the triangle
        avg_c_s = mean(fvals);

        % Add to the total integral
        integral_c_s = integral_c_s + avg_c_s * area;
    end
end

fprintf('The integral of c_s - c_s0 over the (x,y) space is: %.6f\n', integral_c_s);