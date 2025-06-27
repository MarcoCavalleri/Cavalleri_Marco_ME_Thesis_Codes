%% HOMOGENEOUS STEADY STATES BRUSSELLATOR

% To change value of chi, e.g. from 1 to 2, use Ctrl+F and substitute ALL chi=1 with chi=2

clear; clc; close all;

%% Parameters
n_x = 3; % X
n_r = 4; % Reactions
n_y = 4; % Y

mu0 = [0; 0; 0]; % μ⁰
s = ones(n_r, 1); % s
mu_star = [log(2); log(5); log(0.01); log(0.1)]; % μ* 

% L matrix
chi = 1;
L = [0, chi, 0;
     chi, 0, chi;
     0, chi, 0];

% X
nu = [1, 3, 0, 0;
      0, 0, 1, 0;
      0, 0, 0, 0]; % ν

bar_nu = [0, 2, 1, 1;
          0, 1, 0, 0;
          0, 0, 0, 0]; % ν̄ 

% Y
nu_star = [0, 0, 0, 0;
           0, 0, 0, 0;
           0, 0, 0, 1;
           0, 0, 1, 0]; % ν*

bar_nu_star = [1, 0, 0, 0;
               0, 0, 1, 0; 
               0, 0, 0, 0;
               0, 0, 0, 0]; % ν̄*

S = bar_nu - nu;

% Pack all parameters
params.n_x = n_x;
params.n_r = n_r;
params.n_y = n_y;
params.mu0 = mu0;
params.s = s;
params.mu_star = mu_star;
params.L = L;
params.S = S;
params.nu = nu;
params.nu_star = nu_star;
params.bar_nu = bar_nu;
params.bar_nu_star = bar_nu_star;

%% Newton method

% Search ranges: concentrations are labeled by x (different for the two species)
range_x1 = 0.1:0.1:1.5;
range_x2 = 0.1:0.1:5.1;

[x1_comb, x2_comb] = ndgrid(range_x1, range_x2);

x_combinations = [x1_comb(:), x2_comb(:)];

% Homogeneous steady states
ss = [];

fprintf('Starting Newton iterations...\n');

total_combinations = size(x_combinations, 1);

for idx = 1:total_combinations
    x = x_combinations(idx, :)'; % initial guess
    x_nr = 1; % fixed value for non reactive
    
    max_iter = 50; % Newton iterations
    tol = 1e-8; % Tolerance
    
    control = false;
    for iter = 1:max_iter
        % g and its Jacobian J at current x (H is labeled as J)
        [g_val, J_val] = compute_g_and_J(x, x_nr, params);
        
        % Norm of the residual
        norm_g = norm(g_val);
        
        if norm_g < tol
            control = true;
            break;
        end
        
        %  J * delta = -g
        if rcond(J_val) < 1e-12 % J singular
            break;
        end
        delta = - J_val \ g_val(1:n_x-1);
        
        % Ensure concentrations positivity
        step = 1.0;
        while any(x + step*delta <= 0)
            step = step / 2;
            if step < 1e-8
                break; % Skip if too small
            end
        end
        
        % Update
        x = x + step * delta;
    end

    x = real(x);

    if control
        % Check positivity and uniqueness
        if all(x > 0)
            unique = true;
            for s = 1:size(ss, 1)
                diff = abs(ss(s, 1:n_x-1) - x');
                if all(diff < tol * 1e4)
                    unique = false;
                    break;
                end
            end
    
            if unique
                ss = [ss; x', x_nr];
            end
        end
    end

    fprintf('Progress: %.2f%%\n', (idx / total_combinations) * 100);
end

fprintf('Unique solutions found:\n');
disp(ss);

scriptFullPath = mfilename('fullpath'); 
[scriptDir, ~, ~] = fileparts(scriptFullPath);

filePath = fullfile(scriptDir, 'ss_fig2_chi=1.txt');

fileID = fopen(filePath, 'w');
fprintf(fileID, '%.6f\t%.6f\t%.6f\t%.6f\n', ss');
fclose(fileID);

%% Function to compute g and Jacobian J
function [g, J] = compute_g_and_J(c, x_nr, params)
    % Unpack parameters
    n_x = params.n_x;
    n_r = params.n_r;
    n_y = params.n_y;
    mu0 = params.mu0;
    s = params.s;
    mu_star = params.mu_star;
    L = params.L;
    S = params.S;
    nu = params.nu;
    nu_star = params.nu_star;
    bar_nu = params.bar_nu;
    bar_nu_star = params.bar_nu_star;
    
    dx = 1e-8; % Finite difference step size
    
    c_full = [c; x_nr];
    
    g = comp_g(c_full, params);
    J = zeros(n_x - 1, n_x - 1);
    
    for i = 1:(n_x - 1)
        c_perturbed = c;
        c_perturbed(i) = c_perturbed(i) + dx;
        c_full_perturbed = [c_perturbed; x_nr];
        g_perturbed = comp_g(c_full_perturbed, params);
        J(:, i) = (g_perturbed(1:n_x-1) - g(1:n_x-1)) / dx; % Finite differences
    end
end

function g = comp_g(c, params)
    n_x = params.n_x;
    n_r = params.n_r;
    n_y = params.n_y;
    mu0 = params.mu0;
    s = params.s;
    mu_star = params.mu_star;
    L = params.L;
    S = params.S;
    nu = params.nu;
    nu_star = params.nu_star;
    bar_nu = params.bar_nu;
    bar_nu_star = params.bar_nu_star;

    % Chemical potentials
    mu = mu0 + log(c) + L * c; 
    
    % Reaction terms (no diffusivity fo homogeneous)
    Exp_pl = zeros(n_r, 1);
    Exp_mi = zeros(n_r, 1);
    
    for r = 1:n_r
        Exp_pl(r) = sum(mu_star .* nu_star(:, r)) + sum(mu .* nu(:, r));
        Exp_mi(r) = sum(mu_star .* bar_nu_star(:, r)) + sum(mu .* bar_nu(:, r));
    end
    
    g = zeros(n_x, 1);
    
    for x = 1:n_x
        for r = 1:n_r
            g(x) = g(x) + S(x, r) * s(r) * (exp(Exp_pl(r)) - exp(Exp_mi(r)));
        end
    end
end