%% LINEAR STABILITY ANALYSIS BRUSSELLATOR

% To change value of chi, e.g. from 1 to 2, use Ctrl+F and substitute ALL chi=1 with chi=2

clear; clc; close all;

%% Parameters
n_x = 3; % X
n_r = 4; % Reactions
n_y = 4; % Y

mu0 = [0; 0; 0]; % μ⁰
s = ones(n_r, 1); % s
mu_star = [log(2); log(5); log(0.01); log(0.1)]; % μ* 
D = [0.5; 10; 1]; % D

% L matrix
chi=1;
L = [0, chi, 0;
     chi, 0, chi;
     0, chi, 0];

% K matrix
k1 = 1;
k2 = 0.5;
K = [k1, k2, 0;
     k2, k1, k2;
     0, k2, k1];

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

scriptFullPath = mfilename('fullpath'); 
[scriptDir, ~, ~] = fileparts(scriptFullPath);

% Load homogeneous steady states
ss_file_path = fullfile(scriptDir, 'ss_fig2_chi=1.txt');
c_values = load(ss_file_path);

% q range
q_values = linspace(0, 5, 2000);

% Saving files
plots_folder = fullfile(scriptDir, 'plots_fig2_chi=1');
if ~exist(plots_folder, 'dir')
    mkdir(plots_folder);
end

q0_file_path = fullfile(scriptDir, 'q0_fig2.txt_chi=1');
fileID = fopen(q0_file_path, 'w');

image_tol = [5e1; 5e2; 0];

% Loop over steady states
for row_idx = 1:size(c_values, 1)
    c = c_values(row_idx, :)';
    
    % Chemical potentials
    mu = mu0 + log(c) + L * c;
    
    % Eigenvalues
    lambda = zeros(size(q_values));
    lambda_d = zeros(size(q_values));
    
    % Determinants
    detM = zeros(size(q_values));
    detB = zeros(size(q_values));
    detBM = zeros(size(q_values));
    
    % Imaginary eigenvalues
    im = 0;
    im_d = 0;

    % Loop over q values
    for idx = 1:length(q_values)
        q = q_values(idx);
        
        % A matrix
        A = diag(D .* c);
        
        % C matrix
        C = zeros(n_x, n_x);
        for i = 1:n_x
            for j = 1:n_x 
                sum_rho = 0;
                for rho = 1:n_r
                    term1 = nu(j, rho) * exp(sum(mu .* nu(:, rho)) + sum(mu_star .* nu_star(:, rho)));
                    term2 = bar_nu(j, rho) * exp(sum(mu .* bar_nu(:, rho)) + sum(mu_star .* bar_nu_star(:, rho)));
                    sum_rho = sum_rho + S(i, rho) * s(rho) * (term1 - term2);
                end
                C(i, j) = sum_rho;
            end
        end
        
        % B matrix
        B = -q^2 * A + C;
        
        % B_d matrix
        B_d = -q^2 * A;
        
        % M matrix
        inv_c = [1/c(1); 1/c(2); 1/c(3)];
        M = L + q^2*K + diag(inv_c);
        
        % Maximum real eigenvalue (BM)
        vals = eig(B*M);
        lambda(idx) = vals(1);
        for j = 1:length(vals)
            if imag(vals(j)) ~= 0
                im = im + 1;      
                if abs(lambda(idx)) > real(vals(j))
                    lambda(idx) = real(vals(j));
                end
            else
                if lambda(idx) < vals(j)
                    lambda(idx) = vals(j);
                end
            end
        end

        % Maximum real eigenvalue (B_dM)
        vals_d = eig(B_d*M);
        lambda_d(idx) = vals_d(1);
        for j = 1:length(vals_d)
            if imag(vals_d(j)) ~= 0
                im_d = im_d + 1;      
                if abs(lambda_d(idx)) > real(vals_d(j))
                    lambda_d(idx) = real(vals_d(j));
                end
            else
                if lambda_d(idx) < vals_d(j)
                    lambda_d(idx) = vals_d(j);
                end
            end
        end
        
        % Determinants
        detM(idx)   = prod(eig(M));
        detB(idx)   = prod(eig(B));
        detBM(idx)  = prod(eig(B*M)); 
    end
    
    % First q0 where lambda is negative
    q0_idx = find(lambda < 0, 1, 'first');
    q0 = q_values(q0_idx);
    
    % First q0_d where lambda_d is negative
    q0_d_idx = find(lambda_d < 0, 1, 'first');
    q0_d = q_values(q0_d_idx);
    
    fprintf(fileID, '%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', c(1), c(2), c(3), q0, q0_d);
    
    fprintf('For c = [%.4f, %.4f, %.4f]:\n', c(1), c(2), c(3));
    fprintf('First q0 where lambda is negative: %.4f\n', q0);
    fprintf('First q0_d where lambda_d is negative: %.4f\n', q0_d);
    fprintf('Percentage of imaginary eigenvalues lambda: %.2f%%\n', im/length(q_values)*100/3);
    fprintf('Percentage of imaginary eigenvalues lambda_d: %.2f%%\n', im_d/length(q_values)*100/3);
    
    color1 = [203 219 055]/255;
    color2 = [228 156 149]/255;
    color3 = [146 171 211]/255;
    
    % Eigenvalue vs. q
    figure('Visible','on');
    plot(q_values, lambda, '-', 'Color', color3, 'LineWidth', 1.8);
    hold on;
    plot(q_values, lambda_d, '-', 'Color', color1, 'LineWidth', 1.8);
    plot(q_values, zeros(size(q_values)), '--', 'Color', color2, 'LineWidth', 1.8);
    xlabel('$q$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Eigenvalue', 'Interpreter', 'latex', 'FontSize', 14);
    legend('$\lambda(q)$', '$\lambda_d(q)$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex', 'FontSize', 14);
    hold off
    grid on;
    
    % Determinants vs. q (semilogy)
    figure('Visible','off');
    semilogy(q_values, abs(detM), '-', 'Color', color1, 'LineWidth', 1.5);
    hold on;
    semilogy(q_values, abs(detB), '-', 'Color', color2, 'LineWidth', 1.5);
    semilogy(q_values, abs(detBM), '-', 'Color', color3, 'LineWidth', 1.5);
    xlabel('q');
    ylabel('|det|');
    legend('|det M|', '|det B|', '|det BM|');
    grid on;
    hold off;  
end

fclose(fileID);