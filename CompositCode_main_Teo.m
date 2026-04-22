%==========================================================================
% Laminate Composite Analysis
%==========================================================================
clc; close all; clear;

%--------------------------------------------------------------------------
% Material Properties Data
%--------------------------------------------------------------------------
% Each column corresponds to a different material type.
% Lamina    1           2           3             4             5             6           7
% Fibre     E-glass     Boron       Carbon (HT)   Carbon (IM)   Kevlar 49     Carbon      Carbon
Data = [
            0.127       0.220       0.129         0.129         0.132         0.127       0.127;    % 1: Thickness [mm]
            40000       210000      136000        151000        75000         147000      181000;   % 2: E1 [MPa]
            9800        20000       10000         9400          6000          9000        10300;    % 3: E2 [MPa]
            2800        6000        5200          4800          2000          3300        7170;     % 4: G12 [MPa]
            0.3         0.3         0.3           0.31          0.34          0.31        0.28;     % 5: nu12
            1100        1400        1800          2260          1400          2260        1500;     % 6: Sigma1t [MPa]
            600         2800        1200          1200          280           1200        1500;     % 7: Sigma1c [MPa]
            20          80          40            50            30            50          40;       % 8: Sigma2t [MPa]
            140         280         220           190           140           190         246;      % 9: Sigma2c [MPa]
            70          120         80            100           60            100         68;       % 10: Tau12 [MPa]
            0.028       0.007       0.013         0.015         0.019         0.015       0.0083;   % 11: Epsilon1t
            0.015       0.013       0.009         0.008         0.004         0.008       0.0083;   % 12: Epsilon1c
            0.002       0.004       0.004         0.005         0.005         0.005       0.0039;   % 13: Epsilon2t
            0.014       0.014       0.021         0.020         0.023         0.021       0.0239;   % 14: Epsilon2c
            0.014       0.020       0.015         0.022         0.030         0.030       0.0095    % 15: Gamma12
];

%--------------------------------------------------------------------------
% Laminate Setup & Loads
%--------------------------------------------------------------------------
Laminate = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3];    % Choice of lamina materials (starts from bottom)
theta    = [0, 0, 90, 45, -45, -45, 45, 90, 0, 0]; % [02 / 90 / +45 / -45]_s [deg]
N        = [100, 100, 50];    % Normal forces [Nx, Ny, Nxy]
M        = [0, 0, 0];    % Bending moments [Mx, My, Mxy]

% Calculate nu21 (not provided directly in Data matrix)
nu21 = poisson(Data(5, Laminate), Data(2, Laminate), Data(3, Laminate));

% Organize laminate data (each column is a lamina):
% Rows -> [E1, E2, nu12, G12, thickness, nu21, theta]'
LaminateData = [Data([2, 3, 5, 4, 1], Laminate); nu21; theta]; 

%--------------------------------------------------------------------------
% Core Calculations
%--------------------------------------------------------------------------
% 1. Laminate Stiffness Matrices
[A, B, D] = laminateStiffnessM(LaminateData)

% 2. Mid-surface strain and curvature
[epsilon0, kappa] = strain_curvature(A, B, D, N, M)

% 3. Deformation for each lamina in global frame (x,y) at bottom and top interfaces
[epsilon_bot, epsilon_top] = laminaStrains(epsilon0, kappa, LaminateData(5, :))

% 4. Stresses in lamina local frame (1,2) at bottom and top interfaces
[sigma_bot, sigma_top] = stresses(LaminateData, epsilon_bot, epsilon_top)

% 5. Stress ratios (f) for each lamina (maximum absolute value of top and bottom)
f_bot = stressRatios(sigma_bot, Data(6:10, Laminate));
f_top = stressRatios(sigma_top, Data(6:10, Laminate));
f = maxAbs(f_bot, f_top);

% 6. Strain ratios (g) for each lamina (maximum absolute value of top and bottom)
g = strainRatios(LaminateData, epsilon_bot, epsilon_top, Data(11:15, Laminate));

% 7. Fracture Criteria Index (IF) per lamina (IF >= 1 means failure)
IF_MaxStress = max(abs(f), [], 1);
IF_MaxStrain = max(abs(g), [], 1);

th_bot = failureTsaiHill(sigma_bot, Data(6:10, Laminate));
th_top = failureTsaiHill(sigma_top, Data(6:10, Laminate));
IF_TsaiHill = max(th_bot, th_top);

% Display results
fprintf('\n--- Failure Analysis Results ---\n');
fprintf('Lamina | Max Stress IF | Max Strain IF | Tsai-Hill IF\n');
fprintf('-------|---------------|---------------|-------------\n');
for i = 1:length(Laminate)
    fprintf('  %d    |     %.4f    |     %.4f    |    %.4f\n', ...
        i, IF_MaxStress(i), IF_MaxStrain(i), IF_TsaiHill(i));
end

% 8. Homogenized Elastic Constants (Ex, Gxy)
h = sum(LaminateData(5, :)); % Total thickness
a_matrix = inv(A);           % Compliance matrix (inverse of A)
Ex = 1 / (h * a_matrix(1, 1));
Gxy = 1 / (h * a_matrix(3, 3));

fprintf('\n--- Homogenized Elastic Constants ---\n');
fprintf('Ex:  %.2f MPa\n', Ex);
fprintf('Gxy: %.2f MPa\n', Gxy);

fprintf('\n--- D Matrix (Bending Stiffness) [N*mm] ---\n');
disp(D);


%==========================================================================
% Helper Functions
%==========================================================================

function g = strainRatios(data, eps_bot, eps_top, limit_strains)
    % Calculates the ratio of current strain to maximum allowable strain
    % limit_strains format: [Epsilon1t; Epsilon1c; Epsilon2t; Epsilon2c; Gamma12]
    
    num_laminae = size(data, 2);
    g = zeros(3, num_laminae); 
    
    for i = 1:num_laminae
        T = rotationM(data(7, i));
        
        % Transform global strains to local 1-2 frame
        eps12_bot = T' * eps_bot(:, i);
        eps12_top = T' * eps_top(:, i);
        
        % Compute ratios for both interfaces
        r_bot = calcRatio(eps12_bot, limit_strains(:, i));
        r_top = calcRatio(eps12_top, limit_strains(:, i));
        
        % Take the worst-case ratio (by absolute magnitude)
        g(:, i) = maxAbs(r_bot, r_top);
    end
end

function r = calcRatio(eps, limits)
    r = zeros(3, 1);
    if eps(1) >= 0; r(1) = eps(1) / limits(1); else; r(1) = eps(1) / limits(2); end
    if eps(2) >= 0; r(2) = eps(2) / limits(3); else; r(2) = eps(2) / limits(4); end
    r(3) = eps(3) / limits(5);
end

function f = stressRatios(sigma, strengths)
    % Calculates the ratio of current stress to maximum allowable stress
    % strengths format: [Sigma1t; Sigma1c; Sigma2t; Sigma2c; Tau12]
    
    num_laminae = size(sigma, 2);
    f = zeros(3, num_laminae); 
    
    for i = 1:num_laminae
        % Direction 1
        if sigma(1, i) >= 0 % Tensile
            f(1, i) = sigma(1, i) / strengths(1, i);
        else                % Compressive
            f(1, i) = sigma(1, i) / strengths(2, i);
        end
        
        % Direction 2
        if sigma(2, i) >= 0 % Tensile
            f(2, i) = sigma(2, i) / strengths(3, i);
        else                % Compressive
            f(2, i) = sigma(2, i) / strengths(4, i);
        end
        
        % Shear direction
        f(3, i) = sigma(3, i) / strengths(5, i);
    end
end

function [sigma_bot, sigma_top] = stresses(data, eps_bot, eps_top)
    % Calculates normal and shear stresses in the local systems for bottom and top interfaces
    num_laminae = size(data, 2);
    sigma_bot = zeros(3, num_laminae); 
    sigma_top = zeros(3, num_laminae);
    
    for i = 1:num_laminae
        T = rotationM(data(7, i));
        
        % Extract properties: [E1, E2, nu12, G12, nu21]'
        properties = [data(1:4, i); data(6, i)]; 
        Ql = stiffnessM(properties); 
        
        % Stress in local frame
        sigma_bot(:, i) = Ql * T' * eps_bot(:, i);
        sigma_top(:, i) = Ql * T' * eps_top(:, i);
    end
end

function [epsilon_bot, epsilon_top] = laminaStrains(epsilon0, kappa, thickness)
    % Calculates strains in laminate axis system (x,y) at bottom and top interfaces of each ply
    num_laminae = length(thickness);
    epsilon_bot = zeros(3, num_laminae); 
    epsilon_top = zeros(3, num_laminae);
    
    h = sum(thickness);
    z = -h / 2;
    
    for i = 1:num_laminae 
        epsilon_bot(:, i) = epsilon0 + z * kappa;
        z_next = z + thickness(i);
        epsilon_top(:, i) = epsilon0 + z_next * kappa;
        
        z = z_next;
    end
end

function nu21 = poisson(nu12, E1, E2)
    % Calculates minor Poisson's ratio
    nu21 = nu12 .* E2 ./ E1;
end

function [epsilon0, kappa] = strain_curvature(A, B, D, N, M) 
    % Calculates the mid-surface strain (epsilon0) and curvature (kappa)
    ABD = [A, B; B, D];
    loads = [N'; M'];
    
    result = ABD \ loads; 
    
    epsilon0 = result(1:3);
    kappa = result(4:6);
end

function [A, B, D] = laminateStiffnessM(data)
    % Calculates the Laminate Stiffness Matrices (A, B, D)
    thickness = data(5, :); 
    
    A = zeros(3, 3);
    B = zeros(3, 3);
    D = zeros(3, 3);
    
    h = sum(thickness); 
    z = -h / 2;
    
    for i = 1:length(thickness)
        % Extract properties: [E1, E2, nu12, G12, nu21, theta]'
        properties = [data(1:4, i); data(6:7, i)];
        Qi = transStiffnessM(properties); 
        
        t = thickness(i);
        z_next = z + t;
        
        A = A + Qi * t;
        B = B + Qi * (z_next^2 - z^2) / 2;
        D = D + Qi * (z_next^3 - z^3) / 3;
        
        z = z_next;
    end
end

function Qi = transStiffnessM(data)
    % Calculates transformed reduced stiffness matrix of one lamina
    theta = data(6);
    T = rotationM(theta);
    
    % properties: [E1, E2, nu12, G12, nu21]'
    Ql = stiffnessM(data(1:5)); 
    
    % Rotated stiffness matrix in global (x,y) coordinate system
    Qi = T * Ql * T'; 
end

function Ql = stiffnessM(data)
    % Calculates local stiffness matrix (Q) of one lamina before rotation
    E1   = data(1);
    E2   = data(2);
    nu12 = data(3);
    G12  = data(4);
    nu21 = data(5);
    
    nu_factor = 1 / (1 - nu12 * nu21);
    
    Ql = nu_factor * [
        E1,        nu21 * E1, 0;
        nu12 * E2, E2,        0;
        0,         0,         G12
    ];
end

function T = rotationM(theta)
    % Calculates rotation matrix for stress/strain transformation
    c = cosd(theta);
    s = sind(theta);
    
    T = [
        c^2,   s^2,  -2 * s * c;
        s^2,   c^2,   2 * s * c;
        s * c, -s * c,  c^2 - s^2
    ]; 
end

function C = maxAbs(A, B)
    % Returns the element from A or B that has the largest absolute value, retaining its original sign.
    idx = abs(A) > abs(B);
    C = B;
    C(idx) = A(idx);
end

function th = failureTsaiHill(sigma, strengths)
    % Calculates the Tsai-Hill failure index for each lamina interface
    % strengths format: [Sigma1t; Sigma1c; Sigma2t; Sigma2c; Tau12]
    
    num_laminae = size(sigma, 2);
    th = zeros(1, num_laminae);
    
    for i = 1:num_laminae
        s1 = sigma(1, i);
        s2 = sigma(2, i);
        t12 = sigma(3, i);
        
        % Select strengths based on stress sign (tension/compression)
        if s1 >= 0; X = strengths(1, i); else; X = strengths(2, i); end
        if s2 >= 0; Y = strengths(3, i); else; Y = strengths(4, i); end
        S = strengths(5, i);
        
        % Tsai-Hill criterion formula
        th(i) = (s1/X)^2 - (s1*s2)/(X^2) + (s2/Y)^2 + (t12/S)^2;
    end
end
