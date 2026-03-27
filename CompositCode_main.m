clc, close all, clear all

%-----------

Data = [0.127	0.220	0.129	0.129	0.132	0.127	0.127		% Thickness
	40	210	136	151	75	147	181		% E1
	9.8	20	10	9.4	6	9	10.3		% E2
	2.8	6	5.2	4.8	2	3.3	7.17		% G12
	0.3	0.3	0.3	0.31	0.34	0.31	0.28		% nu12
	1100	1400	1800	2260	1400	2260	1500		% Sigma1t
	600	2800	1200	1200	280	1200	1500		% Sigma1c
	20	80	40	50	30	50	40		% Sigma2t
	140	280	220	190	140	190	246		% Sigma2c
	70	120	80	100	60	100	68		% Tau12
	0.028	0.007	0.013	0.015	0.019	0.015	0.0083		% Epsilon1t
	0.015	0.013	0.009	0.008	0.004	0.008	0.0083		% Epsilon1c
	0.002	0.004	0.004	0.005	0.005	0.005	0.0039		% Epsilon2t
	0.014	0.014	0.021	0.020	0.023	0.021	0.0239		% Epsilon2c
	0.014	0.020	0.015	0.022	0.030	0.030	0.0095];	% Gamma12
%----------------------------

% Coice of laminae and rotation of them
Laminate = [1 2 3]; % Choice of laminae for laminate (starts with the lower one)
theta = [45,-45, 0]; % The rotation of the respective laminae
N = [4,0,0]; % Normal force
M = [0,0,0]; % Bending moment 

% Calc. of nu21 That is not given in the data
nu21 = poisson(Data(5,Laminate),Data(2,Laminate),Data(3,Laminate));

% Creation of the important data from Data matrix
LaminateData = [Data([2,3,5,4,1],Laminate);nu21;theta]; % [E1,E2,nu12,G12,thickness,nu21,theta]' Column -> laminae
 
[A,B,D] = laminateStiffnessM(LaminateData) % A,B and D matrices CORRECT

[epsilon0, kappa] = strain_curvature(A,B,D,N,M) % Mid-surface strain and curvature (Row vectors)

epsilon = laminaStrains(epsilon0,kappa,LaminateData(5,:)) % Deformation for each laminae in global frame (1,2) Column -> laminae

sigma = stresses(LaminateData,epsilon) % Stresses in laminae local frame Column -> laminae

function sigmaM = stresses(data,epsilon)
% normal and shear stresses in the local systems for the different ply
    sigmaM = [];
    for i = 1:length(data(1,:))
        T = rotationM(data(7,i));
        Ql = stiffnessM([data(1:4,i);data(6,i)]); % [E1,E2,nu12,G12,nu21]'
        sigma = Ql*T'*epsilon(:,i); % A bit overcomplicated, could just do Q*epsilon in the global frame (1,2)
        sigmaM = [sigmaM,sigma];
    end

end


function epsilon = laminaStrains(epsilon0,kappa,thickness)
% Local strains in 1,2 coordinate system
h = sum(thickness);
epsilon = [];
z = -h/2;
for i = 1:length(thickness) 
    result = epsilon0 + (z+thickness(i)/2)*kappa;
    epsilon = [epsilon,result];
    z = z + thickness(i);
end
end


function nu21 = poisson(nu12,E1,E2)
nu21 = nu12.*E2./E1;
end


function [epsilon0, kappa] = strain_curvature(A,B,D,N,M) % deformation and curvature of laminate
% Calculates the Mid-surface strain and the curvature
result = [A,B;B,D]\[N';M']; 
epsilon0 = result(1:3);
kappa = result(4:6);
end


function [A,B,D] = laminateStiffnessM(data)
% Calculation of the Laminate Stiffness Matrices
% all input are data for each laninae 
thickness = data(5,:); % Extracts thickness from data
h = sum(thickness); % Whole thickness of laminate
z = - h/2;

A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

for i = 1:length(thickness)
    Qi = transStiffnessM([data(1:4,i);data(6:7,i)]); % [E1,E2,nu12,G12,nu21,theta]'
    A = A + Qi*thickness(i);
    B = B + Qi*((z+thickness(i))^2 -z^2)/2;
    D = D + Qi*((z+thickness(i))^3 - z^3)/3;
    z = z + thickness(i);
end
end


function Qi = transStiffnessM(data2)
% Calculation of transformed reduced stiffness matrix of one ply

theta = data2(6);
T = rotationM(theta);

Ql = stiffnessM(data2(1:5)); % [E1,E2,nu12,G12,nu21]'
Qi = T*Ql*T'; % Rotated stiffness matrix in 1,2 coordinate system
end


function Ql = stiffnessM(data3)
% Calculation of local stiffness matrix of one laminae before rotation
E1 = data3(1);
E2 = data3(2);
nu12 = data3(3);
G12 = data3(4);
nu21 = data3(5);

Ql = (1/(1-nu12*nu21))*[
    E1,nu21*E1,0;
    nu12*E2,E2,0
    0,0,G12
    ];
end

function T = rotationM(theta)
T = [
    cosd(theta)^2, sind(theta)^2, -2*sind(theta)*cosd(theta);
    sind(theta)^2, cosd(theta)^2, 2*sind(theta)*cosd(theta);
    sind(theta)*cosd(theta), -sind(theta)*cosd(theta), cosd(theta)^2-sind(theta)^2
    ]; % Rotation matrix for rotation of laminae
end
