clear; clc; close all;

% 유한차분
dx = 0.1; dy = 0.1; dz = 0.1; dt = 0.01;
Nx = 20; Ny = 20; Nz = 20; Nt = 100;
Lx = dx*Nx; Ly = dy*Ny; Lz = dz*Nz; Lt = dt*Nt;
x_grid = -dx:dx:Lx; y_grid = -dy:dy:Ly; z_grid = -dz:dz:Lz; t = dt:dt:Lt;
[X, Y, Z] = ndgrid(x_grid, y_grid, z_grid);

% 상수
R_u = 8.31446261815324; % J/(mol*K)
mu0 = 1.716e-5; % Pa*s
T0 = 273.15; % K
S = 110.4; % K
Pr = 0.71;
G = 6.67430e-11; % m^3/(kg*s^2)
m = 5.972e24; % kg
r0 = [-Lx/2 6.371e6 -Lz/2];
x_species = [0.78084, 0.20946, 0.00934];
x_species = x_species / sum(x_species);
M_species = [28.0134, 31.9988, 39.948] * 1e-3;
M_air = x_species * M_species';
R = R_u / M_air;
[A_air] = cal_A(x_species);

% 초기 설정값
u = zeros(Nx+2, Ny+2, Nz+2, 3);
T_initial = 303.0; % K
p_initial = 101325.0; % Pa

% 변수 초기화
T = ones(Nx+2, Ny+2, Nz+2) * T_initial;
p = ones(Nx+2, Ny+2, Nz+2) * p_initial;
rho = p ./ (R .* T);
[C_p, h] = air_properties(T, A_air, R);
e = h - R .* T;
E = rho .* (e + 0.5 .* sum(u.^2, 4));
gamma = C_p ./ (C_p - R);
mu = mu0 * (T / T0).^(3/2) .* (T0 + S) ./ (T + S);
kappa = zeros(Nx+2, Ny+2, Nz+2);
lambda = kappa - 2/3 * mu;
[dudx, dudy, dudz] = gradient(u, dx, dy, dz, 1);
J = cat(5, dudx, dudy, dudz);
divu = dudx(:,:,:,1) + dudy(:,:,:,2) + dudz(:,:,:,3);
I3 = reshape(eye(3), 1, 1, 1, 3, 3);
tau = mu .* (J + permute(J, [1 2 3 5 4])) + lambda .* divu .* I3;
k = mu .* C_p / Pr;
q = -k .* gradient(T, dx, dy, dz);
rx = X + r0(1);
ry = Y + r0(2);
rz = Z + r0(3);
r_norm = sqrt(rx.^2 + ry.^2 + rz.^2);
g = -G * m ./ r_norm.^3 .* cat(4, rx, ry, rz);

for i = 1:Nt
    [rho, u, E] = FDM(rho, u, E, p, tau, g, q, dx, dy, dz, dt);
    [T, p]
end

function [rho, u, E] = FDM(rho, u, E, p, tau, g, q, dx, dy, dz, dt)
    
end

function [A_air] = cal_A(x_species)
    % NASA 9 coefficients: N2, O2, Ar
    A = [
        2.210371497e4, -3.818461820e2, 6.082738360, -8.530914410e-3, 1.384646189e-5, -9.625793620e-9, 2.519705809e-12, 7.108460860e2, -1.076003744e1;
    -3.425563420e4,  4.847000970e2, 1.119010961,  4.293889240e-3, -6.836300520e-7, -2.023372700e-9, 1.039040018e-12, -3.391454870e3, 1.849699470e1;
        0, 0, 2.5, 0, 0, 0, 0, -7.453750000e2, 4.37967491
    ];

    % 혼합공기 NASA 계수 미리 계산
    A_air = x_species * A;
end
    
function [Cp, h] = air_properties(T, a, R)
    assert(all(T(:) >= 200 & T(:) <= 400), 'Temperature must be between 200 K and 400 K.');

    Cp_R = a(1)./T.^2 + a(2)./T + a(3) + a(4).*T + a(5).*T.^2 + a(6).*T.^3 + a(7).*T.^4;

    h_RT = -a(1)./T.^2 + a(2).*log(T)./T + a(3) + a(4).*T./2 + a(5).*T.^2./3 + a(6).*T.^3./4 + a(7).*T.^4./5 + a(8)./T;

    Cp = Cp_R .* R;
    h = h_RT .* R .* T;
end