clear; clc; close all;

% 유한차분
dx = 0.1; dy = 0.1; dz = 0.1; dt = 0.01;
Nx = 20; Ny = 20; Nz = 20; Nt = 100;
Lx = dx*Nx; Ly = dy*Ny; Lz = dz*Nz; Lt = dt*Nt;
x = dx:dx:Lx; y = dy:dy:Ly; z = dz:dz:Lz; t = dt:dt:Lt; % m

% 상수
R_u = 8.31446261815324; %J/(mol*K) 일반기체상수
M = 0.0289647; %kg/mol 공기분자량
R = R_u/M; %J/(kg*K) 기체상수
mu0 = 1.716e-5; %Pa*s 점성계수
T0 = 273.15; %K 기준온도
S = 110.4; %K Sutherland 상수
Pr = 0.71; %프란틀 수
G = 6.67430e-11; %m^3/(kg*s^2) 중력상수
m = 5.972e24; %kg 지구질량
r0 = [-Lx/2 6.371e6 -Lz/2]; %[m m m] 기준 위치벡터

% 초기 설정값
u = zeros(Nx+2, Ny+2, Nz+2, 3); % m/s
T_initial = 303.0;      % K
p_initial = 101325.0;   % Pa

% 변수 초기화
T = ones(Nx+2, Ny+2, Nz+2) * T_initial;
p = ones(Nx+2, Ny+2, Nz+2) * p_initial;
rho = p ./ (R * T); % kg/m^3
C_p = cp(T, p); % J/(kg*K)
gamma = C_p ./ (C_p - R);
E = p ./ (gamma - 1) + 0.5 * rho .* sum(u.^2, 4); % J/m^3
r = [0 0 0]; % m
g = -G * m * (r + r0) ./ vecnorm(r + r0, 2, 2).^3; % m/s^2

function Cp_R = cp(T)

    if any(T(:) < 200.0) || any(T(:) > 400.0)
        error('Temperature must be between 200 K and 400 K.');
    end

    % Molar masses [kg/mol]
    M_N2 = 28.0134e-3;
    M_O2 = 31.9988e-3;
    M_Ar = 39.948e-3;

    x_N2 = 0.78084;
    x_O2 = 0.20946;
    x_Ar = 0.00934;

    x_sum = x_N2 + x_O2 + x_Ar;

    x_N2 = x_N2 / x_sum;
    x_O2 = x_O2 / x_sum;
    x_Ar = x_Ar / x_sum;

    M_air = ...
          x_N2 * M_N2 ...
        + x_O2 * M_O2 ...
        + x_Ar * M_Ar;

    N2 = [ ...
         2.210371497e4, ...
        -3.818461820e2, ...
         6.082738360, ...
        -8.530914410e-3, ...
         1.384646189e-5, ...
        -9.625793620e-9, ...
         2.519705809e-12];

    O2 = [ ...
        -3.425563420e4, ...
         4.847000970e2, ...
         1.119010961, ...
         4.293889240e-3, ...
        -6.836300520e-7, ...
        -2.023372700e-9, ...
         1.039040018e-12];

    invT  = 1.0 ./ T;
    invT2 = invT .* invT;
    T2 = T .* T;
    T3 = T2 .* T;
    T4 = T2 .* T2;

    CpR_N2 = ...
          N2(1) .* invT2 ...
        + N2(2) .* invT ...
        + N2(3) ...
        + N2(4) .* T ...
        + N2(5) .* T2 ...
        + N2(6) .* T3 ...
        + N2(7) .* T4;

    CpR_O2 = ...
          O2(1) .* invT2 ...
        + O2(2) .* invT ...
        + O2(3) ...
        + O2(4) .* T ...
        + O2(5) .* T2 ...
        + O2(6) .* T3 ...
        + O2(7) .* T4;

    CpR_Ar = 2.5;

    Cp_molar = (x_N2 .* CpR_N2 + x_O2 .* CpR_O2 + x_Ar .* CpR_Ar);
    Cp = Cp_molar ./ M_air;
    Cp_R = Cp;
end