function navier_stokes_solver(C_p)
    clear; clc; close all;

    % 초깃값
    dx = 0.1; dy = 0.1; dz = 0.1; dt = 0.01;
    Nx = 21; Ny = 21; Nz = 21; Nt = 101;
    Lx = dx*(Nx-1); Ly = dy*(Ny-1); Lz = dz*(Nz-1); Lt = dt*(Nt-1);
    x = 0:dx:(Nx-1)*dx; y = 0:dy:(Ny-1)*dy; z = 0:dz:(Nz-1)*dz; t = 0:dt:(Nt-1)*dt;

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

    % 변수
    rho = ones(Nx,Ny,Nz) * 1.225; %kg/m^3 초기밀도
    u = zeros(Nx,Ny,Nz,3); %m/s 초기속도벡터
    p = ones(Nx,Ny,Nz) * 101325; %Pa 초기압력
    T = p./(rho*R); %K 초기온도

    C_p
    
end