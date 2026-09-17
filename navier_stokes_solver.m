clear; clc; close all;

% 유한차분
dx = 0.1; dy = 0.1; dz = 0.1; dt = 1e-5;
Nx = 1e3; Ny = 1e3; Nz = 1e3; Nt = 1e5;
Lx = dx*Nx; Ly = dy*Ny; Lz = dz*Nz; Lt = dt*Nt;
x_grid = -dx:dx:Lx; y_grid = -dy:dy:Ly; z_grid = -dz:dz:Lz; t = dt:dt:Lt;
[X, Y, Z] = ndgrid(x_grid, y_grid, z_grid);

% 상수
mu0 = 1.716e-5; % Pa*s
T0 = 273.15; % K
S = 110.4; % K
Pr = 0.71;

G = 6.67430e-11; % m^3/(kg*s^2)
m = 5.972e24; % kg
r0 = [-Lx/2 6.371e6 -Lz/2];
rx = X + r0(1);
ry = Y + r0(2);
rz = Z + r0(3);
r_norm = sqrt(rx.^2 + ry.^2 + rz.^2);
g = -G * m ./ r_norm.^3 .* cat(4, rx, ry, rz);

x_species = [0.78084, 0.20946, 0.00934];
x_species = x_species / sum(x_species);
M_species = [28.0134, 31.9988, 39.948] * 1e-3;
R = 287.107; % J/(kg*K)
A_air = cal_A(x_species);

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
[tau, q] = update(u, T, A_air, R, Pr, mu0, T0, S, dx, dy, dz, Nx, Ny, Nz);

[T_table, e_table] = make_table_T(R, A_air);

for i = 1:Nt
    [rho2, u2, E2] = FDM(rho, rho .* u, u, E, p, tau, g, q, dx, dy, dz, dt, 1);
    [T, p] = update_T_p(rho2, E2, u2, R, T_table, e_table, A_air);
    [tau, q] = update(u2, T, A_air, R, Pr, mu0, T0, S, dx, dy, dz, Nx, Ny, Nz);
    [rho3, u3, E3] = FDM(rho2, rho2 .* u2, u2, E2, p, tau, g, q, dx, dy, dz, dt, -1);
    rho = 1/2 * (rho + rho3);
    u = 1/2 * (u + u3);
    E = 1/2 * (E + E3);
    [T, p] = update_T_p(rho, E, u, R, T_table, e_table, A_air);
    [rho, u, T, p, E] = apply_bc(rho, u, T, p, A_air, R);
    [tau, q] = update(u, T, A_air, R, Pr, mu0, T0, S, dx, dy, dz, Nx, Ny, Nz);
end

function [rho2, u2, E2] = FDM(rho, rhou, u, E, p, tau, g, q, dx, dy, dz, dt, s)
    rho2 = rho - dt * divergence(rhou, s, dx, dy, dz);
    rhouu = reshape(rhou, size(rhou,1), size(rhou,2), size(rhou,3), 3, 1) .* reshape(u, size(u,1), size(u,2), size(u,3), 1, 3);
    I3 = reshape(eye(3), 1, 1, 1, 3, 3);
    flux = tau - rhouu - p.*I3;
    u2 = 1./rho2 .* (rhou + dt * (divergence(flux, s, dx, dy, dz) + rho .* g));
    u5 = reshape(u, size(u,1), size(u,2), size(u,3), 1, 3);   % (.,.,.,1,j)
    tau_u = sum(tau .* u5, 5);                                 % (tau*u)_i = tau_{i,j} u_j
    flux = (E + p) .* u - tau_u + q; 
    E2 = E - dt * (divergence(flux, s, dx, dy, dz) +  rho .* sum(u .* g, 4));
end

function [T, p] = update_T_p(rho, E, u, R, T_table, e_table, A_air)
    e = E ./ rho - 0.5 .* sum(u.^2, 4);
    T =  interp1(e_table, T_table, e, 'linear');
    for iter = 1:3
        [C_p, ~] = air_properties(T, A_air, R);
        e_ = h - R .* T;
        C_v = C_p - R;
        T = T - (e_ - e) ./ C_v;
    end
    p = rho .* R .* T;
end

function [tau, q] = update(u, T, A_air, R, Pr, mu0, T0, S, dx, dy, dz, Nx, Ny, Nz)
    [C_p, ~] = air_properties(T, A_air, R);
    mu = mu0 * (T / T0).^(3/2) .* (T0 + S) ./ (T + S);
    kappa = zeros(Nx+2, Ny+2, Nz+2);
    lambda = kappa - 2/3 * mu;
    J = grad(u, dx, dy, dz);
    divu = J(:,:,:,1,1) + J(:,:,:,2,2) + J(:,:,:,3,3);
    I3 = reshape(eye(3), 1, 1, 1, 3, 3);
    tau = mu .* (J + permute(J, [1 2 3 5 4])) + lambda .* divu .* I3;
    k = mu .* C_p / Pr;
    q = -k .* grad(T, dx, dy, dz);
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

function [T_table, e_table] = make_table_T(R, a)
    T_table = 200:0.1:400;
    e_table = R *(-a(1)./T_table + a(2).*log(T_table) + (a(3)-1).*T_table + a(4).*T_table.^2./2 + a(5).*T_table.^3./3 + a(6).*T_table.^4./4 + a(7).*T_table.^5./5 + a(8));
end

function G = grad(A, dx, dy, dz)
    hs = [dx dy dz];
    d = ndims(A);
    G = zeros([size(A) 3]);
    idxOut = repmat({':'}, 1, d+1);
    for k = 1:3
        idxOut{d+1} = k;
        G(idxOut{:}) = centralDiff1(A, k, hs(k));
    end
end

function div = divergence(A, s, dx, dy, dz)
    assert(ndims(A) >= 4, 'divergence3: 4번째 차원이 성분축이어야 합니다.');
    hs = [dx dy dz];
    sz = size(A); sz(4) = [];
    idx = repmat({':'}, 1, ndims(A));
    div = zeros([sz 1]);
    for k = 1:3
        idx{4} = k;                                          % k번째 성분만
        Ak = reshape(A(idx{:}), [sz 1]);                      % 성분축 제거
        div = div + oneSidedDiff1(Ak, k, s, hs(k));           % k방향으로 미분
    end
end
function D = centralDiff1(A, k, h)
    d = ndims(A);
    n = size(A, k);
    D = zeros(size(A));
 
    im1 = repmat({':'}, 1, d); im1{k} = 1:n-2;    % i-1
    ip1 = im1;                 ip1{k} = 3:n;      % i+1
    ic  = im1;                 ic{k}  = 2:n-1;    % 내부(i)
    D(ic{:}) = (A(ip1{:}) - A(im1{:})) / (2*h);
 
    i0 = repmat({':'}, 1, d); i0{k} = 1;
    i1 = i0; i1{k} = 2;
    i2 = i0; i2{k} = 3;
    D(i0{:}) = (-3*A(i0{:}) + 4*A(i1{:}) - A(i2{:})) / (2*h);   % 앞쪽 경계
 
    jN  = repmat({':'}, 1, d); jN{k}  = n;
    jN1 = jN; jN1{k} = n-1;
    jN2 = jN; jN2{k} = n-2;
    D(jN{:}) = (3*A(jN{:}) - 4*A(jN1{:}) + A(jN2{:})) / (2*h);  % 뒤쪽 경계
end

function D = oneSidedDiff1(A, k, s, h)
    d = ndims(A);
    n = size(A, k);
    diffk = diff(A, 1, k);
    id = repmat({':'}, 1, d);
    if s > 0                                      % 후진차분: D(i)=A(i)-A(i-1)
        id{k} = 1;
        D = cat(k, diffk(id{:}), diffk);
    else                                           % 전진차분: D(i)=A(i+1)-A(i)
        id{k} = n-1;
        D = cat(k, diffk, diffk(id{:}));
    end
    D = D / h;
end

function [rho, u, T, p, E] = apply_bc(rho, u, T, p, A_air, R)
    % 격자: 1번째/끝 인덱스가 고스트층. y가 커지는 방향이 고도(위쪽).
    % 바닥(y=1) : 고정벽, no-slip     -> 속도 반사(u=0 @ 벽면), T/rho/p는 단열(zero-gradient)
    % 나머지 5면(x 양끝, y=끝(위), z 양끝) : 자유유출(zero-gradient, 안쪽 값 복제)
 
    % --- 자유유출: x 양끝 ---
    u(1,:,:,:)   = u(2,:,:,:);     T(1,:,:)   = T(2,:,:);     rho(1,:,:)   = rho(2,:,:);     p(1,:,:)   = p(2,:,:);
    u(end,:,:,:) = u(end-1,:,:,:); T(end,:,:) = T(end-1,:,:); rho(end,:,:) = rho(end-1,:,:); p(end,:,:) = p(end-1,:,:);
 
    % --- 자유유출: z 양끝 ---
    u(:,:,1,:)   = u(:,:,2,:);     T(:,:,1)   = T(:,:,2);     rho(:,:,1)   = rho(:,:,2);     p(:,:,1)   = p(:,:,2);
    u(:,:,end,:) = u(:,:,end-1,:); T(:,:,end) = T(:,:,end-1); rho(:,:,end) = rho(:,:,end-1); p(:,:,end) = p(:,:,end-1);
 
    % --- 자유유출: y 위쪽 끝 ---
    u(:,end,:,:) = u(:,end-1,:,:); T(:,end,:)  = T(:,end-1,:); rho(:,end,:) = rho(:,end-1,:); p(:,end,:) = p(:,end-1,:);
 
    % --- 바닥(y=1) : no-slip 고정벽 ---
    u(:,1,:,:) = -u(:,2,:,:);      % 벽면에서 속도(법선+접선 모두)가 0이 되도록 반사
    T(:,1,:)   = T(:,2,:);         % 단열벽 (zero-gradient)
    rho(:,1,:) = rho(:,2,:);       % 압력/밀도는 법선방향 zero-gradient
    p(:,1,:)   = p(:,2,:);
 
    % 강제로 맞춘 T,rho,u와 일관되도록 E 재계산
    [~, h] = air_properties(T, A_air, R);
    e = h - R .* T;
    E = rho .* (e + 0.5 .* sum(u.^2, 4));
end