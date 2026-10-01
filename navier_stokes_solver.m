clear; clc; close all;

% 유한차분
dx = 0.1; dy = 0.1; dz = 0.1; dt = 1e-5;
Nx = 1e2; Ny = 1e2; Nz = 1e2; Nt = 1e3;
Lx = dx*Nx; Ly = dy*Ny; Lz = dz*Nz;
x_grid = -dx:dx:Lx; y_grid = -dy:dy:Ly; z_grid = -dz:dz:Lz;
n_report = 100; % 진행 상황 출력 간격 (스텝)

% 상수
par.h = [dx dy dz];
par.dt = dt;
par.mu0 = 1.716e-5; % Pa*s
par.T0 = 273.15; % K
par.S = 110.4; % K
par.Pr = 0.71;

G = 6.67430e-11; % m^3/(kg*s^2)
m = 5.972e24; % kg
r_earth = 6.371e6; % m
par.g = reshape([0, -G*m/r_earth^2, 0], 1, 1, 1, 3); % 영역이 작아 균일하다고 본다

x_species = [0.78084, 0.20946, 0.00934];
x_species = x_species / sum(x_species);
M_species = [28.0134, 31.9988, 39.948] * 1e-3; % kg/mol
par.R = 8.314462618 / (x_species * M_species.'); % J/(kg*K)
par.a = cal_A(x_species);

% e -> T 역변환용 테이블 (NASA 계수의 유효 범위)
par.T_table = 200:0.1:1000;
[~, h_table] = air_properties(par.T_table, par.a, par.R);
par.e_table = h_table - par.R .* par.T_table;

% 초기 설정값
T_initial = 303.0; % K
p_initial = 101325.0; % Pa, y = 0에서의 압력

% 변수 초기화: 정지 상태, 정역학 평형(dp/dy = rho*g_y)
u = zeros(Nx+2, Ny+2, Nz+2, 3);
T = ones(Nx+2, Ny+2, Nz+2) * T_initial;
p = ones(Nx+2, Ny+2, Nz+2) .* p_initial .* exp(par.g(2) .* y_grid ./ (par.R * T_initial));
rho = p ./ (par.R .* T);
[~, h] = air_properties(T, par.a, par.R);
E = rho .* (h - par.R .* T + 0.5 .* sum(u.^2, 4));
[rho, u, E] = apply_bc(rho, u, E, par);

for i = 1:Nt
    % 예측 단계: 후진 차분
    [rho1, rhou1, E1] = FDM(rho, u, E, par, 1);
    [rho1, u1, E1] = apply_bc(rho1, rhou1 ./ rho1, E1, par);

    % 수정 단계: 전진 차분 후 보존량(rho, rho*u, E)을 평균
    [rho2, rhou2, E2] = FDM(rho1, u1, E1, par, -1);
    rhou = 1/2 * (rho .* u + rhou2);
    rho = 1/2 * (rho + rho2);
    E = 1/2 * (E + E2);
    [rho, u, E] = apply_bc(rho, rhou ./ rho, E, par);

    if mod(i, n_report) == 0 || i == Nt
        [T, p] = update_T_p(rho, u, E, par);
        speed = sqrt(sum(u.^2, 4));
        fprintf('step %d/%d  t = %.3e s  max|u| = %.3e m/s  T = [%.3f, %.3f] K\n', ...
            i, Nt, i*dt, max(speed(:)), min(T(:)), max(T(:)));
    end
end

save('navier_stokes_result.mat', 'x_grid', 'y_grid', 'z_grid', 'rho', 'u', 'T', 'p');

% z 중앙 단면
kz = round((Nz+2)/2);
figure;
subplot(1, 2, 1);
imagesc(x_grid, y_grid, T(:,:,kz).'); axis xy equal tight; colorbar;
xlabel('x [m]'); ylabel('y [m]'); title('T [K]');
subplot(1, 2, 2);
imagesc(x_grid, y_grid, speed(:,:,kz).'); axis xy equal tight; colorbar;
xlabel('x [m]'); ylabel('y [m]'); title('|u| [m/s]');

function [rho2, rhou2, E2] = FDM(rho, u, E, par, s)
    % 보존량 (rho, rho*u, E)를 dt만큼 전진. s > 0 후진 차분, s < 0 전진 차분
    [T, p] = update_T_p(rho, u, E, par);
    [tau, q] = viscous_flux(u, T, par);

    rhou = rho .* u;
    u_j = permute(u, [1 2 3 5 4]);                  % (.,.,.,1,j)
    I3 = reshape(eye(3), 1, 1, 1, 3, 3);
    flux_rhou = rhou .* u_j + p .* I3 - tau;        % rho*u_i*u_j + p*delta_ij - tau_ij
    flux_E = (E + p) .* u - sum(tau .* u_j, 5) + q; % (tau*u)_i = tau_{i,j} u_j

    rho2 = rho - par.dt * div_onesided(rhou, s, par.h);
    rhou2 = rhou - par.dt * (div_onesided(flux_rhou, s, par.h) - rho .* par.g);
    E2 = E - par.dt * (div_onesided(flux_E, s, par.h) - rho .* sum(u .* par.g, 4));
end

function [T, p] = update_T_p(rho, u, E, par)
    e = E ./ rho - 0.5 .* sum(u.^2, 4);
    T = temperature(e, par);
    p = rho .* par.R .* T;
end

function T = temperature(e, par)
    % 테이블 보간으로 초기값을 잡고 Newton 반복으로 e(T) = e를 푼다
    T = interp1(par.e_table, par.T_table, e, 'linear');
    if any(isnan(T), 'all')
        error('navier_stokes_solver:temperatureRange', ...
            '온도가 %g~%g K 범위를 벗어났습니다. 계산이 발산했을 수 있습니다.', ...
            par.T_table(1), par.T_table(end));
    end
    for iter = 1:3
        [C_p, h] = air_properties(T, par.a, par.R);
        T = T - (h - par.R .* T - e) ./ (C_p - par.R);
    end
end

function [tau, q] = viscous_flux(u, T, par)
    [C_p, ~] = air_properties(T, par.a, par.R);
    mu = par.mu0 * (T / par.T0).^(3/2) .* (par.T0 + par.S) ./ (T + par.S);
    lambda = -2/3 * mu; % 체적점성 0 (Stokes 가정)
    J = grad(u, par.h);
    divu = J(:,:,:,1,1) + J(:,:,:,2,2) + J(:,:,:,3,3);
    I3 = reshape(eye(3), 1, 1, 1, 3, 3);
    tau = mu .* (J + permute(J, [1 2 3 5 4])) + lambda .* divu .* I3;
    k = mu .* C_p / par.Pr;
    q = -k .* grad(T, par.h);
end

function [A_air] = cal_A(x_species)
    % NASA 9 coefficients (200~1000 K): N2, O2, Ar
    A = [
         2.210371497e4, -3.818461820e2, 6.082738360, -8.530914410e-3,  1.384646189e-5, -9.625793620e-9, 2.519705809e-12,  7.108460860e2, -1.076003744e1;
        -3.425563420e4,  4.847000970e2, 1.119010961,  4.293889240e-3, -6.836300520e-7, -2.023372700e-9, 1.039040018e-12, -3.391454870e3,  1.849699470e1;
         0, 0, 2.5, 0, 0, 0, 0, -7.453750000e2, 4.37967491
    ];

    % 혼합공기 NASA 계수 미리 계산
    A_air = x_species * A;
end

function [Cp, h] = air_properties(T, a, R)
    Cp_R = a(1)./T.^2 + a(2)./T + a(3) + a(4).*T + a(5).*T.^2 + a(6).*T.^3 + a(7).*T.^4;

    h_RT = -a(1)./T.^2 + a(2).*log(T)./T + a(3) + a(4).*T./2 + a(5).*T.^2./3 + a(6).*T.^3./4 + a(7).*T.^4./5 + a(8)./T;

    Cp = Cp_R .* R;
    h = h_RT .* R .* T;
end

function G = grad(A, h)
    d = ndims(A);
    G = zeros([size(A) 3]);
    idxOut = repmat({':'}, 1, d+1);
    for k = 1:3
        idxOut{d+1} = k;
        G(idxOut{:}) = centralDiff1(A, k, h(k));
    end
end

function div = div_onesided(A, s, h)
    % div = sum_k d(A_k)/dx_k. A의 4번째 차원이 미분 방향 k
    div = 0;
    for k = 1:3
        div = div + oneSidedDiff1(A(:,:,:,k,:), k, s, h(k));
    end
    div = permute(div, [1 2 3 5 4]); % 성분축 제거
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

function [rho, u, E] = apply_bc(rho, u, E, par)
    % 격자: 1번째/끝 인덱스가 고스트층. y가 커지는 방향이 고도(위쪽).
    % 바닥(y=1) : 고정벽, no-slip     -> 속도 반사(u=0 @ 벽면), T는 단열(zero-gradient)
    % 나머지 5면(x 양끝, y=끝(위), z 양끝) : 자유유출(zero-gradient, 안쪽 값 복제)
    % 단, y 양끝의 압력은 zero-gradient 대신 정역학 평형(dp/dy = rho*g_y)을 따른다.

    % --- 자유유출: x 양끝 ---
    rho([1 end],:,:) = rho([2 end-1],:,:); u([1 end],:,:,:) = u([2 end-1],:,:,:); E([1 end],:,:) = E([2 end-1],:,:);

    % --- 자유유출: z 양끝 ---
    rho(:,:,[1 end]) = rho(:,:,[2 end-1]); u(:,:,[1 end],:) = u(:,:,[2 end-1],:); E(:,:,[1 end]) = E(:,:,[2 end-1]);

    % --- y 양끝 ---
    % 고스트 셀은 안쪽 셀과 T, |u|가 같아 E/rho도 같으므로, 압력비만큼 rho와 E를 함께 바꾼다
    in = [2, size(rho,2)-1];
    e = E(:,in,:) ./ rho(:,in,:) - 0.5 .* sum(u(:,in,:,:).^2, 4);
    T = temperature(e, par);
    p_ratio = exp(par.g(2) .* [-par.h(2), par.h(2)] ./ (par.R .* T));
    rho(:,[1 end],:) = rho(:,in,:) .* p_ratio;
    E(:,[1 end],:) = E(:,in,:) .* p_ratio;
    u(:,end,:,:) = u(:,end-1,:,:); % 위쪽: 자유유출
    u(:,1,:,:) = -u(:,2,:,:);      % 바닥: 벽면에서 속도(법선+접선 모두)가 0이 되도록 반사
end
