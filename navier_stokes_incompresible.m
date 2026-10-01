clear; clc; close all;

tic

% 유한차분
dx = 0.1; dy = 0.1; dz = 0.1; dt = 1e-2;
Nx = 1e2; Ny = 1e2; Nz = 1e2; Nt = 1e1;
Lx = dx*Nx; Ly = dy*Ny; Lz = dz*Nz; Lt = dt*Nt;
x_grid = -dx:dx:Lx; y_grid = -dy:dy:Ly; z_grid = -dz:dz:Lz; t = dt:dt:Lt;
[X, Y, Z] = ndgrid(x_grid, y_grid, z_grid);

% 상수
nu = 1.716e-5; 
T0 = 273.15; % K
S = 110.4; % K
rho = 1.225; % kg/m^3

%G = 6.67430e-11; % m^3/(kg*s^2)
%m = 5.972e24; % kg
%r0 = [-Lx/2 6.371e6 -Lz/2];
%rx = X + r0(1);
%ry = Y + r0(2);
%rz = Z + r0(3);
%r_norm = sqrt(rx.^2 + ry.^2 + rz.^2);
%g = -G * m ./ r_norm.^3 .* cat(4, rx, ry, rz);
g = reshape([0, -9.81, 0], [1,1,1,3]);

ex = ones(Nx,1);
ey = ones(Ny,1);
ez = ones(Nz,1);

Tx = spdiags([-ex, 2*ex, -ex], -1:1, Nx, Nx);
Ty = spdiags([-ey, 2*ey, -ey], -1:1, Ny, Ny);
Tz = spdiags([-ez, 2*ez, -ez], -1:1, Nz, Nz);

Ix = speye(Nx);
Iy = speye(Ny);
Iz = speye(Nz);

A = kron(Iz, kron(Iy, Tx)) ...
  + kron(Iz, kron(Ty, Ix)) ...
  + kron(Tz, kron(Iy, Ix));
elapsed = toc
tic
A = decomposition(A,"chol");
elapsed = toc
R = 287.107; % J/(kg*K)

% 초기 설정값
u = zeros(Nx+2, Ny+2, Nz+2, 3);
T = 303.15; % K

% 변수 초기화
T = ones(Nx+2, Ny+2, Nz+2) * T;
w = R .* T; % Pa*m^3/kg

U = zeros([Nt,size(u)]);
W = zeros([Nt,size(w)]);
for i = 1:Nt
    [u, w] = FDM(u, w, nu, g, dx, dy, dz, dt, A);
    U(:,:,:,:,i) = u;
end

function [u2, w2] = FDM(u, w, nu, g, dx, dy, dz, dt, A)
    tic
    duu = (-uu(u, dx, dy, dz) - G(w, dx, dy, dz) + nu .* L(u, dx, dy, dz) + g) * dt;
    du = zeros(size(u),"like",u);
    du(2:end-1,2:end-1,2:end-1,:) = duu;
    u2 = u + du;
    [u2, ~] = bc(u2, w);
    f = D(u2, dx, dy, dz)/dt;
    phi = reshape(A \ -f(:), size(f));
    p = zeros(size(w),"like",w);
    p(2:end-1,2:end-1,2:end-1) = phi;
    u2(2:end-1,2:end-1,2:end-1,:) = u2(2:end-1,2:end-1,2:end-1,:) - G(p, dx, dy, dz);
    w2 = w + p;
    [u2, w2] = bc(u2, w2);
    toc
end

function out = G(A, dx, dy, dz)
    h = [dx, dy, dz];
    n = [size(A,1)-2, size(A,2)-2, size(A,3)-2];
    nc = size(A,4);
    assert(nc == 1 || nc == 3);

    if nc == 1
        out = zeros([n, 3], 'like', A);
        for j = 1:3
            out(:,:,:,j) = fd1(A, h(j), j);
        end
    else
        out = zeros([n, 3, 3], 'like', A);
        for i = 1:3
            for j = 1:3
                out(:,:,:,i,j) = fd1(A(:,:,:,i), h(j), j);
            end
        end
    end
end

function out = D(A, dx, dy, dz)
    out = fd1(A(:,:,:,1), dx, 1) ...
        + fd1(A(:,:,:,2), dy, 2) ...
        + fd1(A(:,:,:,3), dz, 3);
end

function out = L(A, dx, dy, dz)
    h = [dx, dy, dz];
    n = [size(A,1)-2, size(A,2)-2, size(A,3)-2];
    nc = size(A,4);
    assert(nc == 1 || nc == 3);

    out = zeros([n, nc], 'like', A);
    for i = 1:nc
        Ai = A(:,:,:,i);
        for j = 1:3
            out(:,:,:,i) = out(:,:,:,i) + fd2(Ai, h(j), j);
        end
    end

    if nc == 1
        out = reshape(out, n);
    end
end

function out = fd1(A, h, dim)
    % 물리 셀의 1차 중앙차분
    c = {2:size(A,1)-1, 2:size(A,2)-1, 2:size(A,3)-1};
    p = c; m = c;
    p{dim} = p{dim}+1;
    m{dim} = m{dim}-1;
    out = (A(p{:}) - A(m{:})) / (2*h);
end

function out = fd2(A, h, dim)
    % 물리 셀의 2차 중앙차분
    c = {2:size(A,1)-1, 2:size(A,2)-1, 2:size(A,3)-1};
    p = c; m = c;
    p{dim} = p{dim}+1;
    m{dim} = m{dim}-1;
    out = (A(p{:}) - 2*A(c{:}) + A(m{:})) / h^2;
end

function out = uu(u, dx, dy, dz)
    % out = ∇·[(u·∇)u], 비압축성 유동
    J = G(u, dx, dy, dz);

    n = [size(u,1)-2, size(u,2)-2, size(u,3)-2];
    out = zeros(n, 'like', u);

    for i = 1:3
        for j = 1:3
            out = out + J(:,:,:,i,j) .* J(:,:,:,j,i);
        end
    end
end

function [u, w] = bc(u, w)
    % 격자: 1번째/끝 인덱스가 고스트층. y가 커지는 방향이 고도(위쪽).
    % 바닥(y=1) : 고정벽, no-slip     -> 속도 반사(u=0 @ 벽면), T/rho/w는 단열(zero-gradient)
    % 나머지 5면(x 양끝, y=끝(위), z 양끝) : 자유유출(zero-gradient, 안쪽 값 복제)
 
    % --- 자유유출: x 양끝 ---
    u(1,:,:,:)   = u(2,:,:,:);     w(1,:,:)   = w(2,:,:);
    u(end,:,:,:) = u(end-1,:,:,:); w(end,:,:) = w(end-1,:,:);
 
    % --- 자유유출: z 양끝 ---
    u(:,:,1,:)   = u(:,:,2,:);     w(:,:,1)   = w(:,:,2);
    u(:,:,end,:) = u(:,:,end-1,:); w(:,:,end) = w(:,:,end-1);
 
    % --- 자유유출: y 위쪽 끝 ---
    u(:,end,:,:) = u(:,end-1,:,:); w(:,end,:) = w(:,end-1,:);
 
    % --- 바닥(y=1) : no-slip 고정벽 ---
    u(:,1,:,:) = -u(:,2,:,:);      w(:,1,:)   = w(:,2,:);
end