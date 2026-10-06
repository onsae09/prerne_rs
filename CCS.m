clear; close all;

% 유한차분
% (b,r,xi,eta)
% b ∈ [1:6], -1 <= xi <= 1, -1 <= eta <= 1
db = 1; dr = 1e-3; dxi = 2e-2; deta = 2e-2; dt = 1e-2;
Nb = 6; Nr = 1e2; Nxi = 2/dxi; Neta = 2/deta; Nt = 1e1;
b = [1:db:db*Nb]; r = [2e-2:dr:dr*Nr]; xi = [-1:dxi:1]; eta = [-1:deta:1]; t = dt:dt:dt*Nt;
%[B, R, XI, ETA] = ndgrid(b, r, xi, eta);

% 상수
nu = 1.56e-5; % m^2/s
g = [0, -9.81, 0]; % m/s^2

%{
ex = ones(size(XI),1);
ey = ones(size(ETA),1);
ez = ones(size(ETA),1);

Tx = spdiags([-ex, 2*ex, -ex], -1:1, size(XI,1), size(XI,1))/dx^2;
Ty = spdiags([-ey, 2*ey, -ey], -1:1, size(ETA,1), size(ETA,1))/dy^2;
Tz = spdiags([-ez, 2*ez, -ez], -1:1, size(ETA,1), size(ETA,1))/dz^2;

% x 양면: phi = 0
Tx(1,1)     = 3/dx^2;
Tx(end,end) = 3/dx^2;

% y 아래: dphi/dn = 0
Ty(1,1) = 1/dy^2;

% y 위: phi = 0
Ty(end,end) = 3/dy^2;

% z 양면: phi = 0
Tz(1,1)     = 3/dz^2;
Tz(end,end) = 3/dz^2;

Ix = speye(size(XI,1));
Iy = speye(size(ETA,1));
Iz = speye(size(ETA,1));

A = kron(Iz, kron(Iy, Tx)) ...
  + kron(Iz, kron(Ty, Ix)) ...
  + kron(Tz, kron(Iy, Ix));

tic
A = decomposition(A, 'chol')
s = whos("A").bytes/2^30
elapsed = toc

tic
A = chol(A);
elapsed = toc
s = whos("A").bytes/2^30
tic
save('cholesky.mat', 'A', '-v7.3')
elapsed = toc

%}

R = 287.107; % J/(kg*K)

% 초기 설정값
u.r = zeros(Nb, Nr+3, Nxi, Neta);
u.xi = zeros(Nb, Nr+3, Nxi, Neta);
u.eta = zeros(Nb, Nr+3, Nxi, Neta);
T = 293.15; % K

% 변수 초기화
w = ones(Nb, Nr+2, Nxi, Neta) .* R .* T; % Pa*m^3/kg

U = zeros([Nt,size(u)]);
W = zeros([Nt,size(w)]);
for i = 1:1
    [u, w] = FDM(u, w, nu, g, dx, dy, dz, dt, A);
    W(i,:,:,:) = w;
end

function [u2, w2] = FDM(u, w, nu, g, dx, dy, dz, dt, A)
    tic
    duu = (-uu(u, dx, dy, dz) - G(w, dx, dy, dz) + nu .* L(u, dx, dy, dz) + g) * dt;
    du = zeros(size(u),"like",u);
    du(2:end-1,2:end-1,2:end-1,:) = duu;
    u2 = u + du;
    [u2, ~] = bc(u2, w);
    f = D(u2, dx, dy, dz)/dt;
    p = reshape(A \ -f(:), size(f));
    phi = zeros(size(w),"like",w);
    phi(2:end-1,2:end-1,2:end-1) = p;
    phi = phi_bc(phi);
    w2 = w + phi - nu * D(u2, dx, dy, dz);
    u2(2:end-1,2:end-1,2:end-1,:) = u2(2:end-1,2:end-1,2:end-1,:) - dt * G(phi, dx, dy, dz);
    [u2, w2] = bc(u2, w2);
    elapsed = toc
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

    out = zeros(size(u), 'like', u);

    for i = 1:3
        for j = 1:3
            out = out + J(:,:,:,i,j) .* J(:,:,:,j,i);
        end
    end
end

function [u, w] = bc(u, w)
    % 격자: 1번째/끝 인덱스가 고스트층. y가 커지는 방향이 고도(위쪽).
    % 바닥(y=1) : 고정벽, no-slip     -> 속도 반사(u=0 @ 벽면), w는 단열(zero-gradient)
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
function phi = phi_bc(phi)

    % x 양면: 실제 경계면에서 phi=0
    phi(1,:,:)   = -phi(2,:,:);
    phi(end,:,:) = -phi(end-1,:,:);

    % z 양면: 실제 경계면에서 phi=0
    phi(:,:,1)   = -phi(:,:,2);
    phi(:,:,end) = -phi(:,:,end-1);

    % y 위: 실제 경계면에서 phi=0
    phi(:,end,:) = -phi(:,end-1,:);

    % y 아래 고체벽: dphi/dn=0
    phi(:,1,:) = phi(:,2,:);
end