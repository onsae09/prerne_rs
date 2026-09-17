function solve_navier_stokes()
    % 1. 격자 및 영역 설정
    Lx = 1.0; Ly = 1.0;          % 영역 크기
    Nx = 41;  Ny = 41;           % 격자점 수
    dx = Lx / (Nx - 1);
    dy = Ly / (Ny - 1);
    x = linspace(0, Lx, Nx);
    y = linspace(0, Ly, Ny);
    [X, Y] = meshgrid(x, y);

    % 2. 물성치 및 계산 조건
    Re = 100;                    % 레이놀즈 수
    nu = 1 / Re;                 % 동점성 계수
    rho = 1.0;                   % 밀도
    dt = 0.001;                  % 시간 간격
    nt = 500;                    % 시간 반복 횟수
    nit = 50;                    % 압력 포아송 방정식 수렴 반복 수

    % 3. 변수 초기화
    u = zeros(Ny, Nx);           % x방향 유속
    v = zeros(Ny, Nx);           % y방향 유속
    p = zeros(Ny, Nx);           % 압력
    b = zeros(Ny, Nx);           % 압력 방정식 우변항

    % 4. 시간 루프 (Time-stepping Loop)
    for n = 1:nt
        un = u;
        vn = v;

        % [Step 1] 압력 포아송 방정식의 우변(b) 계산
        for i = 2:Ny-1
            for j = 2:Nx-1
                du_dx = (un(i, j+1) - un(i, j-1)) / (2*dx);
                dv_dy = (vn(i+1, j) - vn(i-1, j)) / (2*dy);
                du_dy = (un(i+1, j) - un(i-1, j)) / (2*dy);
                dv_dx = (vn(i, j+1) - vn(i, j-1)) / (2*dx);

                b(i,j) = rho * (1/dt * (du_dx + dv_dy) - ...
                    (du_dx^2 + 2*du_dy*dv_dx + dv_dy^2));
            end
        end

        % [Step 2] 압력 포아송 방정식 풀이 (자코비 반복법)
        for it = 1:nit
            pn = p;
            p(2:Ny-1, 2:Nx-1) = ((pn(2:Ny-1, 3:Nx) + pn(2:Ny-1, 1:Nx-2)) * dy^2 + ...
                                 (pn(3:Ny, 2:Nx-1) + pn(1:Ny-2, 2:Nx-1)) * dx^2 - ...
                                 b(2:Ny-1, 2:Nx-1) * dx^2 * dy^2) / (2 * (dx^2 + dy^2));

            % 압력 경계조건 (Neumann Boundary Condition)
            p(:, Nx) = p(:, Nx-1); % 우측
            p(:, 1)  = p(:, 2);    % 좌측
            p(1, :)  = p(2, :);    % 하단
            p(Ny, :) = 0;          % 상단 (기준 압력)
        end

        % [Step 3] 운동량 방정식을 통한 속도장 업데이트
        for i = 2:Ny-1
            for j = 2:Nx-1
                % 대류항 (Advection)
                u_adv = un(i,j)*(un(i,j) - un(i,j-1))/dx + vn(i,j)*(un(i,j) - un(i-1,j))/dy;
                v_adv = un(i,j)*(vn(i,j) - vn(i,j-1))/dx + vn(i,j)*(vn(i,j) - vn(i-1,j))/dy;

                % 확산항 (Viscous Diffusion)
                u_diff = nu * ((un(i,j+1) - 2*un(i,j) + un(i,j-1))/dx^2 + (un(i+1,j) - 2*un(i,j) + un(i-1,j))/dy^2);
                v_diff = nu * ((vn(i,j+1) - 2*vn(i,j) + vn(i,j-1))/dx^2 + (vn(i+1,j) - 2*vn(i,j) + vn(i-1,j))/dy^2);

                % 압력 구배항 (Pressure Gradient)
                dp_dx = (p(i, j+1) - p(i, j-1)) / (2*dx);
                dp_dy = (p(i+1, j) - p(i-1, j)) / (2*dy);

                % 속도 갱신
                u(i,j) = un(i,j) + dt * (-u_adv - dp_dx/rho + u_diff);
                v(i,j) = vn(i,j) + dt * (-v_adv - dp_dy/rho + v_diff);
            end
        end

        % [Step 4] 속도 경계조건 적용 (No-slip Condition)
        u(1, :)  = 0;       % 하단 벽면
        u(:, 1)  = 0;       % 좌측 벽면
        u(:, Nx) = 0;       % 우측 벽면
        u(Ny, :) = 1.0;     % 상단 뚜껑 (x방향 속도 u = 1.0으로 이동)

        v(1, :)  = 0;
        v(:, 1)  = 0;
        v(:, Nx) = 0;
        v(Ny, :) = 0;
    end

    % 5. 시각화 (Visualization)
    visualize_flow(X, Y, u, v, p);
end

function visualize_flow(X, Y, u, v, p)
    vel_mag = sqrt(u.^2 + v.^2); % 속도 크기

    figure('Color', 'w', 'Position', [100, 100, 1000, 400]);

    % 1) 속도 크기 분포 및 유선(Vector Field)
    subplot(1, 2, 1);
    contourf(X, Y, vel_mag, 25, 'LineColor', 'none');
    hold on;
    quiver(X(1:2:end, 1:2:end), Y(1:2:end, 1:2:end), ...
           u(1:2:end, 1:2:end), v(1:2:end, 1:2:end), 1.2, 'k');
    colormap(gca, 'jet');
    colorbar;
    title('속도 크기 (Velocity Magnitude) & 벡터장');
    xlabel('x'); ylabel('y');
    axis equal tight;
    hold off;

    % 2) 압력장 분포 (Pressure Field)
    subplot(1, 2, 2);
    contourf(X, Y, p, 25, 'LineColor', 'none');
    colormap(gca, 'parula');
    colorbar;
    title('압력 분포 (Pressure Field)');
    xlabel('x'); ylabel('y');
    axis equal tight;
end