function result = navier_stokes_solver(T_initial)
    arguments
        T_initial (1,1) double {mustBeFinite, mustBePositive} = 300.0
    end

    % NASA CEA dry-air polynomial을 MATLAB 안에서 직접 호출한다.
    C_p = cp_air_nasa(T_initial); % J/(kg*K)

    % 초깃값
    dx = 0.1; dy = 0.1; dz = 0.1; dt = 0.01;
    Nx = 21; Ny = 21; Nz = 21; Nt = 101;

    % 상수
    R_u = 8.31446261815324; %J/(mol*K) 일반기체상수
    M = 0.0289647; %kg/mol 공기분자량
    R = R_u/M; %J/(kg*K) 기체상수
    mu0 = 1.716e-5; %Pa*s 점성계수
    T0 = 273.15; %K 기준온도
    S = 110.4; %K Sutherland 상수
    Pr = 0.71; %프란틀 수
    % 이상기체식으로 초기 밀도를 계산한다.
    p = ones(Nx,Ny,Nz) * 101325; %Pa 초기압력
    T = ones(Nx,Ny,Nz) * T_initial; %K 초기온도
    rho = p./(R*T); %kg/m^3 초기밀도

    % Sutherland 식으로 초기온도에서의 점성계수를 구한다.
    mu = mu0 * (T_initial/T0)^(3/2) * (T0 + S)/(T_initial + S);

    % Python에서 전달된 C_p가 실제로 사용되는 부분이다.
    k = mu * C_p / Pr;          % W/(m*K), 열전도율
    alpha = k ./ (rho .* C_p);  % m^2/s, 열확산계수

    % 에너지 방정식의 열전도 항을 계산한다. 예제로 x=0 면을 초기
    % 온도보다 50 K 높은 벽으로 두고 내부 온도장을 시간 적분한다.
    T_hot = T_initial + 50.0;
    T(1,:,:) = T_hot;

    for n = 2:Nt
        T_old = T;
        laplacian_T = ...
            (T_old(3:Nx,2:Ny-1,2:Nz-1) - 2*T_old(2:Nx-1,2:Ny-1,2:Nz-1) + T_old(1:Nx-2,2:Ny-1,2:Nz-1))/dx^2 + ...
            (T_old(2:Nx-1,3:Ny,2:Nz-1) - 2*T_old(2:Nx-1,2:Ny-1,2:Nz-1) + T_old(2:Nx-1,1:Ny-2,2:Nz-1))/dy^2 + ...
            (T_old(2:Nx-1,2:Ny-1,3:Nz) - 2*T_old(2:Nx-1,2:Ny-1,2:Nz-1) + T_old(2:Nx-1,2:Ny-1,1:Nz-2))/dz^2;

        % rho*C_p*dT/dt = k*nabla^2(T)
        T(2:Nx-1,2:Ny-1,2:Nz-1) = T_old(2:Nx-1,2:Ny-1,2:Nz-1) + ...
            dt * (k ./ (rho(2:Nx-1,2:Ny-1,2:Nz-1) .* C_p)) .* laplacian_T;

        % 단열(영 법선 온도구배) 경계와 고온 벽 경계조건
        T(Nx,:,:) = T(Nx-1,:,:);
        T(:,1,:) = T(:,2,:);
        T(:,Ny,:) = T(:,Ny-1,:);
        T(:,:,1) = T(:,:,2);
        T(:,:,Nz) = T(:,:,Nz-1);
        T(1,:,:) = T_hot;
    end

    result = struct( ...
        'C_p', C_p, ...
        'mu', mu, ...
        'k', k, ...
        'alpha', mean(alpha, 'all'), ...
        'mean_temperature', mean(T, 'all'), ...
        'temperature', T);
end

function C_p = cp_air_nasa(T)
    % NASA CEA thermo.inp의 고정 조성 건조 공기(Air) 계수.
    % 유효 온도 범위는 300~6000 K이다.
    if T < 300.0 || T > 6000.0
        error('navier_stokes_solver:TemperatureOutOfRange', ...
            'Temperature must be between 300 K and 6000 K.');
    end

    if T <= 1000.0
        coefficients = [ ...
            1.009950160e4, -1.968275610e2, 5.009155110, ...
            -5.761013730e-3, 1.066859930e-5, ...
            -7.940297970e-9, 2.185231910e-12];
    else
        coefficients = [ ...
            2.415214430e5, -1.257874600e3, 5.144558670, ...
            -2.138541790e-4, 7.065227840e-8, ...
            -1.071483490e-11, 6.577800150e-16];
    end

    exponents = [-2, -1, 0, 1, 2, 3, 4];
    molecular_weight = 28.9651159; % kg/kmol
    R_universal = 8314.46261815324; % J/(kmol*K)
    C_p = (R_universal / molecular_weight) * ...
        sum(coefficients .* T.^exponents);
end
