function navier_stokes_solver()

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
    T_initial = 300; %K 초기온도
    p = ones(Nx,Ny,Nz) * 101325; %Pa 초기압력
    T = ones(Nx,Ny,Nz) * T_initial; %K 초기온도
    rho = p./(R*T); %kg/m^3 초기밀도
    C_p = C_p(T); %J/(kg*K) 비열

    % Sutherland 식으로 초기온도에서의 점성계수를 구한다.
    mu = mu0 * (T_initial/T0)^(3/2) * (T0 + S)/(T_initial + S);
end

function C_p = C_p(T)
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
