function navier_stokes_solver()

    % 초깃값
    dx = 0.1; dy = 0.1; dz = 0.1; dt = 0.01;
    Nx = 21; Ny = 21; Nz = 21; Nt = 101;

    T_initial = 293.15; %K 초기온도 (20 °C)
    p = ones(Nx,Ny,Nz) * 101325; %Pa 초기압력
    T = ones(Nx,Ny,Nz) * T_initial; %K 초기온도
    properties = air_properties(T_initial, p(1));
    rho = ones(Nx,Ny,Nz) * properties.rho; %kg/m^3 초기밀도
    C_p = properties.Cp; %J/(kg*K) 비열
    mu = properties.mu; %Pa*s 점성계수
    k = properties.k; %W/(m*K) 열전도도
    Pr = C_p * mu / k; %프란틀 수

    fprintf(['Initial air properties at %.2f K, %.0f Pa: Cp=%.3f J/(kg*K), ', ...
        'rho=%.6f kg/m^3, mu=%.9g Pa*s, k=%.8g W/(m*K), Pr=%.5f\n'], ...
        T_initial, p(1), C_p, properties.rho, mu, k, Pr);
end

function properties = air_properties(T, P)
    % 건조 공기 물성. CoolProp은 200~2000 K에 사용한다.
    % 2000~6000 K에서는 NASA CEA Air 계수와 기존 근사식을 사용한다.
    if ~isscalar(T) || ~isscalar(P) || ~isfinite(T) || ~isfinite(P)
        error('navier_stokes_solver:InvalidState', ...
            'Temperature and pressure must be finite scalar values.');
    end

    if T < 200.0 || T > 6000.0
        error('navier_stokes_solver:TemperatureOutOfRange', ...
            'Temperature must be between 200 K and 6000 K.');
    end

    if T <= 2000.0
        properties.Cp = coolprop_value('Cpmass', T, P);
        properties.rho = coolprop_value('Dmass', T, P);
        properties.mu = coolprop_value('V', T, P);
        properties.k = coolprop_value('L', T, P);
        return;
    end

    % CoolProp Air의 상한(2000 K) 밖: NASA CEA Air Cp를 사용한다.
    % 이 고온 근사는 기존 솔버의 이상기체/Sutherland 가정을 유지한다.
    R_universal = 8314.46261815324; % J/(kmol*K)
    molecular_weight = 28.9651159; % kg/kmol
    R = R_universal / molecular_weight;
    T0 = 273.15; % K
    S = 110.4; % K
    mu0 = 1.716e-5; % Pa*s
    properties.Cp = nasa_cea_air_cp(T);
    properties.rho = P / (R * T);
    properties.mu = mu0 * (T/T0)^(3/2) * (T0 + S)/(T + S);
    properties.k = properties.Cp * properties.mu / 0.71;
end

function cp = nasa_cea_air_cp(T)
    % NASA CEA thermo.inp의 고정 조성 건조 공기(Air) Cp 계수: 300~6000 K.
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
    cp = (R_universal / molecular_weight) * ...
        sum(coefficients .* T.^exponents);
end

function value = coolprop_value(property_name, T, P)
    % CoolProp 8의 nanobind 함수를 MATLAB에서 직접 호출하지 않고,
    % Python 내부에서 실행한 뒤 표준 float 값만 MATLAB으로 가져온다.
    coolprop_setup();
    code = sprintf([ ...
        'from CoolProp.CoolProp import PropsSI\n', ...
        'result = float(PropsSI(''%s'', ''T'', %.17g, ''P'', %.17g, ''Air''))'], ...
        property_name, T, P);
    try
        value = double(pyrun(code, "result"));
    catch exception
        error('navier_stokes_solver:CoolPropCalculationFailed', ...
            'CoolProp failed for %s at T=%g K, P=%g Pa: %s', ...
            property_name, T, P, exception.message);
    end
end

function coolprop_setup()
    persistent initialized
    if ~isempty(initialized) && initialized
        return;
    end

    source_dir = fileparts(mfilename('fullpath'));
    
    % 운영체제에 맞춰 가상환경 파이썬 실행 파일 경로 설정
    if ispc
        python_executable = fullfile(fileparts(source_dir), '.venv', 'Scripts', 'python.exe');
    else
        python_executable = fullfile(fileparts(source_dir), '.venv', 'bin', 'python');
    end

    if ~isfile(python_executable)
        error('navier_stokes_solver:PythonMissing', ...
            'The project Python environment is missing: %s', python_executable);
    end

    environment = pyenv;
    if string(environment.Status) == "NotLoaded"
        pyenv('Version', python_executable);
    end

    try
        pyrun("from CoolProp.CoolProp import PropsSI");
    catch exception
        error('navier_stokes_solver:CoolPropImportFailed', ...
            'Could not import CoolProp from the project venv: %s', exception.message);
    end
    initialized = true;
end