% 초깃값
dx = 0.1; dy = 0.1; dz = 0.1; dt = 0.01;
Nx = 20; Ny = 20; Nz = 20; Nt = 100;

T_initial = 303.0; % K 초기온도 (30 °C)
p = ones(Nx,Ny,Nz) * 101325; % Pa 초기압력 (3D 텐서)
T = ones(Nx,Ny,Nz) * T_initial; % K 초기온도 (3D 텐서)

% 3D 텐서를 직접 입력
properties = air_properties(T, p);

rho = properties.rho % kg/m^3 (3D 텐서)
C_p = properties.Cp;  % J/(kg*K) (3D 텐서)
mu  = properties.mu;  % Pa*s (3D 텐서)
k   = properties.k;   % W/(m*K) (3D 텐서)
Pr  = C_p .* mu ./ k; % 프란틀 수 (3D 텐서 요소별 연산)

function properties = air_properties(T, P)
    % 건조 공기 물성치 계산 (스칼라 및 3D 텐서 지원)
    if any(T(:) < 200.0) || any(T(:) > 2000.0)
        error('navier_stokes_solver:TemperatureOutOfRange', ...
            'Temperature must be between 200 K and 2000 K.');
    end

    % 입력이 배열/텐서인 경우: 1차원 변환 후 일괄 계산 및 Reshape
    grid_size = size(T);
    
    % MATLAB 1차원 배열 -> NumPy 배열 전달로 속도 최적화
    T_np = py.numpy.array(T(:)');
    P_np = py.numpy.array(P(:)');

    properties.Cp  = coolprop_array('Cpmass', T_np, P_np, grid_size);
    properties.rho = coolprop_array('Dmass',  T_np, P_np, grid_size);
    properties.mu  = coolprop_array('V',      T_np, P_np, grid_size);
    properties.k   = coolprop_array('L',      T_np, P_np, grid_size);
end

function value_3d = coolprop_array(property_name, T_np, P_np, grid_size)
    coolprop_setup();
    try
        % CoolProp C++ 엔진으로 일괄 연산 실행
        result_np = py.CoolProp.CoolProp.PropsSI(property_name, 'T', T_np, 'P', P_np, 'Air');
        
        % 파이썬 결과를 double로 변환 후 원래 3D 격자 모양으로 복원
        value_3d = reshape(double(result_np), grid_size);
    catch exception
        error('navier_stokes_solver:CoolPropCalculationFailed', ...
            'CoolProp failed for %s: %s', property_name, exception.message);
    end
end

function coolprop_setup()
    persistent initialized
    if ~isempty(initialized) && initialized
        return;
    end

    source_dir = fileparts(mfilename('fullpath'));
    
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
        pyenv('Version', python_executable, 'ExecutionMode', 'InProcess');
    end

    try
        py.importlib.import_module('CoolProp.CoolProp');
        py.importlib.import_module('numpy');
    catch exception
        error('navier_stokes_solver:CoolPropImportFailed', ...
            'Could not import CoolProp or NumPy from the project venv: %s', exception.message);
    end
    initialized = true;
end