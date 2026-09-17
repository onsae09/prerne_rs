# 개발 환경 설정

`navier_stokes_solver.m`은 200~2000 K의 건조 공기 물성을 CoolProp에서 읽고,
그보다 높은 온도에서는 저장소에 포함된 NASA CEA 계수를 사용한다.

아래 가이드는 macOS, Python 3.11, CoolProp 8.0.0 기준이다.

## 1. 필요한 폴더 구조

솔버는 저장소의 부모 폴더에 있는 `.venv`와 저장소 내부의 `.vendor`를 사용한다.

```text
prerne/
├── .venv/
└── prerne_rs/
    ├── .vendor/
    └── navier_stokes_solver.m
```

이 저장소를 `prerne_rs`라는 이름으로 clone한 다음, 저장소의 부모 폴더에서
아래 명령을 실행한다.

## 2. Python 가상환경 생성

```bash
cd /path/to/prerne
python3.11 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
```

Python 버전을 확인한다.

```bash
.venv/bin/python --version
```

예상 결과:

```text
Python 3.11.x
```

## 3. CoolProp 설치

CoolProp은 저장소 내부의 `.vendor`에 프로젝트 전용으로 설치한다.

```bash
.venv/bin/python -m pip install --target prerne_rs/.vendor CoolProp==8.0.0
```

설치를 확인한다.

```bash
.venv/bin/python -c "import sys; sys.path.insert(0, 'prerne_rs/.vendor'); from CoolProp.CoolProp import PropsSI; print(PropsSI('Cpmass', 'T', 293.15, 'P', 101325, 'Air'))"
```

약 `1006.144`가 출력되면 정상이다.

## 4. MATLAB에서 실행

MATLAB을 새로 실행하고 명령 창에서 다음을 입력한다. `/path/to/prerne`은 실제
프로젝트 경로로 바꾼다.

```matlab
cd('/path/to/prerne/prerne_rs')
pyenv('Version', '/path/to/prerne/.venv/bin/python')
clear functions
rehash
navier_stokes_solver
```

정상 실행 예시:

```text
Initial air properties at 293.15 K, 101325 Pa: Cp=1006.144 J/(kg*K), rho=1.204575 kg/m^3, mu=1.82056752e-05 Pa*s, k=0.025873828 W/(m*K), Pr=0.70796
```

`navier_stokes_solver.m`도 Python이 아직 로드되지 않은 경우 부모 폴더의
`.venv/bin/python`을 자동으로 선택한다. 위에서 `pyenv`를 명시하면 잘못된 Python이
먼저 로드되는 문제를 예방할 수 있다.

## 5. 문제 해결

### `No module named 'CoolProp'`

설치 위치와 파일 존재 여부를 확인한다.

```bash
ls prerne_rs/.vendor/CoolProp
```

없다면 3단계의 설치 명령을 다시 실행한다.

### MATLAB이 다른 Python을 이미 로드한 경우

현재 상태를 확인한다.

```matlab
pyenv
```

가장 확실한 해결 방법은 MATLAB을 재시작한 뒤, 다른 Python 코드를 실행하기 전에
4단계의 `pyenv` 명령부터 실행하는 것이다.

사용 중인 MATLAB 버전과 실행 모드에서 종료가 지원된다면 다음과 같이 다시
설정할 수도 있다.

```matlab
terminate(pyenv)
pyenv('Version', '/path/to/prerne/.venv/bin/python')
clear functions
navier_stokes_solver
```

### 코드를 수정했는데 이전 오류가 계속 나오는 경우

MATLAB 함수 캐시를 비운다.

```matlab
clear functions
rehash
navier_stokes_solver
```

## 물성 모델 범위

- 200~2000 K: CoolProp `Air` 모델로 `Cp`, 밀도, 점도, 열전도도 계산
- 2000~6000 K: NASA CEA `Air` 계수와 기존 고온 근사 사용
- 200 K 미만 또는 6000 K 초과: 범위 오류 발생

