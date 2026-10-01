# 💨🐥

## 연기세 코드

공기 유동(Navier–Stokes) 유한차분 솔버와 건조 공기 비열 계산 도구 모음.

## 구성

| 파일 | 내용 |
| --- | --- |
| `navier_stokes_solver.m` | 3차원 압축성 솔버. MacCormack 방식(전진/후진 차분 예측–수정), NASA 9계수로 `Cp`·엔탈피 계산, Sutherland 점성. 유효 온도 범위 200~400 K |
| `navier_stokes_incompressible.m` | 3차원 비압축성 솔버(작업 중). 압력 포아송 행렬의 Cholesky 분해를 `cholesky.mat`으로 저장 |
| `solve_navier_stokes.m` | 2차원 비압축성 참고용 예제 (Re = 100) |
| `cpcal/` | NASA CEA `thermo.inp`에서 건조 공기 `Cp`를 계산하는 Python 패키지 |
| `cpcal_usage.py` | `cpcal` 사용 예시 |

## 실행

MATLAB에서 저장소 폴더로 이동한 뒤 스크립트 이름을 입력한다.

```matlab
navier_stokes_solver
```

`cpcal`은 표준 라이브러리만 사용한다. 온도 단위는 K.

```bash
python cpcal/main.py 300
```

```python
from cpcal import cp_air
cp_air(300.0)  # J/(kg*K)
```

`cholesky.mat` 같은 `*.mat` 생성물은 git에서 추적하지 않는다.
