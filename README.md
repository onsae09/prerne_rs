# Pre R&E

## prerne 코드

공기 유동(Navier–Stokes) 유한차분 솔버와 건조 공기 비열 계산 도구 모음.

## 구성

| 파일 | 내용 |
| --- | --- |
| `navier_stokes_solver.m` | 3차원 압축성 솔버. MacCormack 방식(전진/후진 차분 예측–수정), NASA 9계수로 `Cp`·엔탈피 계산, Sutherland 점성. 유효 온도 범위 200~1000 K. 초기 상태는 정역학 평형이며, 끝나면 `navier_stokes_result.mat`을 저장하고 z 중앙 단면을 그린다 |
| `navier_stokes_incompressible.m` | 3차원 비압축성 솔버(작업 중). 압력 포아송 행렬의 Cholesky 분해를 `cholesky.mat`으로 저장 |
| `solve_navier_stokes.m` | 2차원 비압축성 참고용 예제 (Re = 100) |
