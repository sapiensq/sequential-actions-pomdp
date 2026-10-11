# Sequential Action and Beliefs under Partially Observable DSGE Environments

[English](README.md) | **한국어**

다음 논문의 MATLAB 코드입니다.

> Kim, Seong-Hoon. 2012. "Sequential Action and Beliefs Under Partially Observable DSGE Environments." *Computational Economics* 40 (3): 219–244. [doi:10.1007/s10614-012-9323-1](https://doi.org/10.1007/s10614-012-9323-1)

## 논문 소개

이 논문은 인공지능 분야에서 빌려 온 두 가지 질문으로 DSGE 모형을 분류합니다. 경제주체가 경제의 상태를 관측할 수 있는가, 그리고 자신의 행동이 상태의 변화에 영향을 주는가? 그러면 세 부류가 나옵니다.

| 부류 | 상태 관측? | 행동이 상태에 영향? | 예 |
|---|---|---|---|
| MDP (마르코프 결정과정) | 예 | 예 | 표준 RBC 모형 |
| HMM (은닉 마르코프 모형) | 아니오 | 아니오 | 외생 충격에 대한 신호추출 |
| **POMDP** (부분관측 MDP) | 아니오 | 예 | 이 논문 |

선형 합리적 기대 모형에서 POMDP의 해는 세 부분으로 이루어집니다. 믿음(belief)에 기초한 행동 규칙, 실제 상태의 전이, 그리고 칼만 필터 형태의 믿음 갱신입니다. 논문은 믿음 갱신 이득(gain)을 "관측한 것"을 "믿는 것"으로 보내는 연산자의 고정점으로 구합니다. MDP와 HMM은 특수한 경우로 포함됩니다.

예제는 신고전파 성장모형으로, 경제주체는 산출은 관측하지만 자본과 기술의 영구·일시 성분은 관측하지 못합니다. 1% "인식 충격"(자본을 과대평가)이 발생하면 소비가 급증한 뒤 투자가 줄고, 회복에 100분기 넘게 걸립니다.

## 파일

| 파일 | 역할 |
|---|---|
| `compecon_kim_run.m` | **실행 진입점.** 모형을 풀고, POMDP와 HMM의 믿음 고정점을 반복 계산하고, 인식 충격을 시뮬레이션해 그림을 그립니다. |
| `compecon_kim_model.m` | 균형 조건을 기호(symbolic)로 정의. |
| `compecon_kim_ss.m` | 캘리브레이션(α = 0.36, β = 0.99, δ = 0.025, ρ = 0.99, N = 1/3)과 정상상태. |
| `anal_deriv.m`, `num_eval.m`, `gx_hx.m`, `gxx_hxx.m`, `gss_hss.m` | Schmitt-Grohé와 Uribe의 섭동(perturbation) 해법 (아래 참조). |

## 실행

MATLAB과 **Symbolic Math Toolbox**가 필요합니다.

```matlab
compecon_kim_run
```

각 믿음 고정점의 반복 횟수를 출력하고(논문: POMDP 221회, HMM 12회) 그림 세 개를 그립니다. 첫 번째 그림이 논문 그림 1·2의 6패널 충격반응입니다. 그림 1(40분기)에 맞추려면 `irdate = 41`, 그림 2(200분기)에 맞추려면 `irdate = 201`로 바꾸세요.

표 2의 이득 행렬과 공분산(`Hcx`, `Mx`, `Mz`, `covxhat`, `covxifr`)은 작업공간에 남습니다.

논문은 1차 근사해만 사용합니다. 2차 항(`gxx`, `hxx`, `gss`, `hss`)은 계산되지만 쓰이지 않습니다.

## 제3자 코드

`anal_deriv.m`, `num_eval.m`, `gx_hx.m`, `gxx_hxx.m`, `gss_hss.m`의 저작권은 Stephanie Schmitt-Grohé와 Martín Uribe에게 있으며, "Solving Dynamic General Equilibrium Models Using a Second-Order Approximation to the Policy Function" (*JEDC* 28, 2004, 755–775)의 공개 코드입니다. `gx_hx.m`은 Paul Klein의 `solab.m`을 바탕으로 했고 Ryan Chahrour가 수정했습니다. 원본은 <https://www.columbia.edu/~mu2166/2nd_order.htm>에서 받을 수 있습니다. 사용 시 해당 논문을 인용해 주세요.
