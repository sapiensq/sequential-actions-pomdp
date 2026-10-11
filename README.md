# Sequential Action and Beliefs under Partially Observable DSGE Environments

**English** | [한국어](README.ko.md)

MATLAB code for

> Kim, Seong-Hoon. 2012. "Sequential Action and Beliefs Under Partially Observable DSGE Environments." *Computational Economics* 40 (3): 219–244. [doi:10.1007/s10614-012-9323-1](https://doi.org/10.1007/s10614-012-9323-1)

## What the paper does

The paper classifies DSGE models by two questions borrowed from the AI literature: can agents observe the state of the economy, and do their own actions change how that state evolves? This yields three classes:

| Class | State observed? | State affected by actions? | Example |
|---|---|---|---|
| MDP (Markov decision process) | yes | yes | standard RBC model |
| HMM (hidden Markov model) | no | no | signal extraction about exogenous shocks |
| **POMDP** (partially observable MDP) | no | yes | this paper |

For linear rational-expectations models, the POMDP solution has three parts: an action rule based on beliefs, the true state transition, and a Kalman-filter-type belief update. The paper finds the belief-updating gain as a fixed point of an operator that maps what agents observe into what they believe. MDP and HMM fall out as special cases.

The example is a neoclassical growth model in which agents observe output but not capital or the permanent and transitory technology components. A 1% "perception shock" (agents overvalue capital) triggers a consumption spree, then disinvestment, and a recovery that takes more than 100 quarters.

## Files

| File | Role |
|---|---|
| `compecon_kim_run.m` | **Entry point.** Solves the model, iterates the POMDP and HMM belief fixed points, simulates the perception shock, and plots the results. |
| `compecon_kim_model.m` | Equilibrium conditions written symbolically. |
| `compecon_kim_ss.m` | Calibration (α = 0.36, β = 0.99, δ = 0.025, ρ = 0.99, N = 1/3) and steady state. |
| `anal_deriv.m`, `num_eval.m`, `gx_hx.m`, `gxx_hxx.m`, `gss_hss.m` | Perturbation solver by Schmitt-Grohé and Uribe (see below). |

## Running

Requires MATLAB with the **Symbolic Math Toolbox**.

```matlab
compecon_kim_run
```

The script prints the number of iterations for each belief fixed point (the paper reports 221 for POMDP and 12 for HMM) and draws three figures. The first is the six-panel impulse response of the paper's Figs. 1 and 2. Set `irdate = 41` to match Fig. 1 (40 quarters) or `irdate = 201` to match Fig. 2 (200 quarters).

The gain matrices and covariances of Table 2 (`Hcx`, `Mx`, `Mz`, `covxhat`, `covxifr`) are left in the workspace.

The paper uses the first-order solution only. The second-order terms (`gxx`, `hxx`, `gss`, `hss`) are computed but not used.

## Third-party code

`anal_deriv.m`, `num_eval.m`, `gx_hx.m`, `gxx_hxx.m`, and `gss_hss.m` are © Stephanie Schmitt-Grohé and Martín Uribe, from the code accompanying "Solving Dynamic General Equilibrium Models Using a Second-Order Approximation to the Policy Function" (*JEDC* 28, 2004, 755–775). `gx_hx.m` is itself adapted from Paul Klein's `solab.m` and was modified by Ryan Chahrour. The original code is available at <https://www.columbia.edu/~mu2166/2nd_order.htm>. Please cite their paper if you use it.
