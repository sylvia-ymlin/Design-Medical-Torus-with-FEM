# Design a Medical Torus

**Authors:** Yangmei Lin, Haote Liu

**Date:** January 4, 2024

**Applied FEM report (PDF):** [AppliedFEM.pdf](AppliedFEM.pdf)

---

## Introduction

This project is about designing a medical torus that releases hormones through diffusion. The objective is to identify appropriate parameters to shape the torus, control diffusion velocity, and ensure the delivery of a specified dose within a predetermined time frame.

The project consists of three parts. Part A is a simpler 1D model that introduces adaptive mesh refinement based on a posteriori error estimation. Part B involves implementing a piece-wise linear finite element approximation in two dimensions and studying the convergence analysis of a finite element method. Finally, in Part C, we will solve the 3D problem using FEniCS.

Through this project, we can deepen our understanding of utilizing finite element methods for solving physics and numerical models.

---

## Part A

In part A, we consider a simpler 1D model and perform adaptive mesh refinement based on a posteriori error estimation.

### Problem 1

Set $e = u - u_h$, then we have:

$$
||e'||^2_{L^2(I)} = \int_I e'^2 \, dx
$$

With Galerkin orthogonality:

$$
\int_I (u - u_h)' v' = 0, \quad \forall v \in \mathbf{V}_{h,0} \\
\int_I e' \pi e' dx = 0, \quad \pi e \in V_{h,0}
$$

Then we can obtain:

$$
\begin{aligned}
||e'||^2_{L^2(I)} &= \int_I e'^2 - e' \pi e' \, dx \\
&= \sum^n_{i=1} \int^{x_i}_{x_{i-1}} e'(e - \pi e)' \, dx \\
&= \text{[IBP]} \\
&= \sum^n_{i=1} \left( \int^{x_i}_{x_{i-1}} -e''(e - \pi e) \, dx + [e'(e - \pi e)]_{x_{i-1}}^{x_i} \right)
\end{aligned}
$$

Since $e$ and $\pi e$ coincide at the nodes:

$$
||e'||^2_{L^2(I)} = \sum^n_{i=1} \int^{x_i}_{x_{i-1}} -e''(e - \pi e) \, dx \quad (1)
$$

For $-e''$ on $I_i$ ($\alpha$ is the diffusion coefficient, so $\alpha > 0$):

$$
-\alpha e'' = -\alpha(u - u_h)'' = f + \alpha u_h''
$$

Plug this equation into (1):

$$
\begin{aligned}
\alpha ||e'||^2_{L^2(I)} &= \sum^n_{i=1} \int^{x_i}_{x_{i-1}} (f + \alpha u''_h)(e - \pi e) dx \\
&\le \sum^n_{i=1} ||f + \alpha u''_h||_{L^2(I_i)} ||e - \pi e||_{L^2(I_i)} \quad \text{(continuous Cauchy-Schwarz inequality)} \\
&\le \sum^n_{i=1} ||f + \alpha u''_h||_{L^2(I_i)} C_0 h_i ||e'||_{L^2(I_i)} \\
&\le C_0 \left( h^2_i \sum^n_{i=1} ||f + \alpha u''_h||^2_{L^2(I_i)} \right)^{\frac{1}{2}} \left( \sum^n_{i=1} ||e'||^2_{L^2(I_i)} \right)^{\frac{1}{2}} \quad \text{(discrete Cauchy-Schwarz inequality)} \\
&= C_0 \left( h^2_i \sum^n_{i=1} ||f + \alpha u''_h||^2_{L^2(I_i)} \right)^{\frac{1}{2}} ||e'||_{L^2(I)}
\end{aligned}
$$

Dividing both sides by $\alpha ||e'||_{L^2(I)}$ concludes the proof after squaring:

$$
||e'||^2_{L^2(I)} \le \left( \frac{C_0}{\alpha} \right)^2 \sum^n_{i=1} h_i^2 ||f + \alpha u''_h||^2_{L^2(I_i)} \le C \sum^n_{i=1} \eta_i^2
$$

### Problem 2

<p align="center">
  <img src="figures/Problem%20A_2.png" width="70%">
</p>

The first sub-figure illustrates the approximation solution $u_h$ for this problem. Through the Residual $R(u_h)$, we observe that the numerical solution is accurate in most intervals but exhibits significant error near the segmentation points of $f(x)$ at $x \in \{-0.8, -0.2, 0.2, 0.8\}$. It is also reflected in the mesh-size distribution. Since the algorithm refines the elements with the biggest contribution to the error, most of the new points inserted concentrate on the intervals near the mentioned segmentation points. Ultimately, we achieved an error smaller than TOL when $N = 40$.

---

## Part B

In Part B, we need to implement the piece-wise linear finite element approximation in two spatial dimensions and address the time-dependent problem. Subsequently, we will analyze its convergence.

### Problem 1

We consider the following 2D stationary problem:

$$
\begin{aligned}
-\triangle u(\mathbf{x}) &= f(\mathbf{x}), \quad \mathbf{x} \in \mathcal{B} \\
u(\mathbf{x}) &= u_{\text{exact}}(\mathbf{x}), \quad \mathbf{x} \in \partial\mathcal{B}
\end{aligned}
$$

with $f(\mathbf{x}) = 8\pi^2 \sin(2\pi x_1) \sin(2\pi x_2)$ and $u_{\text{exact}}(\mathbf{x}) = \sin(2\pi x_1) \sin(2\pi x_2)$.

We first derive the variational formulation. We multiply $-\triangle u(\mathbf{x}) = f(\mathbf{x})$ by a test function $v$, and integrate using Green's formula:

$$
\begin{gathered}
\int_{\mathcal{B}} -\triangle u v \, d\mathbf{x} = \int_{\mathcal{B}} f v \, d\mathbf{x} \\
\int_{\mathcal{B}} \nabla u \cdot \nabla v \, d\mathbf{x} - \int_{\partial\mathcal{B}} \partial_n u v \, ds = \int_{\mathcal{B}} f v \, d\mathbf{x}
\end{gathered}
$$

Then we see that $v$ and $\triangle v$ should be square-integrable, and in addition that $v = 0$ on $\partial\mathcal{B}$. Introducing the spaces:

$$
\begin{aligned}
\mathbf{V} &= \{ v : ||v||_{L^2(\mathcal{B})} + ||\nabla v||_{L^2(\mathcal{B})} < \infty \} \\
\mathbf{V}_0 &= \{ v \in \mathbf{V} : v|_{\partial\mathcal{B}} = 0 \} \\
\mathbf{V}_g &= \{ v \in \mathbf{V} : v|_{\partial\mathcal{B}} = u_{\text{exact}} \}
\end{aligned}
$$

We obtain the variational formulation:

Find $u \in \mathbf{V}_g$ such that:

$$
\int_{\mathcal{B}} \nabla u \cdot \nabla v \, d\mathbf{x} = \int_{\mathcal{B}} f v \, d\mathbf{x}, \quad \forall v \in \mathbf{V}_0
$$

Next, we construct the finite element spaces of the continuous piecewise linear polynomial:

$$
\begin{aligned}
\mathbf{V}_h &:= \{ v : v \in C^0(\mathcal{B}), v|_{\mathcal{K}} \in \mathcal{P}^1(\mathcal{K}), \forall \mathcal{K} \in \tilde{\mathcal{L}}_h \}, \quad \mathbf{V}_h \subset \mathbf{V} \\
\mathbf{V}_{h,0} &:= \{ v \in \mathbf{V}_h : v|_{\partial\mathcal{B}} = 0 \} \\
\mathbf{V}_{h,g} &:= \{ v \in \mathbf{V}_h : v|_{\partial\mathcal{B}} = u_{\text{exact}} \}
\end{aligned}
$$

Thus we obtain (GFEM):

Find $u_h \in \mathbf{V}_{h,g}$ such that:

$$
\int_{\mathcal{B}} \nabla u_h \cdot \nabla v \, d\mathbf{x} = \int_{\mathcal{B}} f v \, d\mathbf{x}, \quad \forall v \in \mathbf{V}_{h,0}
$$

By implementing it in Matlab with $h_{\max} = \{ \frac{1}{2}, \frac{1}{4}, \frac{1}{8}, \frac{1}{16}, \frac{1}{32} \}$, we obtain the numerical convergence rates in the energy norm as $\{ 0.7429, 1.9856, 1.4346, 1.1710 \}$. Figure 2 plots $h_{\max}$ versus the energy norm of the error and $h_{\max}$ versus $h^\gamma_{\max}$. We can see that these two lines are almost parallel, indicating good convergence performance.

<p align="center">
  <img src="figures/ProblemB1_TheConvergenceRateAnalysis.png" width="70%">
</p>

The convergence rate $\gamma$ is defined by:

$$
\gamma_i = \frac{\log(e_{i+1}/e_i)}{\log(h_{i+1}/h_i)}
$$

where $e = ||u - u_h||_E = \left( \int_{\mathcal{B}} (\nabla u - \nabla u_h) \cdot (\nabla u - \nabla u_h) \, d\mathbf{x} \right)^{1/2}$.

The solutions obtained using the coarsest and the finest meshes are plotted in Figure 3.

<p align="center">
  <img src="figures/ProblemB1_CoarsestMeshes.png" width="48%">
  <img src="figures/ProblemB1_FinestMeshes.png" width="48%">
</p>
<p align="center">
  <em>Left: Coarsest Meshes, $h_{\max} = 1/2$. Right: Finest Meshes, $h_{\max} = 1/32$.</em>
</p>

---

### Problem 2

The problem is defined as following:

For $t > 0$:

$$
\begin{aligned}
\partial_t u(\mathbf{x}, t) - \alpha \triangle u(\mathbf{x}, t) &= f(\mathbf{x}), \quad (\mathbf{x}, t) \in \mathcal{B} \times (0, T] \\
u(\mathbf{x}, t) &= 0, \quad (\mathbf{x}, t) \in \partial\mathcal{B} \times (0, T] \\
u(\mathbf{x}, 0) &= \begin{cases} \rho, & \mathbf{x} \in \tau \\ 0, & \mathbf{x} \in \mathcal{B} \setminus \tau \end{cases}, \quad \mathbf{x} \in \mathcal{B}
\end{aligned}
$$

Let $\mathbf{V}_0 = \{ v(\mathbf{x}, t) : ||v(\cdot, t)||^2 + ||\nabla v(\cdot, t)||^2 < \infty, v(\cdot, t) = 0 \text{ on } \partial\mathcal{B}, \forall t \ge 0 \}$.

Let $\mathbf{V}_{h,0} = \{ v(\mathbf{x}, t) : v(\mathbf{x}, t) \in \mathcal{C}^0(\tilde{\mathcal{L}}_h), v(\mathbf{x}, t)|_{\mathcal{K}} \in \mathcal{P}^1(\mathcal{K}), v(\mathbf{x}, t) = 0 \text{ on } \partial\mathcal{B}, \forall t \ge 0 \}$.

**(WF)** Find $u(\mathbf{x}, t) \in \mathbf{V}_0$ such that:

$$
\int_{\mathcal{B}} \partial_t u v \, d\mathbf{x} - \alpha \int_{\mathcal{B}} \triangle u v \, d\mathbf{x} = \int_{\mathcal{B}} f v \, d\mathbf{x}, \quad \forall v \in \mathbf{V}_0
$$

$$
\begin{aligned}
\int_{\mathcal{B}} f v \, d\mathbf{x} &= \int_{\mathcal{B}} \partial_t u v \, d\mathbf{x} + \alpha \int_{\mathcal{B}} \nabla u \nabla v \, d\mathbf{x} - \alpha \int_{\partial\mathcal{B}} \partial_n u v \, ds \\
&= \int_{\mathcal{B}} \partial_t u v \, d\mathbf{x} + \alpha \int_{\mathcal{B}} \nabla u \nabla v \, d\mathbf{x}
\end{aligned}
$$

Since $v = 0$ on $\partial\mathcal{B}$, the boundary term vanished.

**(GFEM)** Find $u_h(\mathbf{x}, t) \in \mathbf{V}_{h,0}$, such that:

$$
\int_{\mathcal{B}} \partial_t u_h v \, d\mathbf{x} + \alpha \int_{\mathcal{B}} \nabla u_h \nabla v \, d\mathbf{x} = \int_{\mathcal{B}} f v \, d\mathbf{x}, \quad \forall v \in \mathbf{V}_{h,0}
$$

Take the usual tent functions $\{ \varphi_j(\mathbf{x}) \}$ as a basis of $\mathbf{V}_{h,0}$ with time dependent coefficients.

For $u_h \in \mathbf{V}_{h,0}$, $\exists \{ \xi_j(t) \}$ such that:

$$
u_h(\mathbf{x}, t) = \sum_{N_j \in \mathcal{N}_h \setminus \mathcal{N}_b} \xi_j(t) \varphi_j(\mathbf{x})
$$

Where $\mathcal{N}_h$ is the set of all nodes on $\tilde{\mathcal{L}}_h$, and $\mathcal{N}_b$ is the set of all boundary nodes on $\tilde{\mathcal{L}}_h$. Insert it into (GFEM):

$$
\sum_{N_j \in \mathcal{N}_h \setminus \mathcal{N}_b} \frac{d\xi_j(t)}{dt} \int_{\mathcal{B}} \varphi_j \varphi_i \, d\mathbf{x} + \alpha \sum_{N_j \in \mathcal{N}_h \setminus \mathcal{N}_b} \xi_j(t) \int_{\mathcal{B}} \nabla \varphi_j \nabla \varphi_i \, d\mathbf{x} = \int_{\mathcal{B}} f \varphi_i \, d\mathbf{x}, \quad \forall N_i \in \mathcal{N}_h \setminus \mathcal{N}_b
$$

We obtain the following linear system of ordinary differential equations:

$$
\begin{aligned}
M \frac{d\boldsymbol{\xi}(t)}{dt} + \alpha A \boldsymbol{\xi}(t) &= \mathbf{b} \\
\boldsymbol{\xi}(0) &= \begin{cases} \rho, & \mathbf{x} \in \tau \\ 0, & \mathbf{x} \in \mathcal{B} \setminus \tau \end{cases}
\end{aligned}
$$

where:

$$
\begin{aligned}
M_{i,j} &= \int_{\mathcal{B}} \varphi_j \varphi_i \, d\mathbf{x} \\
A_{i,j} &= \int_{\mathcal{B}} \nabla \varphi_j \nabla \varphi_i \, d\mathbf{x} \\
b_i &= \int_{\mathcal{B}} f(\mathbf{x}) \varphi_i \, d\mathbf{x}
\end{aligned}
$$

Then approximate $\boldsymbol{\xi}(t)$ on discrete points $0 = t_0 < t_1 < \cdots < t_n = T$ using Crank-Nicolson method:

$$
M \frac{\boldsymbol{\xi}^{(n+1)} - \boldsymbol{\xi}^{(n)}}{t_{n+1} - t_n} + \frac{\alpha}{2} A (\boldsymbol{\xi}^{(n+1)} + \boldsymbol{\xi}^{(n)}) = \mathbf{b}
$$

Rearranging the terms:

$$
\left( \frac{1}{t_{n+1}-t_n} M + \frac{\alpha}{2} A \right) \boldsymbol{\xi}^{(n+1)} = \left( \frac{1}{t_{n+1}-t_n}M - \frac{\alpha}{2}A \right) \boldsymbol{\xi}^{(n)} + \mathbf{b}
$$

### Problem 3

In this problem, we implement a numerical integration in 2D to compute the mass loss:

$$
M(t) = \int_{\mathcal{B}} (u_0(\mathbf{x}) - u(\mathbf{x}, t)) \, d\mathbf{x}
$$

For computing it, we use the Trapezoidal rule in each element $\mathbf{K}_i$:

$$
\int_{\mathbf{K}_i} (u_0(\mathbf{x}) - u(\mathbf{x}, t)) \approx \frac{|\mathbf{K}_i|}{3} \sum^3_{i=1} (u_0(N_i) - u(N_i, t))
$$

The control parameters are set as $f(\mathbf{x}) = 0, \rho = 10, R = 0.5, r = 0.3, T = 30$. So we have the following iteration equation:

$$
\left( \frac{1}{t_{n+1}-t_n} M + \frac{\alpha}{2} A \right) \boldsymbol{\xi}^{(n+1)} = \left( \frac{1}{t_{n+1}-t_n} M - \frac{\alpha}{2} A \right) \boldsymbol{\xi}^{(n)}
$$

The amount of emitted hormone versus time is plotted in Figure 4 with $h_{\max} = 1/5$ and $h_{\max} = 1/20$ respectively.

<p align="center">
  <img src="figures/ProblemB3_MassLoss.png" width="70%">
</p>

Figure 5 displays the solutions at the initial time and the final time with $h_{\max} = 1/20$.

<p align="center">
  <img src="figures/ProblemB3_InitialSolution.png" width="48%">
  <img src="figures/ProblemB3_Solution_T_30.png" width="48%">
</p>
<p align="center">
  <em>Left: Initial solution (T=0). Right: Final solution (T=30).</em>
</p>

---

### Problem 4

In this problem, we need to solve a heat equation with a reaction term:

$$
\partial_t u(\mathbf{x}, t) - \alpha \triangle u(\mathbf{x}, t) = \beta(1 - \gamma u)u, \quad (\mathbf{x}, t) \in \mathcal{B} \times (0, T]
$$

**(GFEM)**

$$
\int_{\mathcal{B}} \partial_t u_h v \, d\mathbf{x} + \alpha \int_{\mathcal{B}} \nabla u_h \nabla v \, d\mathbf{x} = \int_{\mathcal{B}} \beta(1 - \gamma u) u v \, d\mathbf{x}
$$

For the right-hand side, we write it as $S(u) = \beta(1 - \gamma u)u$. Since $S(u)$ is nonlinear, we take a linear interpolant of it to compute it using the solution from the previous time step:

$$
S(u) \approx \pi_h S = \sum^{\mathcal{N}}_{j=1} S_j \varphi_j = \sum^{\mathcal{N}}_{j=1} \beta(1 - \gamma u_j) u_j \varphi_j
$$

Plugging them into the (GFEM):

$$
\sum^{\mathcal{N}}_{j=1} \frac{d\xi_j(t)}{dt} \int_{\mathcal{B}} \varphi_j \varphi_i \, d\mathbf{x} + \alpha \sum^{\mathcal{N}}_{j=1} \xi_j(t) \int_{\mathcal{B}} \nabla \varphi_j \nabla \varphi_i \, d\mathbf{x} = \sum^{\mathcal{N}}_{j=1} \beta(1 - \gamma u_j) u_j \int_{\mathcal{B}} \varphi_j \varphi_i \, d\mathbf{x}
$$

Linear system form:

$$
M \frac{d\boldsymbol{\xi}(t)}{dt} + \alpha A \boldsymbol{\xi}(t) = M (\beta \boldsymbol{\xi}^{(n)} - \beta \gamma \boldsymbol{\xi}^{(n)} \circ \boldsymbol{\xi}^{(n)})
$$

where "$\circ$" represents Hadamard product. Using Crank-Nicolson for the left-hand side and the previous time step for the right-hand side:

$$
M \frac{\boldsymbol{\xi}^{(n+1)} - \boldsymbol{\xi}^{(n)}}{t_{n+1}-t_n} + \frac{\alpha}{2} A (\boldsymbol{\xi}^{(n+1)} + \boldsymbol{\xi}^{(n)}) = \beta M \boldsymbol{\xi}^{(n)} - \beta \gamma M \boldsymbol{\xi}^{(n)} \circ \boldsymbol{\xi}^{(n)}
$$

Iteration equation:

$$
\left( \frac{M}{t_{n+1}-t_n} + \frac{\alpha}{2} A \right) \boldsymbol{\xi}^{(n+1)} = \left( \frac{M}{t_{n+1}-t_n} - \frac{\alpha}{2} A \right) \boldsymbol{\xi}^{(n)} + \beta M \boldsymbol{\xi}^{(n)} - \beta \gamma M \boldsymbol{\xi}^{(n)} \circ \boldsymbol{\xi}^{(n)}
$$

With $\beta = 0.2, \gamma = 0.5$:

<p align="center">
  <img src="figures/ProblemB4_beta_0.2_gamma_0.5_MassLoss.png" width="70%">
</p>

<p align="center">
  <img src="figures/ProblemB4_beta_0.2_gamma_0.5_t_0.png" width="48%">
  <img src="figures/ProblemB4_beta_0.2_gamma_0.5_t_30.png" width="48%">
</p>
<p align="center">
  <em>Left: T=0. Right: T=30.</em>
</p>

With $\beta = 1, \gamma = 0.2$:

<p align="center">
  <img src="figures/ProblemB4_beta_1.0_gamma_0.2_MassLoss.png" width="70%">
</p>

<p align="center">
  <img src="figures/ProblemB4_beta_1.0_gamma_0.2_t_0.0.png" width="48%">
  <img src="figures/ProblemB4_beta_1.0_gamma_0.2_t_30.png" width="48%">
</p>
<p align="center">
  <em>Left: T=0. Right: T=30.</em>
</p>

To study how changes in $\beta$ and $\gamma$ affect the result:

<p align="center">
  <img src="figures/ProblemB4_changeInBeta.png" width="48%">
  <img src="figures/ProblemB4_changeInGamma.png" width="48%">
</p>
<p align="center">
  <em>Left: Change in beta. Right: Change in gamma.</em>
</p>

Mass loss speed is positively correlated with $\beta$ and $\gamma$. Total mass loss is positively correlated with $\gamma$ and negatively correlated with $\beta$.

---

## Part C

In this part, we consider the full 3D model and solve the problems using FEniCS.

### Problem 1

The problem is defined as:

$$
\begin{aligned}
\partial_t u(\mathbf{x}, t) - \alpha \triangle u(\mathbf{x}, t) &= f(\mathbf{x}), \quad (\mathbf{x}, t) \in \mathcal{B} \times (0, T] \\
u(\mathbf{x}, t) &= 0, \quad (\mathbf{x}, t) \in \partial\mathcal{B} \times (0, T] \\
u(\mathbf{x}, 0) &= \begin{cases} \rho, & \mathbf{x} \in \mathcal{T} \\ 0, & \mathbf{x} \in \mathcal{B} \setminus \mathcal{T} \end{cases}
\end{aligned}
$$

where $f(\mathbf{x}) = 0$ and $\mathcal{T}$ denotes the torus:

$$
\mathcal{T} := \begin{cases} (R - \sqrt{x_1^2 + x_2^2})^2 + x_3^2 \le r^2, & \mathbf{x} \in \mathbf{R}^3 \\ |R - \sqrt{x_1^2 + x_2^2}| \le r, & \mathbf{x} \in \mathbf{R}^2 \end{cases}
$$

Parameters: $\rho = 10, R = 0.5, r = 0.2, T = 20$. Time discretization using backward difference:

$$
\int_{\mathcal{B}} (uv + \Delta t \nabla u \cdot \nabla v) \, dx = \int_{\mathcal{B}} u^{(n)} v \, dx
$$

#### 2D Mesh Solutions

<p align="center">
  <img src="figures/2D_initial_line1.png" width="48%">
  <img src="figures/2D_initial_line2.png" width="48%">
</p>
<p align="center">
  <em>Solutions for 2D mesh at T=0 (Image on left, plot along diagonal on right).</em>
</p>

<p align="center">
  <img src="figures/2D_final_line1.png" width="48%">
  <img src="figures/2D_final_line2.png" width="48%">
</p>
<p align="center">
  <em>Solutions for 2D mesh at T=20 (Image on left, plot along diagonal on right).</em>
</p>

#### 3D Mesh Solutions

<p align="center">
  <img src="figures/3D_initial_line1.png" width="48%">
  <img src="figures/3D_initial_line2.png" width="48%">
</p>
<p align="center">
  <em>Solutions for 3D mesh at T=0 (Image on left, plot along diagonal on right).</em>
</p>

<p align="center">
  <img src="figures/3D_final_line1.png" width="48%">
  <img src="figures/3D_final_line2.png" width="48%">
</p>
<p align="center">
  <em>Solutions for 3D mesh at T=20 (Image on left, plot along diagonal on right).</em>
</p>

The hottest point at $T = 20$ is approximately $10\%$ of the initial hottest point.

### Problem 2

Mass loss for different $\rho$ values:

<p align="center">
  <img src="figures/Part3_Task2_mass_loss.png" width="70%">
</p>

Higher initial density results in a steeper gradient and more rapid diffusion.

### Problem 3

Optimal parameters for torus design (least squares problem):

$$
\min F(\rho, R, r) = \sum^2_{i=0} (M(T_i) - M_i)^2, \quad \text{s.t. } 0 < r < R
$$

Target mass loss $M_i = [10, 15, 30]$ at $T_i = [5, 7, 30]$.

Optimized parameters: $\rho \approx 40.90, R \approx 0.50, r \approx 0.30$.

<p align="center">
  <img src="figures/Part3_Task2_Massloss.png" width="70%">
</p>

---

## Concluding Discussion

The investigation showed that initial density $\rho$, nonlinear parameters $\beta$ and $\gamma$, and torus dimensions $R, r$ all influence the diffusion process. Optimization successfully identified parameters to meet specific dosage requirements.

---

## Attached Code

### Problem A.2 (Matlab)

```matlab
a = -1;
b = 1;
gl = 0;
gr = 0;
alpha = 0.01;

beta = 0.9;
N = 11;
TOL = 1e-3;
Nmax = 1e4;
h = (b - a) / N;
x = a:h:b;

while true
    N = length(x);
    zeta = my_discrete_Lapl(x, @fStar, gl, gr);
    eta2 = zeros(N, 1);
    for j = 1 : N-1
        h = x(j+1)- x(j);
        a = f(x(j)) + alpha*zeta(j);
        b = f(x(j+1)) + alpha*zeta(j+1);
        t = (a.^2+b.^2)*h/2;
        eta2(j) = h.^2 * t;
    end
    
    if sum(eta2) <= TOL || N > Nmax
        break;
    end
    
    for k = 1 : (length(eta2)-1)
        if eta2(k) > beta*max(eta2)
            x = [x (x(k+1)+x(k))/2];
        end
    end
    x = sort(x);
end

A = my_stiffness_matrix_assemble(x);
B = my_load_vector_assemble(x, @fStar, gl, gr);
xi = A \ B;
Re = f(x) + alpha * zeta';
eta = eta2.^0.5;

figure;
subplot(2, 2, 1);
plot(x, xi);
title('Solution uh');
xlabel('x');
ylabel('uh');

% Plot the residual R(uh)
subplot(2, 2, 2);
plot(x, Re);
title('Residual R(uh)');
xlabel('x');
ylabel('R(uh)');

% Plot the error indicator eta(uh)
subplot(2, 2, 3);
plot(x, eta);
title('Error Indicator eta(uh)');
xlabel('x');
ylabel('eta(uh)');

% Plot the mesh-size distribution
subplot(2, 2, 4);
plot(x(2:end), 1./diff(x));
title('Mesh-size Distribution');
xlabel('x');
ylabel('1/dx');

disp(size(x));


function y = f(x)
    condition = ((x <= 0.8) & (x >= 0.2)) | ((x <= -0.2) & (x >= -0.8));
    y = zeros(size(x));
    y(condition) = 10;
end

function y = fStar(x)
    condition = ((x <= 0.8) & (x >= 0.2)) | ((x <= -0.2) & (x >= -0.8));
    y = zeros(size(x));
    y(condition) = 10 / 0.01;
end

function zeta = my_discrete_Lapl(x, ~, gl, gr)
    A = my_stiffness_matrix_assemble(x);
    B = my_load_vector_assemble(x, @fStar, gl, gr);
    xi = A \ B;
    M = my_mass_matrix_assemble(x);
    zeta = -inv(M) * A * xi; 
end

function A = my_stiffness_matrix_assemble(x)

    N = length(x) - 1;
    A = zeros(N+1, N+1); 
    
    for i = 1 : N
        h = x(i+1) - x(i);
        n = [i i+1];
        A(n,n) = A(n,n) + [1 -1; -1 1]/h;
    end
    
    % adjust for BC
    A(1,1) = 1;
    A(1,2) = 0;
    A(N+1,N) = 0;
    A(N+1,N+1) = 1;

end

function b = my_load_vector_assemble(x, f, gl, gr)
    %
    % Return assembled load vector b
    % Input vector x of node coords
    % 
    N = length(x) - 1;
    b = zeros(N+1, 1);
    for i = 1 : N
        h = x(i+1) - x(i);
        n = [i i+1];
        b(n) = b(n) + [f(x(i)); f(x(i+1))]*h/2;
    end
    b(1) = gl;
    b(N+1) = gr;

end

function M = my_mass_matrix_assemble(x)
    N = length(x)-1;
    M = zeros(N+1, N+1);
    
    for i = 1 : N
        h = x(i+1) - x(i);
        n = [i i+1];
        M(n,n) = M(n,n) + [1/3 1/6; 1/6 1/3]*h;
    end
end
```

### Problem B.1 (Matlab)

```matlab
clc
clear
g = @circleg;

H = [1/2 1/4 1/8 1/16 1/32];

energyNorm = zeros(length(H),1);
EnE = zeros(length(H),1);

for i = 1:length(H)
    hm = H(i);
        
    [p, e, t] = initmesh(g, 'hmax', hm);
    g1 = u(p(1, e(1,:)), p(2, e(1,:)));
    A = StiffnessAssembler2D(p, t, @(x, y) 1);
    b = LoadAssembler2D(p, t, @f);
    I = eye(length(p));
    
    A(e(1,:),:) = I(e(1,:),:);
    b(e(1,:))= g1;
  
    uh = A\b;

    err = u(p(1,:), p(2,:)) - uh';
    EnE(i) = sqrt(err * A * err');
    
    if i == 1
        figure;
        pdeplot(p, e, t, 'XYData', uh, 'ZData', uh);
    elseif i == length(H)
        figure;
        pdeplot(p, e, t, 'XYData', uh, 'ZData', uh);
    end
end

gamma = zeros((length(H)-1),1);

for i=1:length(energyNorm)-1
    gamma(i) = log(EnE(i+1)/EnE(i))/log(H(i+1)/H(i));
    disp(gamma(i));
end

figure;
loglog(H, EnE, 'r--', 'LineWidth', 1.5);
hold on;
loglog(H, H.^mean(gamma), 'b-', 'LineWidth', 1.5);
xlabel('h_{max}');
grid on;
legend('EnE', ['h_{max}^{\gamma}, \gamma = ' num2str(gamma(end))], 'Location', 'southeast');
set(gca, 'FontSize', 14);
set(legend, 'FontSize', 14);
hold off;

function z = f(x, y)
    z = 8 * pi^2 * sin(2 * pi * x) .* sin(2 * pi * y);
end

function z = u(x, y)
    z = sin(2 * pi * x) .* sin(2 * pi * y);
end
```

### Problem B.3 (Matlab)

```matlab
clc
clear

alpha = 0.01;

T = 30; % final time

dt = 0.1; % time step
tn = T / dt; % number of time levels
MLoss = zeros(tn+1, 2);

g = @circleg;

H = [1/5, 1/20];

for index = 1 : length(H)
    hmax = H(index); % mesh size
    [p, e, t] = initmesh(g, 'hmax', hmax);
    nt = size(t, 2); % number of elements
    
    A = StiffnessAssembler2D(p, t, @(x, y) alpha);
    M = MassAssembler2D(p, t);
    I = eye(length(p));
    
    old_xi = u0(p(1,:), p(2,:))';
    for i = 1: tn
        % impose Dirichlet Boundary Condition
        D = ((1/dt).*M + (1/2).*A);
        D(e(1,:),:) = I(e(1,:),:);

        b = ((1/dt).*M - (1/2).*A)*old_xi;
        b(e(1,:),:) = 0;
        
        xi = D \ b;
        % compute the mass loss
        for K = 1:nt
            triangle = t(1:3, K);
            x = p(1, triangle);
            y = p(2, triangle);
            [area, ~, ~] = HatGradients(x, y);
            MLoss(i+1, index) = MLoss(i+1, index) + area/3 * (sum(u0(x,y))- sum(xi(triangle)));
        end
        old_xi = xi;
    end
end

figure;
plot(linspace(0, T, tn + 1), MLoss(:,1), 'DisplayName', 'h_{max} = 1/5', 'LineWidth', 2, 'Color', 'b');
hold on;
plot(linspace(0, T, tn + 1), MLoss(:,2), 'DisplayName', 'h_{max} = 1/20', 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
hold off;
legend('show');
xlabel('Time (day)', 'Interpreter', 'latex');
ylabel('Concentration of Hormone ($mmol/mm^3$)', 'Interpreter', 'latex');

figure;
pdeplot(p, e, t, 'XYData', u0(p(1,:), p(2,:)), 'ZData', u0(p(1,:), p(2,:)));

figure;
pdeplot(p, e, t, 'XYData', xi, 'ZData', xi);
```

### Problem B.4 (Matlab)

```matlab
clc
clear

alpha = 0.01;
beta = 1;
gamma = 0.2;

T = 30; % final time

dt = 0.1; % time step
tn = T / dt; % number of time levels
MLoss = zeros(tn+1, 2);

g = @circleg;

H = [1/5, 1/20];

for index = 1 : length(H)
    hmax = H(index); % mesh size
    [p, e, t] = initmesh(g, 'hmax', hmax);
    nt = size(t, 2); % number of elements
    
    A = StiffnessAssembler2D(p, t, @(x, y) alpha);
    M = MassAssembler2D(p, t);
    I = eye(length(p));
    
    A(e(1,:),:) = I(e(1,:),:);

    old_xi = u0(p(1,:), p(2,:))';
    for i = 1: tn
        D = ((1/dt).*M + (1/2).*A);
        D(e(1,:),:) = I(e(1,:),:);
        b = (((1/dt).*M - (1/2).*A)*old_xi + beta.*M*old_xi - beta*gamma.*M*(old_xi.*old_xi));
        b(e(1,:),:) = 0;
        xi = D \ b;

        % compute the mass loss
        for K = 1:nt
            triangle = t(1:3, K);
            x = p(1, triangle);
            y = p(2, triangle);
            [area, ~, ~] = HatGradients(x, y);
            MLoss(i+1, index) = MLoss(i+1,index) + area/3 * (sum(u0(x,y))- sum(xi(triangle)));
        end
        old_xi = xi;
    end
end

figure;
plot(linspace(0, T, tn + 1), MLoss(:,1), 'DisplayName', 'h_{max} = 1/5', 'LineWidth', 2, 'Color', 'b');
hold on;
plot(linspace(0, T, tn + 1), MLoss(:,2), 'DisplayName', 'h_{max} = 1/20', 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
hold off;
legend('show');
xlabel('Time (day)', 'Interpreter', 'latex');
ylabel('Concentration of Hormone ($mmol/mm^3$)', 'Interpreter', 'latex');

figure;
pdeplot(p, e, t, 'XYData', u0(p(1,:), p(2,:)), 'ZData', u0(p(1,:), p(2,:)));

figure;
pdeplot(p, e, t, 'XYData', xi, 'ZData', xi);
```

### Functions for 2D Problems (Matlab)

```matlab
function y = u0(x1, x2)
    condition = ((x1.^2 + x2.^2) <= 0.64 & (x1.^2 + x2.^2) >= 0.04);
    y = zeros(size(x1));
    y(condition) = 10;
end

function A = StiffnessAssembler2D(p, t, a)
    np = size(p, 2);
    nt = size(t, 2);
    A = sparse(np, np);
    for K = 1:nt
        loc2glb = t(1:3, K);
        x = p(1, loc2glb);
        y = p(2, loc2glb);
        [area, b, c] = HatGradients(x, y);
        xc = mean(x);
        yc = mean(y);
        abar = a(xc, yc);
        AK = abar * ( b*b' + c*c') * area;
        A(loc2glb, loc2glb) = A(loc2glb, loc2glb) + AK;
    end
end

function M = MassAssembler2D(p,t)
    np = size(p,2);
    nt = size(t,2);
    M = sparse(np,np);
    for K = 1:nt
        loc2glb = t(1:3, K);
        x = p(1, loc2glb);
        y = p(2, loc2glb);
        area = polyarea(x,y);
        MK = [2 1 1; 1 2 1; 1 1 2]/12*area;
        M(loc2glb, loc2glb) = M(loc2glb, loc2glb) + MK; 
    end
end

function b = LoadAssembler2D(p,t,f)
    np = size(p,2);
    nt = size(t,2);
    b = zeros(np, 1);
    for K = 1:nt
        loc2glb = t(1:3, K);
        x = p(1, loc2glb);
        y = p(2, loc2glb);
        area = polyarea(x,y);
        bK = [f(x(1), y(1)); f(x(2), y(2)); f(x(3), y(3))]/3*area;
        b(loc2glb) = b(loc2glb) + bK;
    end
end

function [area, b, c] = HatGradients(x, y)
    area = polyarea(x, y);
    b = [y(2)-y(3); y(3)-y(1); y(1)-y(2)] / 2 / area;
    c = [x(3)-x(2); x(1)-x(3); x(2)-x(1)] / 2 / area;
end
```

### Problem C.1 2D Mesh (Python)

```python
from dolfin import *

mesh = Mesh("circle1.xml")
Q = FunctionSpace(mesh, "CG", 1)  

T = 20
h = mesh.hmin()
dt = h
alpha = 0.01

class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

g = Constant(0.0)
bc = DirichletBC(Q, g, DirichletBoundary())

indata = Expression("abs(R-sqrt(pow(x[0], 2) + pow(x[1], 2))) <= r ? rho : 0", R=0.5, r=0.2, rho=10, degree=2)
u0 = Function(Q)
u0 = interpolate(indata, Q)

u = TrialFunction(Q)
v = TestFunction(Q)
a = u*v*dx + alpha*dt*inner(grad(u),grad(v))*dx
L = u0*v*dx

file = File("output_2D/2D_part3_Task1.pvd")
u = Function(Q)
t = 0
num_steps = int(T / dt)
for n in range(num_steps):
    t += dt
    solve(a == L, u, bc)
    file << (u, t)
    u0.assign(u)
```

### Problem C.1 3D Mesh (Python)

```python
from dolfin import *

mesh = Mesh("sphere1.xml")
Q = FunctionSpace(mesh, "CG", 1)  

T = 20
h = mesh.hmin()
dt = h
alpha = 0.01

class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

g = Constant(0.0)
bc = DirichletBC(Q, g, DirichletBoundary())

indata = Expression("pow((R-sqrt(pow(x[0], 2) + pow(x[1], 2))),2) + pow(x[2], 2) <= pow(r,2) ? rho : 0", R=0.5, r=0.2, rho=10, degree=2)
u0 = Function(Q)
u0 = interpolate(indata, Q)

u = TrialFunction(Q)
v = TestFunction(Q)
a = u*v*dx + alpha*dt*inner(grad(u),grad(v))*dx
L = u0*v*dx

file = File("output_3D/3D_part3_Task1.pvd")
u = Function(Q)
t = 0
num_steps = int(T / dt)
for n in range(num_steps):
    t += dt
    solve(a == L, u, bc)
    file << (u, t)
    u0.assign(u)
```

### Problem C.2 (Python)

```python
from dolfin import *

mesh = Mesh("sphere1.xml")
Q = FunctionSpace(mesh, "CG", 1)  

T = 50
h = mesh.hmin()
dt = h
alpha = 0.01

class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

g = Constant(0.0)
bc = DirichletBC(Q, g, DirichletBoundary())

rho = [10, 20, 40]
for i in range(len(rho)):
    indata = Expression("pow((R-sqrt(pow(x[0], 2) + pow(x[1], 2))),2) + pow(x[2], 2) <= pow(r,2) ? rho : 0", R=0.5, r=0.2, rho=rho[i], degree=2)
    u0 = Function(Q)
    u0 = interpolate(indata , Q)

    u = TrialFunction(Q)
    v = TestFunction(Q)
    a = u*v*dx + alpha*dt*inner(grad(u),grad(v))*dx
    L = u0*v*dx

    file = File(f"output_Task2/rho_{rho[i]}.pvd")
    u = Function(Q)
    u_initial = Function(Q)
    u_initial = interpolate(indata , Q)
    M = (u_initial - u) * dx
    mass_loss = []
    t = 0
    num_steps = int(T / dt)
    for n in range(num_steps):
        t += dt
        solve(a==L,u,bc)
        mass = assemble(M)
        mass_loss.append(mass)
        file << (u, t)
        u0.assign(u)

    with open("mass_loss_{}.txt".format(rho[i]), "w") as f:
        for m in mass_loss:
            f.write(str(m) + "\n")
```

### Problem C.3 (Python)

```python
import scipy.optimize as optimize
from scipy.optimize import minimize
from dolfin import *

def func(data):
    mesh = Mesh("sphere1.xml")
    Q = FunctionSpace(mesh, "CG", 1)  
    T = 30
    dt = 0.1
    alpha = 0.01
    R = data[1]
    r = data[2]
    rho = data[0]

    class DirichletBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary

    g = Constant(0.0)
    bc = DirichletBC(Q, g, DirichletBoundary())

    indata = Expression("pow((R-sqrt(pow(x[0], 2) + pow(x[1], 2))),2) + pow(x[2], 2) <= pow(r,2) ? rho : 0", R=R, r=r, rho=rho, degree=1)
    u0 = Function(Q)
    u0 = interpolate(indata, Q)

    u = TrialFunction(Q)
    v = TestFunction(Q)
    a = u*v*dx + alpha*dt*inner(grad(u),grad(v))*dx
    L = u0*v*dx

    u = Function(Q)
    u_initial = Function(Q)
    u_initial = interpolate(indata , Q) 
    M = (u_initial - u) * dx
    Mt = [0,0,0]
    t = 0
    while t < T:
        t = round(t + dt, 2)
        solve(a==L,u,bc)
        mass = assemble(M)
        if t == 5:
            Mt[0] = mass
        elif t == 7:
            Mt[1] = mass
        elif t == 30:
            Mt[2] = mass
        u0.assign(u)
    F = (Mt[0] - 10)**2 + (Mt[1] - 15)**2 + (Mt[2] - 30)**2
    return F

data = [20, 0.5, 0.1]
res = minimize(func , data, method='nelder-mead', options={'xtol':1e-3, 'disp': True})
print(res)
```
