# gravity-inversion-2d-3d-jacobian
2D and 3D gravity anomaly modeling using nonlinear inversion with a linearized Jacobian matrix approach in MATLAB.
Gravity Anomaly Modeling and Density Inversion (2D & 3D)
📌 Repository Description

This repository contains MATLAB programs for 2D and 3D gravity anomaly modeling and non-linear inversion using a linearized Jacobi matrix approach. The primary objective is to estimate subsurface model parameters—particularly density (ρ)—based on synthetic gravity data generated from a spherical mass model.

The codes are designed for educational and research purposes in geophysics, especially for understanding gravity forward modeling and inverse problem techniques.

🎯 Objectives

Simulate 2D and 3D gravity anomalies caused by a buried spherical body

Apply non-linear inversion using Jacobian matrix linearization

Estimate subsurface parameters:

Horizontal position (x)

Depth (z)

Density contrast (ρ)

Analyze misfit convergence during inversion iterations

📂 Repository Structure
├── gravity_inversion_2D_3D.m     % Full 2D & 3D gravity modeling and inversion
├── density_inversion_rho.m      % Program focused only on density (rho) estimation
├── README.md                    % Repository documentation

🧮 Methodology Overview
1. Forward Gravity Modeling

The subsurface is represented by a homogeneous sphere

Gravity response is calculated using Newton’s gravitational law:

𝑔
=
𝐺
⋅
4
3
𝜋
𝑅
3
𝜌
⋅
𝑧
𝑟
3
g=G⋅
3
4
	​

πR
3
ρ⋅
r
3
z
	​


where:

𝐺
G = gravitational constant

𝑅
R = radius of the sphere

𝜌
ρ = density contrast

𝑟
r = distance between observation point and model center

2. Non-Linear Inversion (Jacobian Matrix Approach)

The inversion problem is solved iteratively by linearizing the forward model:

Δ
𝑑
=
𝐽
⋅
Δ
𝑚
Δd=J⋅Δm

where:

𝐽
J = Jacobian matrix (partial derivatives)

Δ
𝑚
Δm = model parameter updates

Δ
𝑑
Δd = data misfit

A damped least-squares solution is applied for stability.

3. Jacobian Matrix Components
🔹 2D Inversion

Model parameters:

𝑥
0
x
0
	​

 : horizontal position

𝑧
0
z
0
	​

 : depth

𝜌
ρ : density

Jacobian matrix:

𝐽
2
𝐷
=
[
∂
𝑔
∂
𝑥
0
	
∂
𝑔
∂
𝑧
0
	
∂
𝑔
∂
𝜌
]
J
2D
	​

=[
∂x
0
	​

∂g
	​

	​

∂z
0
	​

∂g
	​

	​

∂ρ
∂g
	​

	​

]
🔹 3D Inversion

Model parameters:

𝑥
0
,
𝑦
0
,
𝑧
0
x
0
	​

,y
0
	​

,z
0
	​


Jacobian matrix:

𝐽
3
𝐷
=
[
∂
𝑔
∂
𝑥
0
	
∂
𝑔
∂
𝑦
0
	
∂
𝑔
∂
𝑧
0
]
J
3D
	​

=[
∂x
0
	​

∂g
	​

	​

∂y
0
	​

∂g
	​

	​

∂z
0
	​

∂g
	​

	​

]
🧪 Density-Only Inversion Program

The density_inversion_rho.m script is a simplified inversion case where:

Geometry parameters are updated

The main focus is estimating density (ρ)

Suitable for understanding parameter sensitivity and inversion stability

The final estimated density is displayed as:

Estimated density (rho): XX.XX kg/m³

📊 Outputs

2D gravity anomaly curves (observed vs calculated)

3D gravity anomaly surfaces

Subsurface model visualization

Misfit convergence plots (2D & 3D)

Estimated subsurface density value

🛠 Requirements

MATLAB R2018a or newer

No additional toolboxes required

🚀 How to Run

Clone the repository:

git clone https://github.com/your-username/your-repo-name.git


Open MATLAB and set the repository folder as the working directory

Run:

gravity_inversion_2D_3D


or

density_inversion_rho

📚 Applications

Gravity data interpretation

Geophysical inversion learning

Mineral and geothermal exploration modeling

Academic demonstrations of inverse problems

⚠️ Notes

The data used are synthetic

Results depend on initial model parameters

Regularization is essential for inversion stability

👤 Author

M. Ilham
Geophysics | Gravity Modeling & Inversion

⭐ License

This project is open for academic and educational use.
Please cite appropriately if used in publications.
