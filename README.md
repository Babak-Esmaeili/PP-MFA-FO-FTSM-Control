# Predefined Performance-Based Model-Free Adaptive Fractional-Order Fast Terminal Sliding-Mode Control of MIMO Nonlinear Systems

This repository contains MATLAB code for simulating the examples presented in our paper:

📄 **Title**: Predefined Performance-Based Model-Free Adaptive Fractional-Order Fast Terminal Sliding-Mode Control of MIMO Nonlinear Systems  
📰 **Journal**: ISA Transactions, 2023  
🔗 [DOI: 10.1016/j.isatra.2022.05.036](https://doi.org/10.1016/j.isatra.2022.05.036)

---

## 🧠 Abstract

We propose a novel model-free adaptive controller that integrates fractional-order fast terminal sliding-mode control with predefined performance constraints. The controller is designed for uncertain nonlinear MIMO systems and ensures convergence within a finite time without requiring system models. Simulation results validate the effectiveness of the method in both SISO and MIMO systems, including a simulated 2-DOF robotic manipulator.

---

## 📁 Simulation Code Structure

The `codes/` folder contains three simulation environments:

- `codes/example1_siso/` – MATLAB simulation for a SISO nonlinear uncertain system
- `codes/example2_mimo/` – MATLAB simulation for a MIMO nonlinear uncertain system
- `codes/example3_robot_sim/` – MATLAB simulation for a 2-DOF robotic manipulator

Each folder includes `main.m` and supporting scripts to reproduce the results in the paper.

---

## 🛠 Requirements

- MATLAB R2018b or newer
- No commercial solvers are required

---

## 🛠 Usage

To run a simulation:

1. Download the repository from GitHub **or** clone it via Git (see below)
2. Open the desired folder in MATLAB
3. Open and run the `main.m` file
4. (Optional) Adjust control parameters or initial conditions
5. The script will automatically generate plots for tracking and convergence

---

### 🔁 Optional: Clone via Git

If you're familiar with Git, you can clone the repository directly:

```bash
git clone https://github.com/Babak-Esmaeili/pp-mfa-fo-ftsmc.git
cd pp-mfa-fo-ftsmc/codes/example1_siso/
```

Then open MATLAB in that folder and run `main.m`.

---

## 📜 License and Contact Info

This project is licensed under the MIT License – see the LICENSE file for details.  
You can customize the parameters and use them for your specific control system applications.

If you have any questions or encounter issues, please feel free to contact me.

📧 Babak Esmaeili – esmaeil1@msu.edu

---

## 📚 Citation

If you found this repository useful in your research, please cite it as:

```bibtex
@article{esmaeili2022predefined,
  title={Predefined performance-based model-free adaptive fractional-order fast terminal sliding-mode control of MIMO nonlinear systems},
  author={Esmaeili, Babak and Salim, Mina and Baradarannia, Mahdi},
  journal={ISA transactions},
  volume={131},
  pages={108--123},
  year={2022},
  publisher={Elsevier},
  doi={10.1016/j.isatra.2022.05.036}
}
```
