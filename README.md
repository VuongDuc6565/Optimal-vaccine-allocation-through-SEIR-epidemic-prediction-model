# Optimal Vaccine Allocation with a Multi-City SEIR Model (CasADi + IPOPT)

This repository contains a Python implementation of an **optimal control model** that distributes a limited vaccine stock among multiple Vietnamese provinces in order to **minimise the number of infections after _T_ days**.  
The epidemiological dynamics are captured by an **SEIR** (Susceptible–Exposed–Infectious–Recovered) compartmental model, while the optimisation is handled with **CasADi** and solved using the **IPOPT** nonlinear solver.

---

## Data Input

The model reads an Excel file `data-1.xlsx`, sheet **`covid19_provinces_vn`**, which must contain at least the following columns:

| Column                     | Meaning                                   |
|----------------------------|-------------------------------------------|
| `Province`                 | Province name (63 rows expected)          |
| `N`, `S`, `E`, `I`, `R`    | Initial population & compartment values   |
| `Gamma`, `Alpha`, `Beta`, `Nu`, `Mu` | SEIR parameters per province     |
| `Provincial_Hospital`      | Number of provincial hospitals            |
| `P+ICP`                    | Polyclinic & inter-communal posts         |
| `Maternity_home`           | Maternity homes                           |
| `C.H.C`                    | Commune health centres                    |
| `Vaccine_min`              | Minimum cumulative vaccines for province  |

---

## Model Parameters

| Symbol | Description | This code |
|--------|-------------|-----------|
| `T` | Planning horizon in days | `T = 30` |
| `ω` (omega) | Vaccine efficacy | `omega = 0.75` |
| `V_total` | National vaccine budget | `V_total = 20_000_000` |
| `C[c]` | Daily vaccination capacity of province *c*  | computed from health-facility counts |

**Decision variables**

```
S[c,t], E[c,t], I[c,t], R[c,t]   # SEIR states
v[c,t]                           # Vaccine doses delivered to province c on day t
```

**Daily capacity constraint**

```
0 ≤ v[c,t] ≤ C[c]
```

**Minimum provincial allocation**

```
SUM_t v[c,t] ≥ Vaccine_min[c]
```

**National budget**

```
SUM_c SUM_t v[c,t] ≤ V_total
```

**Epidemic dynamics (per province)**

```
S(t+1) = S(t) + ν·N(t) - β·S(t)·I(t)/N(t) - ω·v(t+1) - μ·S(t)
E(t+1) = E(t) + β·S(t)·I(t)/N(t) - (α+μ)·E(t)
I(t+1) = I(t) + α·E(t) - (γ+μ)·I(t)
R(t+1) = R(t) + γ·I(t) + ω·v(t+1) - μ·R(t)
N(t+1) = (1+ν-μ)·N(t)
```

**Objective: Minimize total infections at final day**

```
min  Σ_c  I[c,T]
```

---

## Quick Start

```bash
# 1. Install packages
pip install casadi pandas xlrd xlsxwriter

# 2. Run the optimisation
python casadi_seir.py 

# 3. Results
# SEIR_results.xlsx  (Summary + Daily_Vaccine sheets)

```

---

## Output Files

| File | Content |
|------|---------|
| `SEIR_results.xlsx` | **Summary** sheet: infections & total vaccine per province<br>**Daily_Vaccine** sheet: doses delivered each day |
| `schedule.log` | IPOPT solver iterations, objective value, convergence info |

---

## Model Background

A concise derivation of the constraints and the optimisation model can be found in *“Phân phối tối ưu vaccine thông qua mô hình dự đoán dịch bệnh SEIR”* (PDF included). It details:

* Classic SEIR dynamics (Kermack & McKendrick, 1927)  
* Capacity calculation `C[c] = 5000·Provincial_Hospital + 1500·P+ICP + 1000·Maternity_home + 500·C.H.C`  
* Logical constraints ensuring realistic vaccine roll-out

For a deeper dive, open **`SEIR.pdf`** in the repository.

---

## License

MIT – feel free to adapt the model to other regions or different vaccine parameters.
