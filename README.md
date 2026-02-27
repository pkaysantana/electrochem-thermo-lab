# Experiment 2: Electrochemical Cells Toolkit

**COSHH: Zinc Amalgam electrodes contain toxic mercury. DO NOT abrade with steel wool.**

## Thermodynamic Derivations

This Python toolkit executes thermodynamic calculations across two isolated measurement pathways: Pathway A (Electrical) and Pathway B (Thermal).

### Pathway A: Potentiometry

By measuring cell potential ($E$) as a monotonic function of uniform absolute temperature ($T$):

1. **Gibbs Free Energy ($\Delta G$)**
   $$ \Delta G = -nFE $$
   Where $n$ is the number of electrons transferred, $F$ is Faraday's constant, and $E$ is the mean cell potential.

2. **Entropy ($\Delta S$)**
   $$ \Delta S = nF \left( \frac{\partial E}{\partial T} \right)_P $$
   Where $\frac{\partial E}{\partial T}$ is derived functionally as the slope of the linear regression traversing $E$ against $T$.

3. **Enthalpy ($\Delta H$)**
   $$ \Delta H = \Delta G + T\Delta S $$

### Pathway B: Calorimetry

By establishing a localized macroscopic change in uniform temperature ($\Delta T$) surrounding the solution:

1. **Heat Transferred ($q$)**
   $$ q = mc\Delta T $$
   Where $m$ is the total inclusive mass of the solution and $c$ is the specific isobaric heat capacity.

2. **Molar Enthalpy ($\Delta H_{cal}$)**
   $$ \Delta H_{cal} = -\frac{q}{\text{moles of limiting reactant}} $$

## Getting Started

```bash
# Setup (if necessary)
pip install numpy scipy matplotlib

# Run the Toolkit and Advisory Watchdog
python app.py
```
