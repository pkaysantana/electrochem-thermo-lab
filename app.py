import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from constants import FARADAY_C_MOL, R_J_MOL_K, C_WATER_J_G_K, CELL_N

def pathway_a_potentiometry(T_pot_K, E_pot_V, cell_name="Daniell"):
    """
    Pathway A: Thermodynamic Calculator from Potentiometry.
    ΔG = -nFE
    ΔS = nF(dE/dT)
    ΔH = ΔG + TΔS
    """
    n_elec = CELL_N.get(cell_name, 2)
    
    # Run linear regression: E vs T
    slope, intercept, r_value, p_value, std_err = stats.linregress(T_pot_K, E_pot_V)
    r_squared = r_value ** 2
    
    # Advisory Watchdog
    if r_squared < 0.95:
        print(f"\n[ADVISORY WATCHDOG] POLARIZATION WARNING for {cell_name} (R²={r_squared:.4f}):")
        print("Voltage drop detected. Ensure you are taking readings instantly.")

    # Thermodynamic Sanity Check
    if slope > 0:
        print(f"\n[CRITICAL THERMODYNAMIC WARNING]: Positive slope detected (\u0394S > 0).")
        print("This violates the expected physical chemistry for this cell.")
        print("Verify that your T and E arrays are not inverted and that")
        print("polarization didn't skew your higher-temperature readings.")
        
    # Energy calculations
    T_mean = np.mean(T_pot_K)
    E_mean = slope * T_mean + intercept
    
    delta_G_J_mol = -n_elec * FARADAY_C_MOL * E_mean
    delta_S_J_mol_K = n_elec * FARADAY_C_MOL * slope
    delta_H_pot_J_mol = delta_G_J_mol + T_mean * delta_S_J_mol_K
    
    # Plotting (16:9, High-DPI)
    plt.figure(figsize=(16, 9), dpi=300)
    plt.scatter(T_pot_K, E_pot_V, color='blue', label='Experimental Data', s=100)
    
    T_fit = np.linspace(min(T_pot_K)-5, max(T_pot_K)+5, 100)
    E_fit = slope * T_fit + intercept
    
    # Trendline in Scientific Notation (3 decimal places)
    slope_sci = f"{slope:.3e}"
    intercept_sci = f"{intercept:.3e}"
    trendline_label = f"Trendline: E = {slope_sci} * T + {intercept_sci}\n$R^2$ = {r_squared:.4f}"
    
    plt.plot(T_fit, E_fit, color='red', linestyle='--', label=trendline_label)
    plt.xlabel("Temperature (K)", fontsize=14)
    plt.ylabel("Cell Potential (V)", fontsize=14)
    plt.title(f"Pathway A: {cell_name} Cell Potential vs Temperature", fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.tight_layout()
    plot_filename = f"{cell_name.replace('/', '_')}_pot_plot.png"
    plt.savefig(plot_filename)
    plt.close()
    
    return delta_G_J_mol, delta_S_J_mol_K, delta_H_pot_J_mol, plot_filename

def pathway_b_calorimetry(mass_sol_g, delta_T_cal_K, moles_rxn):
    """
    Pathway B: Thermodynamic Calculator from Calorimetry.
    q = mcΔT
    ΔH_cal = -q / moles_rxn
    """
    q_J = mass_sol_g * C_WATER_J_G_K * delta_T_cal_K
    delta_H_cal_J_mol = -q_J / moles_rxn
    return delta_H_cal_J_mol

def main():
    """
    =======================================================================
                             DUAL DATA FLOW ARCHITECTURE
    =======================================================================
          [ constants.py ]
          + F = 96485.3 C/mol
          + R = 8.314 J/(mol K)
          + c_water = 4.184 J/(g K)
          + cell_n = {"Daniell": 2, "Zn/Ag": 2, "Cu/Ag": 2}
                  |
                  v
    ====================================
               [ app.py ]
    ====================================
                 /  \
                /    \
    [PATHWAY A: POTENTIOMETRY]      [PATHWAY B: CALORIMETRY]
     Inputs: T_pot_K, E_pot_V        Inputs: mass_sol_g, c_cal_J_gK,
             n_elec                          delta_T_cal_K, moles_rxn
        |                               |
        |-- scipy.stats.linregress      |-- q = mcΔT
        |-- R^2 < 0.95 Check -> Warn    |-- ΔH_cal = -q / moles_rxn
        |-- plot(E vs T, 16:9, DPI)     |
        v                               v
     Outputs: ΔG, ΔS, ΔH_pot          Outputs: ΔH_cal
    """
    print("=======================================================================")
    print("                      ELECTROCHEMICAL CELLS TOOLKIT                    ")
    print("                         DUAL DATA FLOW EXECUTION                      ")
    print("=======================================================================")
    
    # --- Mock Daniell Cell Data (chemically accurate: E decreases as T rises, ΔS < 0) ---
    # T in °C: [15, 20, 25, 30, 35, 40] → converted to Kelvin
    T_mock_K = np.array([288.15, 293.15, 298.15, 303.15, 308.15, 313.15])
    E_mock_V = np.array([1.104, 1.102, 1.100, 1.098, 1.095, 1.093])
    
    # 1. Pathway A (Good Data)
    print("\n--- PATHWAY A (POTENTIOMETRY) ---")
    dG, dS, dH_pot, plot_file = pathway_a_potentiometry(T_mock_K, E_mock_V, "Daniell")
    print(f"Calculated ΔG: {dG/1000:,.2f} kJ/mol")
    print(f"Calculated ΔS: {dS:,.2f} J/(mol*K)")
    print(f"Calculated ΔH (Potentiometry): {dH_pot/1000:,.2f} kJ/mol")
    print(f"Plot saved successfully as '{plot_file}'")
    
    # 2. Watchdog Test (Bad Data)
    print("\n--- PATHWAY A (POTENTIOMETRY) - BAD DATA TEST ---")
    E_bad_V = np.array([1.100, 1.050, 1.090, 1.010, 1.105]) # Will trigger R² < 0.95
    pathway_a_potentiometry(T_mock_K, E_bad_V, "Daniell_Bad")
    
    # 3. Pathway B (Calorimetry)
    print("\n--- PATHWAY B (CALORIMETRY) ---")
    mass_mock_g = 100.0        # 100 g solution
    delta_T_mock_K = 5.0       # 5 K rise
    moles_mock = 0.01          # 0.01 moles reacting
    dH_cal = pathway_b_calorimetry(mass_mock_g, delta_T_mock_K, moles_mock)
    print(f"Calculated ΔH (Calorimetry): {dH_cal/1000:,.2f} kJ/mol")
    print("\nExecution complete.")

if __name__ == "__main__":
    main()
