import pandas as pd
import numpy as np

# Load posterior samples from allesfitter CSV file
data = pd.read_csv('allelclass_eccentric_values.csv')

# === CONSTANTS ===
R_star_rsun = 1.705            # Stellar radius in solar radii (given)
M_star_msun = 1.33                # Stellar mass in solar masses (given)

Rsun_to_Rjup = 9.73116          # Conversion: Solar radius to Jupiter radius
Rjup_to_cm = 7.1492e9             # Jupiter radius in cm
Mjup_to_g = 1.899e30              # Jupiter mass in grams
Msun_to_kg = 1.989e30           # Solar mass in kilograms
G_SI = 6.67430e-11                # Gravitational constant (m³/kg/s²)
sec_per_day = 86400               # Seconds per day

# === DERIVED VALUES ===

# 1. Calculate companion radius in Jupiter radii
radius_rjup = data['b_rr'] * R_star_rsun * Rsun_to_Rjup
data['radius_rjup'] = radius_rjup

# 2. Convert radius to centimeters for volume calculation
radius_cm = radius_rjup * Rjup_to_cm
volume_cm3 = (4/3) * np.pi * radius_cm**3  # Volume of sphere in cm³

# 3. Calculate eccentricity from fitted parameters
ecc = np.sqrt(data['b_f_c']**2 + data['b_f_s']**2)

# 4. Calculate sin(i)^2 from cosine of inclination
sin2i = 1 - data['b_cosi']**2

# 5. Semi-major axis calculation:
#    a/R_star (dimensionless) is given by 'b_rsuma'
a_over_rstar = data['b_rsuma']

# Convert stellar radius from solar radii to meters
R_star_m = R_star_rsun * 6.957e8  # Solar radius in meters

# Calculate semi-major axis in meters
a_m = a_over_rstar * R_star_m

# 6. Radial velocity semi-amplitude in m/s (convert from km/s)
K_ms = data['b_K'] * 1000  # km/s to m/s

# 7. Stellar mass in kilograms
M_star_kg = M_star_msun * Msun_to_kg

# 8. Orbital period in seconds
P_days = data['b_period']
P_sec = P_days * sec_per_day

# 9. Calculate M_p * sin(i) in kg using RV formula approximation
#    Formula: K = (2πG / P)^{1/3} * (M_p sin i) / (M_star)^{2/3} / sqrt(1 - e²)
#    Rearranged: M_p sin i = K * sqrt(1 - e²) * (M_star)^{2/3} * P^{1/3} / (2πG)^{1/3}
factor = (K_ms * np.sqrt(1 - ecc**2) * M_star_kg**(2/3) * P_sec**(1/3)) / ( (2 * np.pi * G_SI)**(1/3) )
Mpsini_kg = factor

# 10. Correct for inclination (sin i) to get true companion mass in kg
sin_i = np.sqrt(sin2i)
mass_kg = Mpsini_kg / sin_i

# 11. Convert companion mass to grams for Jupiter mass conversion
mass_g = mass_kg * 1000  # kg to grams

# 12. Convert companion mass to Jupiter masses
mass_mjup = mass_g / Mjup_to_g
data['mass_mjup'] = mass_mjup

# 13. Calculate density in g/cm³ = mass (g) / volume (cm³)
density_gcm3 = mass_g / volume_cm3
data['density_gcm3'] = density_gcm3

# === PRINT SUMMARY STATISTICS ===




import numpy as np

def summarize(arr, name):
    p50 = np.percentile(arr, 50)
    p16 = np.percentile(arr, 16)
    p84 = np.percentile(arr, 84)
    
    # Format with 2 significant figures using general format
    def fmt(x):
        return f"{x:.3f}"
    
    print(f"{name}: {fmt(p50)} +{fmt(p84 - p50)} -{fmt(p50 - p16)}")



summarize(radius_rjup, "Radius [Rjup]")
summarize(mass_mjup, "Mass [Mjup]")
summarize(density_gcm3, "Density [g/cm³]")

# Optional: Save the dataframe with computed columns
# data.to_csv("posterior_with_correct_mass_density.csv", index=False)

