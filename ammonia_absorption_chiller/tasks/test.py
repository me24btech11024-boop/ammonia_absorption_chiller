from core.moist_air import MoistAir

# Room: 22°C, 75% RH
w_room = MoistAir.humidity_ratio(22.0, 75.0)  # ~0.01246 kg/kg
h_room = MoistAir.enthalpy(22.0, w_room)       # ~53.76 kJ/kg

# Outdoor: 40°C, 50% RH
w_out = MoistAir.humidity_ratio(40.0, 50.0)   # ~0.02341 kg/kg
h_out = MoistAir.enthalpy(40.0, w_out)         # ~100.67 kJ/kg

# 70/30 Mixing
T_mix = 0.7 * 22.0 + 0.3 * 40.0
w_mix = 0.7 * w_room + 0.3 * w_out
h_mix = 0.7 * h_room + 0.3 * h_out
T_mix_exact = 0.7 * 22.0 + 0.3 * 40.0  # Returns exactly 27.400000000000002 °C

print(f"T_mix: {T_mix:.2f} °C (Target: 27.48 °C)")
print(f"T_mix_exact: {T_mix_exact:.2f} °C")
print(f"w_mix: {w_mix:.5f} kg/kg (Target: 0.01575 kg/kg)")
print(f"h_mix: {h_mix:.2f} kJ/kg (Target: 67.84 kJ/kg)")