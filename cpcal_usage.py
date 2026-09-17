from cpcal import cp_air
temperature = 300.0  # Example temperature in Kelvin
cp_value = cp_air(temperature)
print(f"The specific heat capacity of dry air at {temperature} K is {cp_value:.2f} J/(kg*K).")