from cpcal import cp_air
import matlab.engine

eng = matlab.engine.start_matlab()
temp = int(input("Enter the temperature in Kelvin: "))
cp = cp_air(temp)

eng.calculate(cp, nargout=0)  # Call the MATLAB function with cp as an argument