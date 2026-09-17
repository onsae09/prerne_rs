from cpcal import cp_air
import matlab.engine

eng = matlab.engine.start_matlab()
temp = int(input("Enter the temperature in Kelvin: "))
C_p = cp_air(temp)

eng.navier_stokes_solver(C_p, nargout=0)  # Call the MATLAB function with cp as an argument