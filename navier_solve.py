from pathlib import Path

from cpcal import cp_air
import matlab.engine

def main() -> None:
    temp = float(input("Enter the temperature in Kelvin: "))
    cp = float(cp_air(temp))

    eng = matlab.engine.start_matlab()
    try:
        # MATLAB does not automatically search the directory containing this
        # Python file when the script is launched from another directory.
        matlab_dir = str(Path(__file__).resolve().parent)
        eng.addpath(matlab_dir, nargout=0)

        result = eng.navier_stokes_solver(cp, float(temp), nargout=1)
        print(f"Cp = {cp:.3f} J/(kg*K)")
        print(f"mu = {float(result['mu']):.6e} Pa*s")
        print(f"k = {float(result['k']):.6e} W/(m*K)")
        print(f"alpha = {float(result['alpha']):.6e} m^2/s")
        print(f"final mean temperature = {float(result['mean_temperature']):.3f} K")
    finally:
        eng.quit()


if __name__ == "__main__":
    main()
