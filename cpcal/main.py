import argparse
import math
from pathlib import Path
from typing import Any, Dict, Union

R_UNIVERSAL = 8314.46261815324  # J/(kmol*K)


def _parse_float(x: str) -> float:
    return float(
        x.strip()
        .replace("D", "E")
        .replace("d", "e")
    )


def load_species_from_thermo(
    path: Union[str, Path], target_species: str
) -> Dict[str, Any]:
    lines = Path(path).read_text(
        encoding="utf-8",
        errors="ignore"
    ).splitlines()

    # Find THERMO header
    start = None

    for i, line in enumerate(lines):
        s = line.strip().lower()

        if s.startswith("ther"):
            # Next line = global temperature ranges
            start = i + 2
            break

    if start is None:
        raise ValueError("THERMO header not found")

    i = start

    while i < len(lines):

        # Skip empty/comment lines
        if not lines[i].strip():
            i += 1
            continue

        if lines[i].lstrip().startswith(("!", "#")):
            i += 1
            continue

        # thermo.inp contains multiple sections.  In particular, the fitted
        # dry-air entry appears after "END PRODUCTS", in the reactants-only
        # section, so an END marker is a section delimiter rather than
        # necessarily the end of the file.
        if lines[i].strip().lower().startswith("end"):
            i += 1
            continue

        # NASA format:
        # (a15, a65)
        name = lines[i][:15].strip()

        if i + 1 >= len(lines):
            break

        meta = lines[i + 1]

        try:
            # NASA:
            # (i2,1x,a6,1x,5(a2,f6.2),i2,f13.5,f15.3)

            n_ranges = int(meta[0:2])

            # Molecular weight field
            #
            # positions:
            # 0:2      ntl
            # 3:9      date
            # 10:50    5*(element + amount)
            # 50:52    phase
            # 52:65    molecular weight
            # 65:80    Hf

            molecular_weight = _parse_float(
                meta[52:65]
            )

        except (ValueError, IndexError):
            i += 1
            continue

        ranges = []

        cursor = i + 2

        for _ in range(n_ranges):

            range_line = lines[cursor]

            # NASA format:
            #
            # (2f11.3, i1, 8f5.1, 2x, f15.3)

            T_low = _parse_float(
                range_line[0:11]
            )

            T_high = _parse_float(
                range_line[11:22]
            )

            n_coeff = int(
                range_line[22:23]
            )

            # Read exponents
            exponents = []

            for j in range(8):
                start_pos = 23 + j * 5
                end_pos = start_pos + 5

                field = range_line[
                    start_pos:end_pos
                ].strip()

                if field:
                    exponents.append(
                        _parse_float(field)
                    )

            # --------------------------------
            # Coefficients
            #
            # NASA format:
            #
            # 5D16.8 /
            # 2D16.8, 16X, 2D16.8
            #
            # Therefore ONLY TWO lines.
            # --------------------------------

            line1 = lines[cursor + 1]
            line2 = lines[cursor + 2]

            coeffs = []

            # a1 ... a5
            for j in range(5):

                field = line1[
                    j * 16:(j + 1) * 16
                ]

                coeffs.append(
                    _parse_float(field)
                )

            # a6, a7
            for j in range(2):

                field = line2[
                    j * 16:(j + 1) * 16
                ]

                coeffs.append(
                    _parse_float(field)
                )

            # b1, b2
            #
            # line2:
            #
            # a6 | a7 | 16 spaces | b1 | b2

            b1 = _parse_float(
                line2[48:64]
            )

            b2 = _parse_float(
                line2[64:80]
            )

            ranges.append({
                "T_low": T_low,
                "T_high": T_high,
                "n_coeff": n_coeff,
                "exponents": exponents[:n_coeff],
                "coeffs": coeffs,
                "b1": b1,
                "b2": b2,
            })

            # One range =
            # 1 range line + 2 coefficient lines
            cursor += 3

        if name.lower() == target_species.lower():

            return {
                "name": name,
                "molecular_weight": molecular_weight,
                "ranges": ranges,
            }

        i = cursor

    raise KeyError(
        f"Species '{target_species}' not found"
    )


def cp_species(T: float, species: Dict[str, Any]) -> float:

    if not math.isfinite(T):
        raise ValueError("Temperature must be a finite number")

    interval = None

    for r in species["ranges"]:

        if r["T_low"] <= T <= r["T_high"]:
            interval = r
            break

    if interval is None:

        Tmin = min(
            r["T_low"]
            for r in species["ranges"]
        )

        Tmax = max(
            r["T_high"]
            for r in species["ranges"]
        )

        raise ValueError(
            f"T={T} K outside valid range "
            f"{Tmin} - {Tmax} K"
        )

    # Cp / R
    cp_over_R = 0.0

    for a, exponent in zip(
        interval["coeffs"],
        interval["exponents"]
    ):
        cp_over_R += (
            a * T**exponent
        )

    # kg/kmol
    M = species["molecular_weight"]

    # J/(kg*K)
    R_specific = R_UNIVERSAL / M

    return R_specific * cp_over_R


# ==========================
# Load database ONCE
# ==========================

THERMO_PATH = Path(__file__).with_name("thermo.inp")

AIR = load_species_from_thermo(
    THERMO_PATH,
    "Air"
)


def cp_air(T: float) -> float:
    """Return fixed-composition dry-air Cp in J/(kg*K) at temperature T [K]."""
    return cp_species(T, AIR)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "cp calculator for dry air.  "
            "Temperature unit is K."
        )
    )
    parser.add_argument(
        "temperature",
        type=float,
        nargs="?",
        help="air temperature in K (300~6000 K)",
    )
    args = parser.parse_args()

    temperature = args.temperature
    if temperature is None:
        try:
            temperature = float(input("air temperature [K]: ").strip())
        except ValueError:
            parser.error("temperature must be a number")

    try:
        cp = cp_air(temperature)
    except ValueError as error:
        parser.error(str(error))

    print(f"T  = {temperature:g} K")
    print(f"Cp = {cp:.3f} J/(kg*K)")
    print(f"   = {cp / 1000:.6f} kJ/(kg*K)")

if __name__ == "__main__":
    main()
