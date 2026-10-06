"""DSSAT daily fields: one separator plus a five-character numeric token."""
import math


def format_wth_value(value, decimals=1):
    """Keep the separator DSSAT's header reader skips; never clip a value.

    Nonfinite source values become the DSSAT missing marker. Reduce decimal
    precision only when rounding would fill the separator column. Values that
    cannot fit even as integers raise rather than corrupt another column.
    """
    if value is None:
        return "   -99"
    number = float(value)
    if not math.isfinite(number) or number == -99:
        return "   -99"
    for precision in range(int(decimals), -1, -1):
        token = f"{number:.{precision}f}"
        if len(token) <= 5:
            return token.rjust(6)
    raise ValueError(f"DSSAT daily value cannot fit a five-character token: {value!r}")


def wind_run(value):
    """Convert provider m/s to DSSAT km/day, preserving missing wind."""
    if value is None:
        return -99.0
    number = float(value)
    return number * 86.4 if math.isfinite(number) and number >= 0 else -99.0
