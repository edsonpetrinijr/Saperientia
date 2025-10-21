from manimlib.utils.bezier import bezier
import numpy as np

def stars_func(t: float, pull_factor: float = 0.9) -> float:
    return bezier([0, 0, pull_factor, pull_factor, 1, 1, 1])(t)

def fast_to_slow(t: float, power: float = 15) -> float:
    return 1 - (1 - t)**power

def slow_to_fast(t: float, power: float = 0.065) -> float:
    return 1 - (1 - t)**power