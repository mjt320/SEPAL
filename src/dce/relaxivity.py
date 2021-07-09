"""Relaxivity models.

Classes: c_to_r_model and derived subclasses:
    c_to_r_linear
"""

from abc import ABC, abstractmethod


class c_to_r_model(ABC):
    """Abstract base class for relaxivity models.

    Subclasses correspond to specific relaxivity models (e.g. linear).
    The main purpose of these classes is to convert tracer concentration to
    relaxation rates.

    Methods
    -------
    R1(R10, c):
        get the R1 relaxation rates and corresponding population fractions for
        each exponential T1 component
    """
    
    @abstractmethod
    def R1(self, R10, c):
        pass
    
    @abstractmethod
    def R2(self, R20, c):
        pass
   


class c_to_r_linear(c_to_r_model):
    
    def __init__(self, r1, r2):
        self.r1 = r1
        self.r2 = r2
    
    def R1(self, R10, c):
        R1 = R10 + self.r1 * c
        return R1
    
    def R2(self, R20, c):
        R2 = R20 + self.r2 * c
        return R2
    
