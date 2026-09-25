"""Small-strain STATIC solution for a bonded concentric elastic sphere.

Compression is negative stress; external_pressure is a positive magnitude.
This is a 3D reference, not a transient ramp solution or a plane-strain formula.
Python standard library only. All defaults match the pressure PFW input.
"""
from __future__ import annotations
import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Sphere:
    a: float = .22
    b: float = .72
    matrix_E: float = 1.
    inclusion_E: float = 4.
    nu: float = .25
    pressure: float = .004

    def coefficients(self):
        if not all(math.isfinite(v) for v in (self.a,self.b,self.matrix_E,self.inclusion_E,self.nu,self.pressure)):
            raise ValueError("Reference parameters must be finite")
        if not 0 < self.a < self.b or min(self.matrix_E, self.inclusion_E) <= 0:
            raise ValueError("Require 0<a<b and positive elastic moduli")
        if not -1 < self.nu < .5 or not math.isfinite(self.pressure):
            raise ValueError("Invalid Poisson ratio or pressure")
        km = self.matrix_E / (3*(1-2*self.nu))
        ki = self.inclusion_E / (3*(1-2*self.nu))
        gm = self.matrix_E / (2*(1+self.nu))
        # u_i=Ai*r, u_m=Am*r+B/r^2. Continuity at a gives Ai=Am+B/a^3.
        # Radial traction at a and sigma_rr(b)=-p determine Am and B.
        ratio = -3*(ki-km)*self.a**3 / (3*ki+4*gm)  # B/Am
        am = -self.pressure / (3*km - 4*gm*ratio/self.b**3)
        bm = ratio*am
        return am+bm/self.a**3, am, bm

    def evaluate(self, r, side=None):
        if r < 0 or r > self.b*(1+1e-12):
            raise ValueError("Radius outside the physical sphere")
        ai, am, bm = self.coefficients()
        ki = self.inclusion_E/(3*(1-2*self.nu))
        km = self.matrix_E/(3*(1-2*self.nu))
        gm = self.matrix_E/(2*(1+self.nu))
        inclusion = r < self.a if side is None else side == "inclusion"
        if inclusion:
            return dict(region="inclusion", displacement_r=ai*r,
                        stress_rr=3*ki*ai, stress_tt=3*ki*ai)
        if r <= 0:
            raise ValueError("Matrix evaluation requires positive radius")
        return dict(region="matrix", displacement_r=am*r+bm/r**2,
                    stress_rr=3*km*am-4*gm*bm/r**3,
                    stress_tt=3*km*am+2*gm*bm/r**3)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("spherical_static_reference.csv"))
    parser.add_argument("--a", type=float, default=.22)
    parser.add_argument("--b", type=float, default=.72)
    parser.add_argument("--matrix-E", type=float, default=1.)
    parser.add_argument("--inclusion-E", type=float, default=4.)
    parser.add_argument("--nu", type=float, default=.25)
    parser.add_argument("--pressure", type=float, default=.004)
    args=parser.parse_args()
    sphere=Sphere(args.a,args.b,args.matrix_E,args.inclusion_E,args.nu,args.pressure)
    radii=[(sphere.a*i/100,"inclusion") for i in range(101)]
    radii += [(sphere.a+(sphere.b-sphere.a)*i/200,"matrix") for i in range(201)]
    with args.output.open("w", newline="") as stream:
        writer=csv.DictWriter(stream, fieldnames=["radius","region","displacement_r","stress_rr","stress_tt"])
        writer.writeheader()
        for r,side in radii:
            writer.writerow(dict(radius=r,**sphere.evaluate(r,side)))
    print("Wrote static reference (NOT simulated data):",args.output)


if __name__ == "__main__":
    main()
