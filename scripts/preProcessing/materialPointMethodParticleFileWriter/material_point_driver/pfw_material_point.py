#!/usr/bin/env python3
"""
Single material-point driver front end for GEOS-MPM constitutive testing.

This module is deliberately independent of the SolidMechanicsMPM solver.  It
implements the loading-path controller, material-direction kinematics, state
bookkeeping, and stress-power energy integration needed for single-element / single
particle tests.  The constitutive update is delegated to a backend so that the same
front end can drive a future compiled GEOS constitutive-point executable without
adding branches to the MPM solver.

The built-in elastic backend is intended for driver verification.  Production MPM
material calibration should use a backend that calls the compiled GEOS constitutive
objects through the same update interface as SolidMechanicsMPM.
"""

from __future__ import annotations

import argparse
import copy
import csv
import getpass
import importlib
import importlib.util
import json
import math
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple, Union

import numpy as np


VOIGT_COMPONENTS = ("xx", "yy", "zz", "yz", "xz", "xy")
DEFAULT_COMPILED_DRIVER_NAME = "geos_mpm_material_point_driver"
VOIGT_INDEX = {name: i for i, name in enumerate(VOIGT_COMPONENTS)}
VOIGT_TENSOR_INDEX = {
    "xx": (0, 0),
    "yy": (1, 1),
    "zz": (2, 2),
    "yz": (1, 2),
    "xz": (0, 2),
    "xy": (0, 1),
}


class MaterialPointInputError(ValueError):
    """Raised when a material-point input file is malformed."""


class BackendError(RuntimeError):
    """Raised when a constitutive backend fails."""



def _as_array(value: Any, shape: Optional[Tuple[int, ...]] = None, name: str = "array") -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if shape is not None and arr.shape != shape:
        raise MaterialPointInputError(f"{name} must have shape {shape}, got {arr.shape}")
    return arr


def _normalize(v: Sequence[float], name: str = "vector") -> np.ndarray:
    arr = _as_array(v, (3,), name)
    norm = float(np.linalg.norm(arr))
    if norm <= 1.0e-14:
        raise MaterialPointInputError(f"{name} has near-zero magnitude")
    return arr / norm


def _choose_tangent(normal: np.ndarray, tangent_hint: Optional[Sequence[float]] = None) -> np.ndarray:
    if tangent_hint is not None:
        hint = np.asarray(tangent_hint, dtype=float)
        hint = hint - float(np.dot(hint, normal)) * normal
        if np.linalg.norm(hint) > 1.0e-12:
            return hint / np.linalg.norm(hint)

    candidates = (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0]))
    best = min(candidates, key=lambda c: abs(float(np.dot(c, normal))))
    tangent = best - float(np.dot(best, normal)) * normal
    return tangent / np.linalg.norm(tangent)


def frame_from_normal(normal: Sequence[float], tangent_hint: Optional[Sequence[float]] = None) -> np.ndarray:
    """Build a right-handed, row-wise 3x3 GEOS material frame from one normal/fiber.

    Row 0 is the supplied direction.  Rows 1 and 2 are orthonormal tangents.
    """
    n = _normalize(normal, "materialDirection normal")
    t1 = _choose_tangent(n, tangent_hint)
    t2 = np.cross(n, t1)
    t2 = t2 / np.linalg.norm(t2)
    return np.vstack((n, t1, t2))


def orthonormalize_frame_rows(frame: np.ndarray) -> np.ndarray:
    """Return a right-handed orthonormal frame, preserving row 0 as much as possible."""
    q0 = _normalize(frame[0], "materialDirection row 0")
    q1 = np.asarray(frame[1], dtype=float) - float(np.dot(frame[1], q0)) * q0
    if np.linalg.norm(q1) <= 1.0e-12:
        q1 = _choose_tangent(q0, frame[2])
    else:
        q1 = q1 / np.linalg.norm(q1)
    q2 = np.cross(q0, q1)
    q2 = q2 / np.linalg.norm(q2)
    return np.vstack((q0, q1, q2))


def parse_material_frame(spec: Optional[Mapping[str, Any]]) -> Tuple[np.ndarray, str]:
    """Parse a GEOS/PFW-style material direction specification.

    Accepted forms:
      {"type": "normal", "value": [nx, ny, nz], "tangentHint": [...]}
      {"type": "fiber",  "value": [fx, fy, fz], "tangentHint": [...]}
      {"type": "matrix", "rows": [[...], [...], [...]]}
      {"type": "matrix", "value": [[...], [...], [...]]}

    The returned frame is row-wise, matching GEOS/PFW materialDirection storage.
    """
    if spec is None:
        return np.eye(3), "fixed"

    kind = str(spec.get("type", "matrix")).lower()
    update = str(spec.get("update", "fixed"))
    if kind in ("normal", "fiber", "vector"):
        value = spec.get("value", spec.get("normal", spec.get("fiber")))
        if value is None:
            raise MaterialPointInputError("materialDirection normal/fiber requires a value")
        frame = frame_from_normal(value, spec.get("tangentHint"))
    elif kind in ("matrix", "frame", "basis", "3x3"):
        value = spec.get("rows", spec.get("value"))
        if value is None:
            raise MaterialPointInputError("materialDirection matrix requires rows or value")
        frame = orthonormalize_frame_rows(_as_array(value, (3, 3), "materialDirection matrix"))
    elif kind == "graphite_theta":
        theta = math.radians(float(spec["theta_deg"]))
        normal = [0.0, math.cos(theta), math.sin(theta)]
        frame = frame_from_normal(normal, spec.get("tangentHint", [1.0, 0.0, 0.0]))
        if "update" not in spec:
            update = "graphite"
    else:
        raise MaterialPointInputError(f"Unsupported materialDirection type {kind!r}")
    return frame, update


def voigt_to_tensor(voigt: Sequence[float]) -> np.ndarray:
    v = _as_array(voigt, (6,), "Voigt tensor")
    return np.array(
        [
            [v[0], v[5], v[4]],
            [v[5], v[1], v[3]],
            [v[4], v[3], v[2]],
        ],
        dtype=float,
    )


def tensor_to_voigt(tensor: Sequence[Sequence[float]]) -> np.ndarray:
    t = _as_array(tensor, (3, 3), "tensor")
    return np.array([t[0, 0], t[1, 1], t[2, 2], t[1, 2], t[0, 2], t[0, 1]], dtype=float)


def strain_increment_from_L(L: np.ndarray, dt: float) -> np.ndarray:
    return np.array(
        [
            L[0, 0] * dt,
            L[1, 1] * dt,
            L[2, 2] * dt,
            (L[1, 2] + L[2, 1]) * dt,
            (L[0, 2] + L[2, 0]) * dt,
            (L[0, 1] + L[1, 0]) * dt,
        ],
        dtype=float,
    )


def L_from_strain_increment(dstrain: Sequence[float], dt: float) -> np.ndarray:
    de = _as_array(dstrain, (6,), "strain increment")
    L = np.zeros((3, 3), dtype=float)
    L[0, 0] = de[0] / dt
    L[1, 1] = de[1] / dt
    L[2, 2] = de[2] / dt
    L[1, 2] = L[2, 1] = 0.5 * de[3] / dt
    L[0, 2] = L[2, 0] = 0.5 * de[4] / dt
    L[0, 1] = L[1, 0] = 0.5 * de[5] / dt
    return L


def cofactor(F: np.ndarray) -> np.ndarray:
    return float(np.linalg.det(F)) * np.linalg.inv(F).T


def polar_rotation(F: np.ndarray) -> np.ndarray:
    U, _, Vt = np.linalg.svd(F)
    R = U @ Vt
    if np.linalg.det(R) < 0.0:
        U[:, -1] *= -1.0
        R = U @ Vt
    return R


def update_material_frame(reference_frame: np.ndarray, F: np.ndarray, R: np.ndarray, mode: str) -> np.ndarray:
    """Update row-wise material frame using GEOS-MPM-friendly kinematics."""
    mode_key = mode.lower()
    Q0 = reference_frame

    if mode_key in ("fixed", "none", "reference"):
        return Q0.copy()

    if mode_key in ("rotation", "rotational"):
        return orthonormalize_frame_rows(Q0 @ R.T)

    if mode_key in ("fiber", "fiberlike", "fiber_like"):
        Q = Q0 @ F.T
        return np.vstack([_normalize(Q[i], f"fiber materialDirection row {i}") for i in range(3)])

    cofF = cofactor(F)

    if mode_key in ("mpmcofactor", "mpm_cofactor", "cofactor"):
        Q = Q0 @ cofF.T
        return np.vstack([_normalize(Q[i], f"cofactor materialDirection row {i}") for i in range(3)])

    if mode_key in ("normal", "normal_like", "graphite", "graphitelike", "graphite_like"):
        normal = _normalize(Q0[0] @ cofF.T, "normal materialDirection row 0")
        tangent = Q0[1] @ F.T
        tangent = tangent - float(np.dot(tangent, normal)) * normal
        if np.linalg.norm(tangent) <= 1.0e-12:
            tangent = Q0[2] @ F.T
            tangent = tangent - float(np.dot(tangent, normal)) * normal
        if np.linalg.norm(tangent) <= 1.0e-12:
            tangent = _choose_tangent(normal)
        else:
            tangent = tangent / np.linalg.norm(tangent)
        binormal = np.cross(normal, tangent)
        binormal = binormal / np.linalg.norm(binormal)
        return np.vstack((normal, tangent, binormal))

    raise MaterialPointInputError(f"Unsupported materialDirection update mode {mode!r}")


def stress_power(stress_voigt: Sequence[float], L: np.ndarray) -> float:
    """Return sigma:L using GEOS Voigt order and a full velocity-gradient tensor."""
    return float(np.tensordot(voigt_to_tensor(stress_voigt), L, axes=2))


def evaluate_schedule(value: Any, time: float, default: Optional[float] = None) -> float:
    """Evaluate a scalar or piecewise-linear schedule.

    A schedule may be a scalar, {"value": scalar}, or {"table": [[t0,v0], ...]}.
    """
    if value is None:
        if default is None:
            raise MaterialPointInputError("Missing required scalar value/schedule")
        return float(default)
    if isinstance(value, Mapping):
        if "value" in value:
            return float(value["value"])
        if "target" in value:
            return evaluate_schedule(value["target"], time, default)
        if "table" in value:
            return evaluate_schedule(value["table"], time, default)
        raise MaterialPointInputError(f"Unsupported schedule dictionary keys: {list(value.keys())}")
    if isinstance(value, (list, tuple)) and value and isinstance(value[0], (list, tuple)):
        table = sorted((float(row[0]), float(row[1])) for row in value)
        if time <= table[0][0]:
            return table[0][1]
        if time >= table[-1][0]:
            return table[-1][1]
        for (t0, v0), (t1, v1) in zip(table[:-1], table[1:]):
            if t0 <= time <= t1:
                alpha = (time - t0) / (t1 - t0)
                return (1.0 - alpha) * v0 + alpha * v1
        return table[-1][1]
    return float(value)


@dataclass
class MaterialPointState:
    stress: np.ndarray = field(default_factory=lambda: np.zeros(6))
    total_strain: np.ndarray = field(default_factory=lambda: np.zeros(6))
    F: np.ndarray = field(default_factory=lambda: np.eye(3))
    material_frame_reference: np.ndarray = field(default_factory=lambda: np.eye(3))
    material_frame: np.ndarray = field(default_factory=lambda: np.eye(3))
    temperature: float = 300.0
    temperature_rate: float = 0.0
    density0: float = 1.0
    density: float = 1.0
    specific_internal_energy: float = 0.0
    accumulated_stress_power: float = 0.0
    length_scale: float = 1.0
    strength_scale: float = 1.0
    jacobian: float = 1.0
    extra: Dict[str, Any] = field(default_factory=dict)

    def clone(self) -> "MaterialPointState":
        return copy.deepcopy(self)


@dataclass
class StepResult:
    step: int
    time: float
    dt: float
    converged: bool
    iterations: int
    residual_norm: float
    L: np.ndarray
    strain_increment: np.ndarray
    stress_power: float
    state: MaterialPointState


@dataclass
class CompiledDriverRunResult:
    """Sidecar files and command line for the optional compiled GEOS backend."""

    case: Path
    material_xml: Path
    load_path_csv: Path
    run_script: Path
    output_csv: Path
    command: List[str]
    log_file: Optional[Path] = None
    returncode: Optional[int] = None


class ConstitutiveBackend:
    """Abstract material update backend."""

    name = "abstract"

    def update(
        self,
        state: MaterialPointState,
        dt: float,
        L: np.ndarray,
        F_new: np.ndarray,
        R_old: np.ndarray,
        R_new: np.ndarray,
        material_frame_new: np.ndarray,
    ) -> MaterialPointState:
        raise NotImplementedError


class LinearElasticBackend(ConstitutiveBackend):
    """Small-strain elastic backend for driver verification.

    Supported material dictionaries:
      - ElasticIsotropic and derived dictionaries with K/G, E/nu, E/G, K/nu, or lambda/G.
      - ElasticTransverseIsotropic and Graphite dictionaries via engineering constants.

    This backend is not a replacement for the compiled GEOS Graphite inelastic model; it is
    included so mixed control, material-frame update, and energy integration can be tested
    without a GEOS build.
    """

    name = "elastic"

    def __init__(self, material: Mapping[str, Any]):
        self.material = material
        self.model = str(material.get("model", "ElasticIsotropic"))
        if self.model in ("ElasticTransverseIsotropic", "Graphite"):
            self.kind = "transverse_isotropic"
            self.C_local = self._build_ti_stiffness(material)
        else:
            self.kind = "isotropic"
            K, G = self._isotropic_constants(material)
            lam = K - 2.0 * G / 3.0
            C = np.zeros((6, 6), dtype=float)
            C[:3, :3] = lam
            for i in range(3):
                C[i, i] += 2.0 * G
            C[3, 3] = G
            C[4, 4] = G
            C[5, 5] = G
            self.C = C

    @staticmethod
    def _isotropic_constants(material: Mapping[str, Any]) -> Tuple[float, float]:
        m = material
        if "defaultBulkModulus" in m and "defaultShearModulus" in m:
            return float(m["defaultBulkModulus"]), float(m["defaultShearModulus"])
        if "defaultYoungModulus" in m and "defaultPoissonRatio" in m:
            E = float(m["defaultYoungModulus"])
            nu = float(m["defaultPoissonRatio"])
            K = E / (3.0 * (1.0 - 2.0 * nu))
            G = E / (2.0 * (1.0 + nu))
            return K, G
        if "defaultYoungModulus" in m and "defaultShearModulus" in m:
            E = float(m["defaultYoungModulus"])
            G = float(m["defaultShearModulus"])
            nu = E / (2.0 * G) - 1.0
            K = E / (3.0 * (1.0 - 2.0 * nu))
            return K, G
        if "defaultBulkModulus" in m and "defaultPoissonRatio" in m:
            K = float(m["defaultBulkModulus"])
            nu = float(m["defaultPoissonRatio"])
            G = 3.0 * K * (1.0 - 2.0 * nu) / (2.0 * (1.0 + nu))
            return K, G
        if "defaultLambda" in m and "defaultShearModulus" in m:
            lam = float(m["defaultLambda"])
            G = float(m["defaultShearModulus"])
            K = lam + 2.0 * G / 3.0
            return K, G
        raise MaterialPointInputError("Elastic backend could not infer isotropic constants from material")

    @staticmethod
    def _build_ti_stiffness(material: Mapping[str, Any]) -> np.ndarray:
        # Local axis 0 is the material direction row 0.  Rows 1 and 2 are the transverse plane.
        if all(k in material for k in ("defaultC11", "defaultC13", "defaultC33", "defaultC44", "defaultC66")):
            # Input convention in GEOS documentation is usually axis-2/3 as the symmetry axis.
            # Re-map to a local axis-0 convention by constructing compliance from engineering
            # constants is not possible from stiffness constants without more convention metadata.
            # Use the common TI stiffness with axis 0 here: C00=C33, C11=C22=C11, C01=C02=C13,
            # C12=C11-2*C66, C02=C13, shear involving axis0=C44, transverse shear=C66.
            c11 = float(material["defaultC11"])
            c13 = float(material["defaultC13"])
            c33 = float(material["defaultC33"])
            c44 = float(material["defaultC44"])
            c66 = float(material["defaultC66"])
            C = np.zeros((6, 6), dtype=float)
            C[0, 0] = c33
            C[1, 1] = C[2, 2] = c11
            C[0, 1] = C[1, 0] = c13
            C[0, 2] = C[2, 0] = c13
            C[1, 2] = C[2, 1] = c11 - 2.0 * c66
            C[3, 3] = c66  # yz transverse-plane shear
            C[4, 4] = c44  # xz axial-transverse shear
            C[5, 5] = c44  # xy axial-transverse shear
            return C

        required = (
            "defaultYoungModulusTransverse",
            "defaultYoungModulusAxial",
            "defaultPoissonRatioTransverse",
            "defaultPoissonRatioAxialTransverse",
            "defaultShearModulusAxialTransverse",
        )
        missing = [k for k in required if k not in material]
        if missing:
            raise MaterialPointInputError(f"TI elastic backend missing constants: {missing}")

        Et = float(material["defaultYoungModulusTransverse"])
        Ea = float(material["defaultYoungModulusAxial"])
        nut = float(material["defaultPoissonRatioTransverse"])
        nuat = float(material["defaultPoissonRatioAxialTransverse"])
        Gat = float(material["defaultShearModulusAxialTransverse"])
        Gt = Et / (2.0 * (1.0 + nut))

        # Compliance in local GEOS order with symmetry axis = component 0.
        S = np.zeros((6, 6), dtype=float)
        S[0, 0] = 1.0 / Ea
        S[1, 1] = S[2, 2] = 1.0 / Et
        S[1, 2] = S[2, 1] = -nut / Et
        S[0, 1] = S[1, 0] = -nuat / Ea
        S[0, 2] = S[2, 0] = -nuat / Ea
        S[3, 3] = 1.0 / Gt
        S[4, 4] = 1.0 / Gat
        S[5, 5] = 1.0 / Gat
        return np.linalg.inv(S)

    def _ti_increment(self, material_frame: np.ndarray, strain_increment_voigt: np.ndarray) -> np.ndarray:
        E_lab = voigt_to_tensor(strain_increment_voigt)
        Q = material_frame
        E_local = Q @ E_lab @ Q.T
        de_local = tensor_to_voigt(E_local)
        ds_local = self.C_local @ de_local
        S_local = voigt_to_tensor(ds_local)
        S_lab = Q.T @ S_local @ Q
        return tensor_to_voigt(S_lab)

    def update(
        self,
        state: MaterialPointState,
        dt: float,
        L: np.ndarray,
        F_new: np.ndarray,
        R_old: np.ndarray,
        R_new: np.ndarray,
        material_frame_new: np.ndarray,
    ) -> MaterialPointState:
        del R_old, R_new
        new_state = state.clone()
        de = strain_increment_from_L(L, dt)
        if self.kind == "isotropic":
            ds = self.C @ de
        else:
            ds = self._ti_increment(material_frame_new, de)
        new_state.stress = state.stress + ds
        new_state.total_strain = state.total_strain + de
        new_state.F = F_new.copy()
        new_state.material_frame = material_frame_new.copy()
        new_state.jacobian = float(np.linalg.det(F_new))
        new_state.density = new_state.density0 / new_state.jacobian
        return new_state


class ExternalJsonBackend(ConstitutiveBackend):
    """Backend adapter for a compiled GEOS material-point executable.

    The external command receives one JSON document on stdin and must return one JSON
    document on stdout.  This keeps the Python control/optimization layer independent
    from the compiled constitutive implementation.
    """

    name = "external-json"

    def __init__(self, command: Sequence[str], material: Mapping[str, Any]):
        if not command:
            raise MaterialPointInputError("external-json backend requires a command")
        self.command = list(command)
        self.material = dict(material)

    def update(
        self,
        state: MaterialPointState,
        dt: float,
        L: np.ndarray,
        F_new: np.ndarray,
        R_old: np.ndarray,
        R_new: np.ndarray,
        material_frame_new: np.ndarray,
    ) -> MaterialPointState:
        request = {
            "material": self.material,
            "dt": dt,
            "L": L.tolist(),
            "FNew": F_new.tolist(),
            "RBeginning": R_old.tolist(),
            "REnd": R_new.tolist(),
            "materialDirection": material_frame_new.tolist(),
            "state": state_to_json_dict(state),
        }
        try:
            proc = subprocess.run(
                self.command,
                input=json.dumps(request),
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
        except OSError as exc:
            raise BackendError(f"Failed to launch external backend {self.command!r}: {exc}") from exc
        if proc.returncode != 0:
            raise BackendError(
                f"External backend returned {proc.returncode}. stderr:\n{proc.stderr}"
            )
        try:
            response = json.loads(proc.stdout)
        except json.JSONDecodeError as exc:
            raise BackendError(f"External backend did not return valid JSON. stdout:\n{proc.stdout}") from exc
        return state_from_json_dict(response.get("state", response))


def state_to_json_dict(state: MaterialPointState) -> Dict[str, Any]:
    return {
        "stress": state.stress.tolist(),
        "totalStrain": state.total_strain.tolist(),
        "F": state.F.tolist(),
        "materialFrameReference": state.material_frame_reference.tolist(),
        "materialFrame": state.material_frame.tolist(),
        "temperature": state.temperature,
        "temperatureRate": state.temperature_rate,
        "density0": state.density0,
        "density": state.density,
        "specificInternalEnergy": state.specific_internal_energy,
        "accumulatedStressPower": state.accumulated_stress_power,
        "lengthScale": state.length_scale,
        "strengthScale": state.strength_scale,
        "jacobian": state.jacobian,
        "extra": copy.deepcopy(state.extra),
    }


def state_from_json_dict(data: Mapping[str, Any]) -> MaterialPointState:
    state = MaterialPointState()
    state.stress = _as_array(data.get("stress", state.stress), (6,), "state.stress")
    state.total_strain = _as_array(data.get("totalStrain", data.get("total_strain", state.total_strain)), (6,), "state.totalStrain")
    state.F = _as_array(data.get("F", state.F), (3, 3), "state.F")
    state.material_frame_reference = _as_array(
        data.get("materialFrameReference", data.get("material_frame_reference", state.material_frame_reference)),
        (3, 3),
        "state.materialFrameReference",
    )
    state.material_frame = _as_array(data.get("materialFrame", data.get("material_frame", state.material_frame)), (3, 3), "state.materialFrame")
    state.temperature = float(data.get("temperature", state.temperature))
    state.temperature_rate = float(data.get("temperatureRate", data.get("temperature_rate", state.temperature_rate)))
    state.density0 = float(data.get("density0", state.density0))
    state.density = float(data.get("density", state.density))
    state.specific_internal_energy = float(data.get("specificInternalEnergy", data.get("specific_internal_energy", state.specific_internal_energy)))
    state.accumulated_stress_power = float(data.get("accumulatedStressPower", data.get("accumulated_stress_power", state.accumulated_stress_power)))
    state.length_scale = float(data.get("lengthScale", data.get("length_scale", state.length_scale)))
    state.strength_scale = float(data.get("strengthScale", data.get("strength_scale", state.strength_scale)))
    state.jacobian = float(data.get("jacobian", state.jacobian))
    state.extra = copy.deepcopy(data.get("extra", {}))
    return state


def load_material(material_spec: Mapping[str, Any]) -> Dict[str, Any]:
    """Load a PFW material dictionary and apply optional parameter overrides.

    A case may provide either ``{"source": "pfw_materials.py", "name": ...}``,
    ``{"inline": {...}}``, or a direct material dictionary containing ``model``.
    The direct form is useful for optimizer-generated cases that have already
    copied and modified parameters from ``pfw_materials.py``.
    """
    if "inline" in material_spec:
        material = copy.deepcopy(material_spec["inline"])
    elif "model" in material_spec and "source" not in material_spec and "module" not in material_spec:
        material = copy.deepcopy({k: v for k, v in material_spec.items() if k != "overrides"})
    else:
        module_name = str(material_spec.get("module", material_spec.get("source", "pfw_materials")))
        if module_name.endswith(".py"):
            module_name = Path(module_name).stem
        material_name = material_spec.get("name")
        if not material_name:
            raise MaterialPointInputError("material block requires name unless inline material is supplied")
        pfw_root = Path(__file__).resolve().parents[1]
        if str(pfw_root) not in sys.path:
            sys.path.insert(0, str(pfw_root))
        module = importlib.import_module(module_name)
        try:
            material_obj = getattr(module, str(material_name))
        except AttributeError as exc:
            raise MaterialPointInputError(f"Could not find material {material_name!r} in module {module_name!r}") from exc
        material = material_obj.copy() if hasattr(material_obj, "copy") else copy.deepcopy(material_obj)

    overrides = material_spec.get("overrides", {})
    for key, value in overrides.items():
        material[key] = value
    if hasattr(material, "refreshMaterialString"):
        material.refreshMaterialString()
    return dict(material)


def make_backend(case: Mapping[str, Any], material: Mapping[str, Any]) -> ConstitutiveBackend:
    backend_spec = case.get("backend", "elastic")
    if isinstance(backend_spec, str):
        backend_name = backend_spec
        backend_options: Mapping[str, Any] = {}
    elif isinstance(backend_spec, Mapping):
        backend_name = str(backend_spec.get("type", "elastic"))
        backend_options = backend_spec
    else:
        raise MaterialPointInputError("backend must be a string or dictionary")

    key = backend_name.lower()
    if key in ("elastic", "linear-elastic", "linear_elastic"):
        return LinearElasticBackend(material)
    if key in ("external-json", "external_json", "geos-json", "geos_json"):
        command = backend_options.get("command")
        if isinstance(command, str):
            command = command.split()
        return ExternalJsonBackend(command, material)
    raise MaterialPointInputError(f"Unsupported backend type {backend_name!r}")


def initialize_state(case: Mapping[str, Any], material: Mapping[str, Any]) -> Tuple[MaterialPointState, str]:
    initial = case.get("initial", {})
    frame, update_mode = parse_material_frame(case.get("materialDirection", initial.get("materialDirection")))
    state = MaterialPointState()
    state.stress = _as_array(initial.get("stress", np.zeros(6)), (6,), "initial.stress")
    state.total_strain = _as_array(initial.get("strain", initial.get("totalStrain", np.zeros(6))), (6,), "initial.strain")
    state.F = _as_array(initial.get("F", np.eye(3)), (3, 3), "initial.F")
    state.material_frame_reference = frame.copy()
    state.material_frame = update_material_frame(frame, state.F, polar_rotation(state.F), update_mode)
    state.temperature = float(initial.get("temperature", case.get("temperature", {}).get("initial", 300.0)))
    state.temperature_rate = float(initial.get("temperatureRate", 0.0))
    density0 = float(initial.get("density0", initial.get("density", material.get("defaultDensity", 1.0))))
    state.density0 = density0
    state.jacobian = float(np.linalg.det(state.F))
    state.density = float(initial.get("density", density0 / state.jacobian))
    state.specific_internal_energy = float(initial.get("specificInternalEnergy", initial.get("internalEnergy", 0.0)))
    state.accumulated_stress_power = float(initial.get("accumulatedStressPower", 0.0))
    state.length_scale = float(initial.get("lengthScale", initial.get("length_scale", 1.0)))
    state.strength_scale = float(initial.get("strengthScale", 1.0))
    state.extra = copy.deepcopy(initial.get("extra", {}))

    # Preserve model-specific initial state fields without needing to know every wrapper name.
    for key, value in initial.items():
        if key not in {
            "stress", "strain", "totalStrain", "F", "materialDirection", "temperature", "temperatureRate",
            "density0", "density", "specificInternalEnergy", "internalEnergy", "accumulatedStressPower",
            "lengthScale", "length_scale", "strengthScale", "extra",
        }:
            state.extra[key] = copy.deepcopy(value)
    return state, update_mode


class ControlProgram:
    def __init__(self, controls: Sequence[Mapping[str, Any]]):
        self.controls = list(controls)
        seen = set()
        for ctrl in self.controls:
            comp = str(ctrl.get("component", "")).lower()
            if comp not in VOIGT_INDEX:
                raise MaterialPointInputError(f"Unsupported control component {comp!r}")
            if comp in seen:
                raise MaterialPointInputError(f"Duplicate control component {comp!r}")
            seen.add(comp)
        missing = [c for c in VOIGT_COMPONENTS if c not in seen]
        if missing:
            # Default missing components to zero strain increments.  This makes simple uniaxial
            # inputs terse while keeping the resulting control mask explicit in the output.
            for comp in missing:
                self.controls.append({"component": comp, "mode": "strain", "value": 0.0})
        self.controls.sort(key=lambda c: VOIGT_INDEX[str(c["component"]).lower()])

    @property
    def stress_controls(self) -> List[Mapping[str, Any]]:
        return [c for c in self.controls if str(c.get("mode", "")).lower() in ("stress", "stresscontrol", "stress_control")]

    def prescribed_strain_increment(self, state: MaterialPointState, dt: float, time_new: float) -> np.ndarray:
        del state
        de = np.zeros(6, dtype=float)
        for ctrl in self.controls:
            comp = str(ctrl["component"]).lower()
            idx = VOIGT_INDEX[comp]
            mode = str(ctrl.get("mode", "strain")).lower()
            if mode in ("stress", "stresscontrol", "stress_control"):
                continue
            if mode in ("strain", "strainincrement", "strain_increment"):
                de[idx] = evaluate_schedule(ctrl.get("value", ctrl.get("increment", 0.0)), time_new)
            elif mode in ("strainrate", "strain_rate"):
                de[idx] = evaluate_schedule(ctrl.get("rate", ctrl.get("value")), time_new) * dt
            elif mode in ("truestrainrate", "true_strain_rate", "logstrainrate", "log_strain_rate"):
                de[idx] = evaluate_schedule(ctrl.get("rate", ctrl.get("value")), time_new) * dt
            elif mode in ("l", "velocitygradient", "velocity_gradient"):
                de[idx] = evaluate_schedule(ctrl.get("value", ctrl.get("rate")), time_new) * dt
            else:
                raise MaterialPointInputError(f"Unsupported control mode {mode!r} for component {comp}")
        return de

    def stress_target_vector(self, time_new: float) -> Tuple[List[int], np.ndarray]:
        indices: List[int] = []
        targets: List[float] = []
        for ctrl in self.stress_controls:
            comp = str(ctrl["component"]).lower()
            indices.append(VOIGT_INDEX[comp])
            targets.append(evaluate_schedule(ctrl.get("target", ctrl.get("value")), time_new))
        return indices, np.asarray(targets, dtype=float)

    def describe_modes(self) -> Dict[str, str]:
        return {str(c["component"]).lower(): str(c.get("mode", "strain")) for c in self.controls}


class MaterialPointDriver:
    def __init__(self, case: Mapping[str, Any]):
        self.case = copy.deepcopy(case)
        self.material = load_material(self.case.get("material", {}))
        self.backend = make_backend(self.case, self.material)
        self.state, self.material_direction_update = initialize_state(self.case, self.material)
        self.controls = ControlProgram(self.case.get("control", []))
        time_block = self.case.get("time", self.case.get("steps", {}))
        self.dt = float(time_block.get("dt", self.case.get("dt", 1.0)))
        self.n_steps = int(time_block.get("nSteps", time_block.get("nsteps", self.case.get("nSteps", 1))))
        self.time = float(time_block.get("initial", 0.0))
        self.max_newton_iterations = int(self.case.get("solver", {}).get("maxIterations", 20))
        self.tolerance = float(self.case.get("solver", {}).get("stressTolerance", 1.0e-8))
        self.finite_difference_step = float(self.case.get("solver", {}).get("finiteDifferenceStrain", 1.0e-7))
        self.energy = self.case.get("energy", {})
        self.temperature = self.case.get("temperature", {})

    def trial_update(self, base_state: MaterialPointState, de: np.ndarray, dt: float) -> Tuple[MaterialPointState, np.ndarray]:
        L = L_from_strain_increment(de, dt)
        # Match explicit MPM-style kinematics: Fdot = L F and forward Euler in time.
        F_old = base_state.F
        F_new = (np.eye(3) + L * dt) @ F_old
        R_old = polar_rotation(F_old)
        R_new = polar_rotation(F_new)
        frame_new = update_material_frame(base_state.material_frame_reference, F_new, R_new, self.material_direction_update)
        new_state = self.backend.update(base_state, dt, L, F_new, R_old, R_new, frame_new)
        return new_state, L

    def solve_step(self, step: int) -> StepResult:
        time_new = self.time + self.dt
        base = self.state.clone()
        de_prescribed = self.controls.prescribed_strain_increment(base, self.dt, time_new)
        stress_indices, stress_targets = self.controls.stress_target_vector(time_new)
        unknown_components = stress_indices
        x = np.zeros(len(unknown_components), dtype=float)

        def residual_and_state(x_trial: np.ndarray) -> Tuple[np.ndarray, MaterialPointState, np.ndarray, np.ndarray]:
            de = de_prescribed.copy()
            for val, comp_index in zip(x_trial, unknown_components):
                de[comp_index] = val
            trial_state, L = self.trial_update(base, de, self.dt)
            if not stress_indices:
                residual = np.zeros(0, dtype=float)
            else:
                residual = trial_state.stress[stress_indices] - stress_targets
            return residual, trial_state, L, de

        converged = False
        residual_norm = 0.0
        iterations = 0
        if not unknown_components:
            residual, accepted_state, L, de = residual_and_state(x)
            converged = True
        else:
            accepted_state = base
            L = np.zeros((3, 3), dtype=float)
            de = de_prescribed.copy()
            for iteration in range(1, self.max_newton_iterations + 1):
                iterations = iteration
                residual, trial_state, trial_L, trial_de = residual_and_state(x)
                residual_norm = float(np.linalg.norm(residual, ord=2))
                if residual_norm <= self.tolerance:
                    accepted_state = trial_state
                    L = trial_L
                    de = trial_de
                    converged = True
                    break

                J = np.zeros((len(residual), len(x)), dtype=float)
                for j in range(len(x)):
                    h = self.finite_difference_step * max(1.0, abs(x[j]))
                    x_fd = x.copy()
                    x_fd[j] += h
                    residual_fd, _, _, _ = residual_and_state(x_fd)
                    J[:, j] = (residual_fd - residual) / h

                try:
                    dx = np.linalg.solve(J, -residual)
                except np.linalg.LinAlgError:
                    dx = np.linalg.lstsq(J, -residual, rcond=None)[0]

                # Damped line search.  The backend state is restored by residual_and_state
                # because every trial starts from the same base clone.
                alpha = 1.0
                best_x = x + dx
                best_norm = math.inf
                for _ in range(10):
                    candidate_x = x + alpha * dx
                    candidate_residual, _, _, _ = residual_and_state(candidate_x)
                    candidate_norm = float(np.linalg.norm(candidate_residual, ord=2))
                    if candidate_norm < best_norm:
                        best_norm = candidate_norm
                        best_x = candidate_x
                    if candidate_norm < residual_norm:
                        break
                    alpha *= 0.5
                x = best_x
            if not converged:
                residual, accepted_state, L, de = residual_and_state(x)
                residual_norm = float(np.linalg.norm(residual, ord=2))

        self._integrate_energy_and_temperature(base, accepted_state, L, self.dt)
        self.state = accepted_state
        self.time = time_new
        return StepResult(
            step=step,
            time=self.time,
            dt=self.dt,
            converged=converged,
            iterations=iterations,
            residual_norm=residual_norm,
            L=L,
            strain_increment=de,
            stress_power=stress_power(accepted_state.stress, L),
            state=accepted_state.clone(),
        )

    def _integrate_energy_and_temperature(self, old_state: MaterialPointState, new_state: MaterialPointState, L: np.ndarray, dt: float) -> None:
        mode = str(self.energy.get("mode", "stressPower")).lower()
        if mode in ("off", "none"):
            return
        if mode not in ("stresspower", "stress_power", "driver", "material"):
            raise MaterialPointInputError(f"Unsupported energy mode {self.energy.get('mode')!r}")
        old_power = stress_power(old_state.stress, L)
        new_power = stress_power(new_state.stress, L)
        rho = float(new_state.density if self.energy.get("densityEvaluation", "current") == "current" else old_state.density)
        if rho <= 0.0:
            raise BackendError("Cannot integrate energy with non-positive density")
        delta_work = 0.5 * dt * (old_power + new_power) / rho
        retention = float(self.energy.get("retentionFactor", 1.0))
        new_state.accumulated_stress_power = old_state.accumulated_stress_power + delta_work
        new_state.specific_internal_energy = old_state.specific_internal_energy + retention * delta_work

        temp_mode = str(self.temperature.get("mode", "prescribed")).lower()
        if temp_mode in ("prescribed", "isothermal"):
            target = self.temperature.get("value", self.temperature.get("initial", old_state.temperature))
            new_state.temperature = evaluate_schedule(target, self.time + dt, old_state.temperature)
            new_state.temperature_rate = (new_state.temperature - old_state.temperature) / dt
        elif temp_mode in ("adiabatic", "coupledtemperature", "coupled_temperature"):
            cv = float(self.temperature.get("heatCapacity", self.energy.get("heatCapacity", 1.0)))
            if cv <= 0.0:
                raise MaterialPointInputError("Adiabatic temperature update requires positive heatCapacity")
            new_state.temperature = old_state.temperature + retention * delta_work / cv
            new_state.temperature_rate = (new_state.temperature - old_state.temperature) / dt
        elif temp_mode in ("frommaterial", "from_material"):
            pass
        else:
            raise MaterialPointInputError(f"Unsupported temperature mode {self.temperature.get('mode')!r}")

    def run(self) -> List[StepResult]:
        results = []
        for step in range(1, self.n_steps + 1):
            result = self.solve_step(step)
            results.append(result)
            if not result.converged and self.case.get("solver", {}).get("failOnNonconvergence", True):
                raise BackendError(
                    f"Material-point stress control failed to converge at step {step}: "
                    f"residual norm {result.residual_norm}"
                )
        return results


def flatten_result(result: StepResult, controls: Optional[ControlProgram] = None) -> Dict[str, Any]:
    s = result.state
    row: Dict[str, Any] = {
        "step": result.step,
        "time": result.time,
        "dt": result.dt,
        "converged": int(result.converged),
        "newtonIterations": result.iterations,
        "stressResidualNorm": result.residual_norm,
        "specificInternalEnergy": s.specific_internal_energy,
        "accumulatedStressPower": s.accumulated_stress_power,
        "stressPower": result.stress_power,
        "temperature": s.temperature,
        "temperatureRate": s.temperature_rate,
        "density": s.density,
        "jacobian": s.jacobian,
        "lengthScale": s.length_scale,
        "strengthScale": s.strength_scale,
    }
    if controls is not None:
        for comp, mode in controls.describe_modes().items():
            row[f"control_{comp}"] = mode
    for i, comp in enumerate(VOIGT_COMPONENTS):
        row[f"sigma_{comp}"] = s.stress[i]
        row[f"eps_{comp}"] = s.total_strain[i]
        row[f"deps_{comp}"] = result.strain_increment[i]
    for i in range(3):
        for j in range(3):
            row[f"F_{i}{j}"] = s.F[i, j]
            row[f"L_{i}{j}"] = result.L[i, j]
            row[f"materialFrame_{i}{j}"] = s.material_frame[i, j]
            row[f"materialFrameReference_{i}{j}"] = s.material_frame_reference[i, j]
    for key, value in sorted(s.extra.items()):
        if isinstance(value, (int, float, np.integer, np.floating, str)):
            row[str(key)] = value
    return row


def write_csv(results: Sequence[StepResult], path: Path, controls: Optional[ControlProgram] = None) -> None:
    rows = [flatten_result(r, controls) for r in results]
    if not rows:
        return
    fieldnames: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)



def _compiled_control_mode(mode: str) -> str:
    key = mode.lower()
    if key in ("stress", "stresscontrol", "stress_control"):
        return "stress"
    if key in ("strain", "strainincrement", "strain_increment", "increment"):
        return "strain"
    if key in ("strainrate", "strain_rate", "engineeringstrainrate", "engineering_strain_rate"):
        return "strainRate"
    if key in ("truestrainrate", "true_strain_rate", "logstrainrate", "log_strain_rate"):
        return "trueStrainRate"
    raise MaterialPointInputError(
        f"Control mode {mode!r} cannot be exported to the compiled GEOS driver; "
        "supported modes are strain, strainRate, trueStrainRate, and stress."
    )


def _control_value_for_compiled_driver(ctrl: Mapping[str, Any], time_new: float) -> float:
    mode = str(ctrl.get("mode", "strain"))
    key = mode.lower()
    if key in ("stress", "stresscontrol", "stress_control"):
        return evaluate_schedule(ctrl.get("target", ctrl.get("value")), time_new)
    if key in ("strain", "strainincrement", "strain_increment", "increment"):
        return evaluate_schedule(ctrl.get("value", ctrl.get("increment", 0.0)), time_new)
    if key in ("strainrate", "strain_rate", "engineeringstrainrate", "engineering_strain_rate"):
        return evaluate_schedule(ctrl.get("rate", ctrl.get("value")), time_new)
    if key in ("truestrainrate", "true_strain_rate", "logstrainrate", "log_strain_rate"):
        return evaluate_schedule(ctrl.get("rate", ctrl.get("value")), time_new)
    raise MaterialPointInputError(f"Unsupported compiled-driver control mode {mode!r}")


def _material_xml_string(material: Mapping[str, Any]) -> str:
    xml = material.get("materialString")
    if not xml:
        pfw_root = Path(__file__).resolve().parents[1]
        if str(pfw_root) not in sys.path:
            sys.path.insert(0, str(pfw_root))
        pfw_materials = importlib.import_module("pfw_materials")
        xml = pfw_materials.generateMaterialString(dict(material))
    text = str(xml).strip()
    if text.startswith("<Constitutive"):
        return text + "\n"
    return "<Constitutive>\n" + text + "\n</Constitutive>\n"


def _join_reals(values: Iterable[float]) -> str:
    return ",".join(f"{float(v):.17g}" for v in values)


def _quote_command(command: Sequence[str]) -> str:
    return " \\\n  ".join(shlex.quote(str(part)) for part in command)


def _is_unset_userdefs_value(value: Any) -> bool:
    if value is None:
        return True
    text = str(value).strip()
    return text == "" or text == "CHANGEME" or text == "_CHANGEME"


def _expand_executable_path(value: Any) -> str:
    return str(Path(os.path.expandvars(os.path.expanduser(str(value)))))


def _load_userdefs_from(path: Path) -> Optional[Any]:
    if not path.is_file():
        return None
    module_name = f"_geos_mpm_material_point_userdefs_{abs(hash(str(path)))}"
    try:
        spec = importlib.util.spec_from_file_location(module_name, path)
        if spec is None or spec.loader is None:
            return None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    except Exception:
        return None


def _iter_userdefs_candidates() -> Iterable[Path]:
    """Yield likely userDefs files without requiring the caller's CWD.

    setupMPM writes userDefs_<user>.py next to the PFW scripts, which are often
    reached through the top-level mpm symlink.  Some wrappers copy userDefs into
    the current run directory, so check both locations before falling back to any
    userDefs_*.py found in the PFW root.
    """
    yielded = set()

    explicit = os.environ.get("GEOS_MPM_USERDEFS") or os.environ.get("PFW_USERDEFS")
    if explicit:
        path = Path(os.path.expandvars(os.path.expanduser(explicit))).resolve()
        yielded.add(path)
        yield path

    user_names = []
    for name in (os.environ.get("USER"), os.environ.get("LOGNAME")):
        if name and name not in user_names:
            user_names.append(name)
    try:
        name = getpass.getuser()
    except Exception:
        name = ""
    if name and name not in user_names:
        user_names.append(name)

    pfw_root = Path(__file__).resolve().parents[1]
    bases = [Path.cwd(), pfw_root]
    for base in bases:
        for user in user_names:
            path = (base / f"userDefs_{user}.py").resolve()
            if path not in yielded:
                yielded.add(path)
                yield path
    for path in sorted(p for p in pfw_root.glob("userDefs_*.py") if p.name != "userDefs_username.py"):
        resolved = path.resolve()
        if resolved not in yielded:
            yielded.add(resolved)
            yield resolved


def default_compiled_driver_executable(driver_name: str = DEFAULT_COMPILED_DRIVER_NAME) -> str:
    """Return the default compiled material-point driver executable.

    Precedence:
      1. GEOS_MPM_MATERIAL_POINT_DRIVER, for explicit overrides;
      2. userDefs_<user>.py materialPointDriverPath, if present;
      3. a driver executable next to userDefs.geosPath;
      4. the bare executable name, allowing PATH lookup.
    """
    env_driver = os.environ.get("GEOS_MPM_MATERIAL_POINT_DRIVER", "").strip()
    if env_driver:
        return _expand_executable_path(env_driver)

    for candidate in _iter_userdefs_candidates():
        userdefs = _load_userdefs_from(candidate)
        if userdefs is None:
            continue

        configured_driver = getattr(userdefs, "materialPointDriverPath", None)
        if not _is_unset_userdefs_value(configured_driver):
            return _expand_executable_path(configured_driver)

        geos_path = getattr(userdefs, "geosPath", None)
        if not _is_unset_userdefs_value(geos_path):
            return str(Path(_expand_executable_path(geos_path)).parent / driver_name)

    return driver_name




def default_compiled_driver_path(driver_name: str = DEFAULT_COMPILED_DRIVER_NAME) -> str:
    """Backward-compatible alias for the compiled material-point driver path.

    Older material-point verification scripts and some package initializers use
    this path-oriented name, while the driver helper that resolves the same
    executable is named ``default_compiled_driver_executable`` in this branch.
    Keep both names so mixed incremental patch application does not break
    imports.
    """
    return default_compiled_driver_executable(driver_name)


def write_compiled_driver_files(
    case: Mapping[str, Any],
    prefix: Path,
    driver_name: Optional[str] = None,
) -> Dict[str, Any]:
    """Write sidecar files for the optional compiled GEOS material-point executable.

    The compiled executable intentionally accepts a minimal GEOS material XML file and
    a path CSV rather than the Python optimizer-oriented JSON.  This helper bridges
    PFW material dictionaries to that executable without involving SolidMechanicsMPM.
    """
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if driver_name is None:
        driver_name = default_compiled_driver_executable()

    material = load_material(case.get("material", {}))
    controls = ControlProgram(case.get("control", []))
    state, update_mode = initialize_state(case, material)
    time_block = case.get("time", case.get("steps", {}))
    dt = float(time_block.get("dt", case.get("dt", 1.0)))
    n_steps = int(time_block.get("nSteps", time_block.get("nsteps", case.get("nSteps", 1))))
    initial_time = float(time_block.get("initial", 0.0))

    json_path = prefix.with_suffix(".json")
    material_xml_path = prefix.with_suffix(".material.xml")
    path_csv_path = prefix.with_suffix(".path.csv")
    run_script_path = prefix.with_suffix(".run.sh")
    output_csv_path = prefix.with_suffix(".csv")

    case_to_write = copy.deepcopy(dict(case))
    case_to_write.setdefault("materialName", material.get("name", case.get("name", "material")))
    with json_path.open("w") as stream:
        json.dump(case_to_write, stream, indent=2, sort_keys=True, default=str)
        stream.write("\n")

    material_xml_path.write_text(_material_xml_string(material))

    with path_csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("dt",) + VOIGT_COMPONENTS)
        for step_index in range(1, n_steps + 1):
            time_new = initial_time + step_index * dt
            values = [_control_value_for_compiled_driver(ctrl, time_new) for ctrl in controls.controls]
            writer.writerow([dt] + values)

    material_name = str(material.get("name", case.get("materialName", case.get("name", "material"))))
    control_arg = ",".join(_compiled_control_mode(str(ctrl.get("mode", "strain"))) for ctrl in controls.controls)
    initial = case.get("initial", {})
    temperature = case.get("temperature", {})
    energy = case.get("energy", {})
    solver = case.get("solver", case.get("newton", {}))

    command: List[str] = [
        "--material-xml", str(material_xml_path),
        "--material-name", material_name,
        "--path", str(path_csv_path),
        "--output", str(output_csv_path),
        "--control", control_arg,
        "--initial-stress", _join_reals(initial.get("stress", np.zeros(6))),
        "--material-direction", _join_reals(state.material_frame_reference.reshape(-1)),
        "--material-direction-update", str(case.get("materialDirection", {}).get("update", update_mode)),
        "--temperature-mode", str(temperature.get("mode", "prescribed")),
        "--energy-mode", str(energy.get("mode", "stressPower")),
        "--initial-temperature", str(float(initial.get("temperature", temperature.get("initial", 300.0)))),
        "--initial-temperature-rate", str(float(initial.get("temperatureRate", 0.0))),
        "--initial-specific-internal-energy", str(float(initial.get("specificInternalEnergy", initial.get("internalEnergy", 0.0)))),
        "--initial-density", str(float(initial.get("density", material.get("defaultDensity", 0.0)))),
        "--initial-length-scale", str(float(initial.get("lengthScale", initial.get("length_scale", 1.0)))),
        "--initial-strength-scale", str(float(initial.get("strengthScale", 1.0))),
        "--retention-factor", str(float(energy.get("retentionFactor", 1.0))),
        "--stress-tolerance", str(float(solver.get("stressTolerance", solver.get("absoluteTolerance", 1.0e-8)))),
        "--fd-epsilon", str(float(solver.get("finiteDifferenceStrain", solver.get("finiteDifferenceScale", 1.0e-7)))),
        "--max-newton-iterations", str(int(solver.get("maxIterations", 20))),
        "--max-line-search-iterations", str(int(solver.get("maxLineSearchIterations", 12))),
        "--max-stress-bracket-iterations", str(int(solver.get("maxStressBracketIterations", solver.get("maxBracketIterations", 32)))),
        "--max-stress-bisection-iterations", str(int(solver.get("maxStressBisectionIterations", solver.get("maxBisectionIterations", 64)))),
        "--stress-bracket-initial-scale", str(float(solver.get("stressBracketInitialScale", solver.get("bracketInitialScale", 0.0)))),
        "--stress-bracket-max-strain", str(float(solver.get("stressBracketMaxStrain", solver.get("bracketMaxStrain", 5.0e-2)))),
        "--stress-bracket-growth", str(float(solver.get("stressBracketGrowth", solver.get("bracketGrowth", 2.0)))),
        "--stress-control-failure-policy", str(solver.get("stressControlFailurePolicy", solver.get("failurePolicy", "error"))),
    ]
    if "heatCapacity" in temperature or "heatCapacity" in energy:
        command.extend(["--heat-capacity", str(float(temperature.get("heatCapacity", energy.get("heatCapacity", 1.0))))])

    core_initial_keys = {
        "stress", "strain", "totalStrain", "F", "materialDirection", "temperature", "temperatureRate",
        "density0", "density", "specificInternalEnergy", "internalEnergy", "accumulatedStressPower",
        "lengthScale", "length_scale", "strengthScale", "extra",
    }
    for key, value in sorted(initial.items()):
        if key in core_initial_keys:
            continue
        if isinstance(value, (int, float, np.integer, np.floating)):
            command.extend(["--initial-field", f"{key}={float(value):.17g}"])
        elif isinstance(value, (list, tuple)) and all(isinstance(v, (int, float, np.integer, np.floating)) for v in value):
            command.extend(["--initial-field", f"{key}={_join_reals(value)}"])

    driver_log_path = prefix.parent / f"{prefix.name}.driver.log"
    run_command = "\"$DRIVER\" \\" + "\n  " + _quote_command(command) + " > \"$LOG_FILE\" 2>&1"
    run_script_path.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        "# Set GEOS_MPM_MATERIAL_POINT_DRIVER to override the executable path.\n"
        "if [ -n \"${GEOS_MPM_MATERIAL_POINT_DRIVER:-}\" ]; then\n"
        "  DRIVER=\"$GEOS_MPM_MATERIAL_POINT_DRIVER\"\n"
        "else\n"
        f"  DRIVER={shlex.quote(str(driver_name))}\n"
        "fi\n"
        f"LOG_FILE={shlex.quote(str(driver_log_path))}\n"
        "echo \"driver log: $LOG_FILE\"\n"
        + run_command
        + "\n"
    )
    os.chmod(run_script_path, 0o755)

    return {
        "case": json_path,
        "materialXml": material_xml_path,
        "pathCsv": path_csv_path,
        "runScript": run_script_path,
        "outputCsv": output_csv_path,
        "command": command,
    }


def kroonblawd_graphite_case(
    theta_deg: float = 75.0,
    pressure_gpa: float = 30.0,
    dt: float = 5.0e-7,
    n_steps: int = 1000,
    strain_rate_us: float = -1.0e3,
    material: str = "graphiteSingleCrystal",
) -> Dict[str, Any]:
    """Return a Kroonblawd-style graphite material-point case.

    The case uses GEOS tension-positive stress convention, so a positive
    background pressure is written as a negative normal stress target.
    """
    pressure = float(pressure_gpa)
    target_stress = -pressure
    case = example_case()
    case["name"] = f"graphite_theta{float(theta_deg):g}_p{pressure:g}_material_point"
    case["material"] = {"source": "pfw_materials.py", "name": material}
    case["initial"]["stress"] = [target_stress, target_stress, target_stress, 0.0, 0.0, 0.0]
    case["materialDirection"] = {
        "type": "graphite_theta",
        "theta_deg": float(theta_deg),
        "tangentHint": [1.0, 0.0, 0.0],
        "update": "graphite",
    }
    case["time"] = {"dt": float(dt), "nSteps": int(n_steps)}
    case["control"][1]["target"] = target_stress
    case["control"][2]["rate"] = float(strain_rate_us)
    case["output"] = {"file": f"{case['name']}.csv"}
    return case


def run_material_point(
    case: Mapping[str, Any],
    executable: Optional[str] = None,
    work_dir: Optional[Union[Path, str]] = None,
    parameter_overrides: Optional[Mapping[str, Any]] = None,
    dry_run: bool = False,
    check: bool = True,
    capture_output: bool = True,
    echo_output: bool = False,
) -> CompiledDriverRunResult:
    """Prepare and optionally run the optional compiled GEOS material-point driver.

    This is the small optimizer-friendly API: callers can copy a PFW material,
    apply parameter overrides, emit isolated driver sidecar files, and run the
    compiled constitutive-point executable without touching SolidMechanicsMPM.
    """
    if executable is None:
        executable = default_compiled_driver_executable()
    executable = str(executable)

    case_copy = copy.deepcopy(dict(case))
    if parameter_overrides:
        material_block = copy.deepcopy(case_copy.get("material", {}))
        overrides = copy.deepcopy(material_block.get("overrides", {}))
        overrides.update(dict(parameter_overrides))
        material_block["overrides"] = overrides
        case_copy["material"] = material_block

    name = str(case_copy.get("name", "material_point"))
    root = Path(work_dir) if work_dir is not None else Path(name)
    root.mkdir(parents=True, exist_ok=True)
    files = write_compiled_driver_files(case_copy, root / name, executable)

    command = [str(executable)] + [str(part) for part in files["command"]]
    log_file = root / f"{name}.driver.log"
    result = CompiledDriverRunResult(
        case=Path(files["case"]),
        material_xml=Path(files["materialXml"]),
        load_path_csv=Path(files["pathCsv"]),
        run_script=Path(files["runScript"]),
        output_csv=Path(files["outputCsv"]),
        command=command,
        log_file=log_file,
    )

    if not dry_run:
        if capture_output:
            proc = subprocess.run(
                command,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            log_text = (
                "# Command\n" + _quote_command(command) +
                "\n\n# stdout\n" + proc.stdout +
                "\n# stderr\n" + proc.stderr
            )
            log_file.write_text(log_text)
            if echo_output or (check and proc.returncode != 0):
                if proc.stdout:
                    print(proc.stdout, end="")
                if proc.stderr:
                    print(proc.stderr, end="", file=sys.stderr)
        else:
            proc = subprocess.run(command, check=False)
        result.returncode = proc.returncode
        if check and proc.returncode != 0:
            detail = f"; driver log: {log_file}" if capture_output else ""
            raise BackendError(f"Compiled GEOS material-point driver returned {proc.returncode}{detail}")
    return result

def example_case() -> Dict[str, Any]:
    """Return a Kroonblawd-style graphite local-driver case template."""
    return {
        "name": "graphite_theta75_p30_material_point",
        "units": "mm_us_mg_GPa_K",
        "backend": "elastic",
        "material": {"source": "pfw_materials.py", "name": "graphiteSingleCrystal"},
        "initial": {
            "stress": [-30.0, -30.0, -30.0, 0.0, 0.0, 0.0],
            "temperature": 300.0,
            "specificInternalEnergy": 0.0,
            "lengthScale": 1.0e-4,
            "strengthScale": 1.0,
            "damage": 0.0,
            "basalPlaneDamage": 0.0,
            "comminutionDamage": 0.0,
        },
        "materialDirection": {
            "type": "graphite_theta",
            "theta_deg": 75.0,
            "tangentHint": [1.0, 0.0, 0.0],
            "update": "graphite",
        },
        "time": {"dt": 5.0e-7, "nSteps": 1000},
        "temperature": {"mode": "isothermal", "initial": 300.0},
        "energy": {"mode": "stressPower", "retentionFactor": 0.0, "outputAccumulatedStressPower": True},
        "solver": {"stressTolerance": 1.0e-7, "maxIterations": 20, "finiteDifferenceStrain": 1.0e-7},
        "control": [
            {"component": "xx", "mode": "strain", "value": 0.0},
            {"component": "yy", "mode": "stress", "target": -30.0},
            {"component": "zz", "mode": "trueStrainRate", "rate": -1.0e3},
            {"component": "yz", "mode": "stress", "target": 0.0},
            {"component": "xz", "mode": "strain", "value": 0.0},
            {"component": "xy", "mode": "strain", "value": 0.0},
        ],
        "output": {"file": "graphite_theta75_p30_material_point.csv"},
    }


def read_case(path: Path) -> Dict[str, Any]:
    with path.open() as f:
        return json.load(f)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Run a standalone GEOS-MPM material-point control path.")
    parser.add_argument("--case", type=Path, help="JSON material-point case file")
    parser.add_argument("--output", type=Path, help="CSV output file; overrides case output.file")
    parser.add_argument("--write-example", type=Path, help="Write an example JSON case and exit")
    parser.add_argument("--backend", help="Override backend type, e.g. elastic or external-json")
    parser.add_argument("--export-compiled-driver-prefix", type=Path,
                        help="Write .material.xml, .path.csv, .run.sh, and normalized .json files for the optional compiled GEOS driver, then exit")
    parser.add_argument("--compiled-driver-name", default=None,
                        help="Executable name/path embedded in the generated .run.sh file. Default: GEOS_MPM_MATERIAL_POINT_DRIVER, then userDefs.geosPath-adjacent driver, then PATH lookup.")
    args = parser.parse_args(argv)

    if args.write_example:
        args.write_example.parent.mkdir(parents=True, exist_ok=True)
        with args.write_example.open("w") as f:
            json.dump(example_case(), f, indent=2)
            f.write("\n")
        return 0

    if not args.case:
        parser.error("--case is required unless --write-example is used")

    case = read_case(args.case)
    if args.backend:
        case["backend"] = args.backend
    if args.export_compiled_driver_prefix:
        files = write_compiled_driver_files(case, args.export_compiled_driver_prefix, args.compiled_driver_name)
        for label, path in files.items():
            if isinstance(path, Path):
                print(f"{label}: {path}")
        return 0
    driver = MaterialPointDriver(case)
    results = driver.run()
    output_file = args.output or Path(case.get("output", {}).get("file", f"{case.get('name', 'material_point')}.csv"))
    write_csv(results, output_file, driver.controls)
    print(f"Wrote {len(results)} material-point steps to {output_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
