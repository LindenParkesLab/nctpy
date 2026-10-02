"""Optimise nodes' decay rates (intrinsic neural timescales) for a state transition.

Implements the method of Kim et al., Nat Commun 16:11639 (2025): each node's decay rate is a learnable parameter,
fitted by gradient descent so that the optimal trajectory from x0 to xf passes close to the reference state at the
midpoint of the time horizon, while every eigenvalue of the fitted system stays negative. Equation numbers refer to
that paper's Methods. Needs PyTorch: ``pip install "nctpy[optimize]"``.
"""

from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt

try:
    import torch
except ImportError as exc:
    raise ImportError(
        "nctpy.optimize needs PyTorch, which is not installed with nctpy by default. Install it with:\n"
        "    pip install 'nctpy[optimize]'"
    ) from exc
from tqdm import tqdm

from nctpy._validation import _check_rho
from nctpy.energies import _as_float, _check_reference, _column, _reference
from nctpy.utils import matrix_normalization


@dataclass(frozen=True)
class DecayRateFit:
    """The result of :func:`optimize_decay_rates`.

    Attributes
    ----------
    decay : (N,) ndarray
        The fitted decay rate of each node, in the convention of :func:`nctpy.utils.matrix_normalization`:
        ``matrix_normalization(A, "continuous", c=c, zero_diagonal=zero_diagonal, decay=fit.decay)`` gives
        ``A_norm`` (to rounding). A value v reported in Kim et al. (2025), Fig. 2B is ``1 - decay``.
    A_norm : (N, N) ndarray
        The fitted, normalised connectivity matrix. Use it directly (e.g. with
        :func:`nctpy.energies.get_control_inputs`) rather than rebuilding it from the decay rates.
    w : (N,) ndarray
        The trained parameter itself: ``A_norm`` is the normalised connectome minus ``diag(w)``, so
        ``decay = 1 + w``. Kim et al.'s code saves this as ``optimized_weights``.
    loss : (n_steps_run,) ndarray
        The loss after each step.
    max_eigenvalue : (n_steps_run,) ndarray
        The largest eigenvalue of the fitted matrix after each step, as used by the stability penalty.
    n_steps_run : int
        Number of steps run.
    stopped_early : bool
        Whether training stopped before `n_steps` because the loss and the largest eigenvalue had stopped changing.
    """

    decay: npt.NDArray[np.float64]
    A_norm: npt.NDArray[np.float64]
    w: npt.NDArray[np.float64]
    loss: npt.NDArray[np.float64]
    max_eigenvalue: npt.NDArray[np.float64]
    n_steps_run: int
    stopped_early: bool


def optimize_decay_rates(
    A: npt.ArrayLike,
    x0: npt.ArrayLike,
    xf: npt.ArrayLike,
    *,
    T: float = 1,
    B: npt.ArrayLike | None = None,
    S: npt.ArrayLike | None = None,
    rho: float = 1,
    xr: npt.ArrayLike | str = "xf",
    c: float = 1,
    zero_diagonal: bool = True,
    init: str = "one",
    n_steps: int = 1000,
    lr: float = 0.01,
    eig_weight: float = 1.0,
    reg_weight: float = 1e-4,
    reg_type: str = "l2",
    eigenvalues: str = "symmetric",
    early_stopping: bool = True,
    device: str | None = None,
    progress: bool = True,
) -> DecayRateFit:
    """Fit each node's decay rate so that the optimal trajectory from x0 to xf passes through xr at time T/2.

    The connectome is normalised for a continuous-time system (``A / (c + λ_max) - I``), and a decay term
    ``diag(w)`` is subtracted from it (Eq. 4), starting from ``w = 1`` (``init="one"``). Adam then minimises the sum
    of three terms:

    - the distance ``||x*(T/2) - xr||`` between the optimal trajectory at the midpoint and the reference state
      (Eq. 10), averaged over transitions; the trajectory comes from the state-costate solution (Eqs. 6-9);
    - a stability penalty, ``mean(-eig_weight / λ_i)`` over the fitted matrix's eigenvalues, which grows without
      bound as any eigenvalue approaches 0 (Eq. 11);
    - a penalty on the size of the diagonal, ``reg_weight * sum(diag^2)`` (l2) or ``sum(|diag|)`` (l1).

    The defaults are those used by Kim et al. (2025).

    Parameters
    ----------
    A : (N, N) array_like
        Structural connectome (not normalised).
    x0, xf : (N,) or (N, k) array_like
        Initial and target states: one transition, or k transitions fitted together as columns.
    T : float, default 1
        Time horizon.
    B, S : (N, N) array_like, optional
        Control set and state-trajectory constraints; default the identity.
    rho : float, default 1
        Weight of the input cost relative to the trajectory cost. Must be positive.
    xr : (N,) or (N, k) array_like, or {'xf', 'x0', 'zero', 'midpoint'}, default 'xf'
        Reference state, used both in the cost and as the midpoint target.
    c : float, default 1
        Normalisation constant.
    zero_diagonal : bool, default True
        Remove self-connections before normalising. The method assumes a connectome without self-connections, so
        that every node starts from the same decay rate; it is supported only for such connectomes. With False, a
        connectome's self-connections are kept, at the user's own discretion.
    init : {'one', 'zero'}, default 'one'
        Starting value of w: 'one' starts every node at decay rate 2, 'zero' at decay rate 1.
    n_steps : int, default 1000
        Maximum number of gradient steps.
    lr : float, default 0.01
        Adam's learning rate.
    eig_weight : float, default 1.0
        Weight of the stability penalty.
    reg_weight : float, default 1e-4
        Weight of the penalty on the diagonal.
    reg_type : {'l2', 'l1'}, default 'l2'
        Form of the penalty on the diagonal.
    eigenvalues : {'symmetric', 'general'}, default 'symmetric'
        How the stability penalty and the tracked largest eigenvalue are computed. 'symmetric' uses the routine for
        symmetric matrices, as Kim et al. (2025) did. It gives the exact eigenvalues of an undirected connectome;
        for a directed one it uses only the lower triangle, so it is an approximation. 'general' uses the real parts
        of the eigenvalues of the matrix as it is, which is exact for directed connectomes too.
    early_stopping : bool, default True
        Stop once the loss and the largest eigenvalue, rounded to 4 decimals, have not varied over the last 20% of
        `n_steps` steps.
    device : str, optional
        PyTorch device, e.g. 'cpu' or 'cuda'. Default: 'cuda' if available, else 'cpu'. Results on different
        devices agree to floating-point precision.
    progress : bool, default True
        Show a progress bar.

    Returns
    -------
    DecayRateFit
        The fitted decay rates and matrix, and the training traces.
    """
    for name, value, options in (
        ("init", init, ("one", "zero")),
        ("reg_type", reg_type, ("l2", "l1")),
        ("eigenvalues", eigenvalues, ("symmetric", "general")),
    ):
        if value not in options:
            raise ValueError(f"{name} must be one of {options}; got {value!r}")
    _check_rho(rho)

    A_norm = matrix_normalization(A, system="continuous", c=c, zero_diagonal=zero_diagonal)
    n_nodes = A_norm.shape[0]
    X0, XF = _column(_as_float(x0)), _column(_as_float(xf))
    XR = _reference(xr, X0, XF, n_nodes)
    _check_reference(XR)
    B = np.eye(n_nodes) if B is None else _as_float(B)
    S = np.eye(n_nodes) if S is None else _as_float(S)

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    model = _DecayRateModel(A_norm, X0, XF, XR, T, B, S, rho, init, eig_weight, reg_weight, reg_type, eigenvalues)
    model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    model.train()

    max_eig = (lambda M: np.max(np.linalg.eigvalsh(M))) if eigenvalues == "symmetric" else _max_real_eigenvalue
    window = int(n_steps * 0.20)
    losses, max_eigenvalues, ws = [], [], []
    stopped_early = False
    loss = model()
    for i in tqdm(range(n_steps), disable=not progress):
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        loss = model()  # the loss after this step, whose graph also gives the next step's gradient
        losses.append(loss.item())
        ws.append(model.w.detach().cpu().numpy().flatten())
        max_eigenvalues.append(max_eig(A_norm - np.diag(ws[-1])))

        if early_stopping and i > window:
            loss_var = np.round(np.var(losses[i - window : i]), 4)
            eig_var = np.round(np.var(max_eigenvalues[i - window : i]), 4)
            if loss_var == 0 and eig_var == 0:
                stopped_early = True
                break

    w = ws[-1]
    return DecayRateFit(
        decay=1 + w,
        A_norm=A_norm - np.diag(w),
        w=w,
        loss=np.array(losses),
        max_eigenvalue=np.array(max_eigenvalues),
        n_steps_run=len(losses),
        stopped_early=stopped_early,
    )


def _max_real_eigenvalue(M: npt.NDArray[np.float64]) -> float:
    return float(np.max(np.linalg.eigvals(M).real))


class _DecayRateModel(torch.nn.Module):
    """The loss of Kim et al. (2025) as a function of the decay parameter w; float64 throughout."""

    def __init__(
        self,
        A_norm: npt.NDArray[np.float64],
        X0: npt.NDArray[np.float64],
        XF: npt.NDArray[np.float64],
        XR: npt.NDArray[np.float64],
        T: float,
        B: npt.NDArray[np.float64],
        S: npt.NDArray[np.float64],
        rho: float,
        init: str,
        eig_weight: float,
        reg_weight: float,
        reg_type: str,
        eigenvalues: str,
    ) -> None:
        super().__init__()
        n_nodes = A_norm.shape[0]
        start = np.ones((n_nodes, 1)) if init == "one" else np.zeros((n_nodes, 1))
        self.w = torch.nn.Parameter(torch.from_numpy(start))
        for name, value in (("A_norm", A_norm), ("X0", X0), ("XF", XF), ("XR", XR), ("B", B), ("S", S)):
            self.register_buffer(name, torch.from_numpy(np.ascontiguousarray(value, dtype=np.float64)))
        self.n_nodes = n_nodes
        self.T, self.rho = T, rho
        self.eig_weight, self.reg_weight, self.reg_type = eig_weight, reg_weight, reg_type
        self.eigenvalues = eigenvalues

    def forward(self) -> Any:
        n, dtype, device = self.n_nodes, torch.float64, self.A_norm.device
        A_sub = self.A_norm - torch.diag(self.w[:, 0])

        # stability penalty (Eq. 11) and the penalty on the diagonal
        if self.eigenvalues == "symmetric":
            eigvals = torch.linalg.eigvalsh(A_sub)
        else:
            eigvals = torch.linalg.eigvals(A_sub).real
        eig_loss = torch.mean(-torch.div(self.eig_weight, eigvals))
        diagonal = torch.diag(A_sub)
        reg = self.reg_weight * (torch.square(diagonal).sum() if self.reg_type == "l2" else torch.abs(diagonal).sum())
        loss = reg + eig_loss

        # the state-costate system (Eq. 6) and its propagator over [0, T/2] and [0, T] (Eqs. 7-9)
        M = torch.concat(
            (
                torch.concat((A_sub, torch.matmul(-self.B, self.B.T) / (2 * self.rho)), dim=1),
                torch.concat((-2 * self.S, -A_sub.T), dim=1),
            ),
            dim=0,
        )
        E_half = torch.matrix_exp(M * (self.T / 2))
        E = torch.matrix_power(E_half, 2)
        r = torch.arange(n, device=device)
        E11, E12 = E[r, :][:, r], E[r, :][:, r + n]
        eye = torch.eye(n, dtype=dtype, device=device)
        state_rows = torch.concat((eye, torch.zeros((n, n), dtype=dtype, device=device)), dim=1)

        n_transitions = self.X0.shape[1]
        for j in range(n_transitions):
            x0, xf, xr = self.X0[:, j : j + 1], self.XF[:, j : j + 1], self.XR[:, j : j + 1]
            ref_input = torch.concat((torch.zeros((n, 1), dtype=dtype, device=device), 2 * torch.matmul(self.S, xr)))
            ref_offset = torch.linalg.solve(M, ref_input)
            p0 = torch.linalg.solve(
                E12, xf - torch.matmul(E11, x0) - torch.matmul(torch.concat((E11 - eye, E12), dim=1), ref_offset)
            )
            z0 = torch.concat((x0, p0), dim=0)
            x_mid = torch.matmul(E_half[r, :], z0) + torch.matmul(E_half[r, :] - state_rows, ref_offset)  # x*(T/2)
            loss = loss + torch.linalg.norm(x_mid - xr) / n_transitions  # Eq. 10
        return loss
