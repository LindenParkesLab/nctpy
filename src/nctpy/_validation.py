"""Input checks shared across nctpy (private).

These reject input that does not define a problem. They never reject ill-conditioned or incomplete
control problems: those return their energies and error terms, as the protocol paper describes.
"""


def _check_system(system, see="matrix_normalization help"):
    """Raise unless `system` is 'continuous' or 'discrete'.

    The exception type (bare Exception) and the messages are part of the frozen API: downstream code may
    catch Exception, and tests compare the messages.
    """
    if system is None:
        raise Exception(
            "Time system not specified. "
            "Please nominate whether you are normalizing A for a continuous-time or a discrete-time system "
            f"(see {see})."
        )
    if system != "continuous" and system != "discrete":
        raise Exception(
            "Incorrect system specification. Please specify either 'system=discrete' or 'system=continuous'."
        )


def _check_rho(rho):
    """Raise ValueError unless rho > 0.

    rho weights the cost of the control inputs and enters the solver as -B B^T / (2 rho), so rho <= 0 does not
    define an optimisation problem. rho = 0 used to return NaN energies and error terms with numpy warnings.
    """
    if not rho > 0:  # also rejects NaN
        raise ValueError(
            f"rho must be positive, got {rho!r}. rho weights the cost of the control inputs against the cost of "
            "the state trajectory; with S all zeros (minimum-energy control), every positive rho gives the same "
            "result."
        )
