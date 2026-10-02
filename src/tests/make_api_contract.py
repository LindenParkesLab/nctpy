"""Generate fixtures/api_contract.json: the frozen public signatures of nctpy and null_models.

The public API is frozen within 1.x: the Nature Protocols paper (and its Supplementary
Information) prints code that must keep running unmodified. This script records every public
function and class signature so that test_api_contract.py can detect any change.

Signatures are read from the source with ``ast`` rather than by importing the modules, because
importing nctpy.plotting fetches fsaverage over the network (a default argument). The JSON uses
inspect's parameter-kind names so the test can compare it directly with ``inspect.signature``.

Defaults are stored as ``{"value": ...}`` when they are literals and ``{"expr": "..."}`` when they
are not (e.g. surface_plot's fsaverage); the test checks literal values exactly and non-literal
defaults for presence only.

Run from src/tests, and only when the API is deliberately changed (which within 1.x means a new
keyword argument whose default reproduces existing behaviour):

    python make_api_contract.py
"""
import ast
import json
from pathlib import Path

SRC = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "fixtures" / "api_contract.json"

MODULES = {
    "nctpy.energies": "nctpy/energies.py",
    "nctpy.metrics": "nctpy/metrics.py",
    "nctpy.pipelines": "nctpy/pipelines.py",
    "nctpy.utils": "nctpy/utils.py",
    "nctpy.plotting": "nctpy/plotting.py",
    "nctpy.optimize": "nctpy/optimize.py",
    "null_models.geomsurr": "null_models/geomsurr.py",
}

# Provenance of each contract symbol:
#   P  printed and called in the protocol paper or its SI
#   Pi imported by the paper's printed code but never called there
#   D  used in the published documentation or the notebooks under scripts/
#   X  called by nct_xr
CONTRACT = {
    "nctpy.energies": {
        "sim_state_eq": "D",
        "get_control_inputs": "P D X",
        "integrate_u": "P D X",
        "gramian": "D",
        "minimum_energy_fast": "D",
    },
    "nctpy.metrics": {
        "ave_control": "P D",
        "modal_control": "D",
    },
    "nctpy.pipelines": {
        "ComputeControlEnergy": "P D",
        "ComputeOptimizedControlEnergy": "P D",
    },
    "nctpy.utils": {
        "matrix_normalization": "P D X",
        "normalize_state": "P D X",
        "normalize_weights": "P D",
        "convert_states_str2int": "P D",
        "expand_states": "D",
        "get_null_p": "P D",
        "get_fdr_p": "P D",
        "get_p_val_string": "D",
        "expm": "D",
    },
    "nctpy.plotting": {
        "set_plotting_params": "D",
        "reg_plot": "D",
        "null_plot": "P D",
        "roi_to_vtx": "Pi D",
        "surface_plot": "P D",
        "add_module_lines": "Pi D",
    },
    "null_models.geomsurr": {
        "geomsurr": "P D",
    },
}


def _default(node):
    try:
        return {"value": ast.literal_eval(node)}
    except ValueError:
        return {"expr": ast.unparse(node)}


def _params(fn, skip_self=False):
    a = fn.args
    params = []
    positional = [(p, "POSITIONAL_ONLY") for p in a.posonlyargs] + \
                 [(p, "POSITIONAL_OR_KEYWORD") for p in a.args]
    # positional defaults align with the end of the positional list
    pos_defaults = [None] * (len(positional) - len(a.defaults)) + list(a.defaults)
    for (p, kind), d in zip(positional, pos_defaults):
        params.append({"name": p.arg, "kind": kind, **({"default": _default(d)} if d is not None else {})})
    if a.vararg:
        params.append({"name": a.vararg.arg, "kind": "VAR_POSITIONAL"})
    for p, d in zip(a.kwonlyargs, a.kw_defaults):
        params.append({"name": p.arg, "kind": "KEYWORD_ONLY", **({"default": _default(d)} if d is not None else {})})
    if a.kwarg:
        params.append({"name": a.kwarg.arg, "kind": "VAR_KEYWORD"})
    if skip_self:
        params = params[1:]
    return params


def _symbols(module, path):
    tree = ast.parse((SRC / path).read_text())
    listed = CONTRACT.get(module, {})
    out = []
    for node in tree.body:
        if not isinstance(node, (ast.FunctionDef, ast.ClassDef)) or node.name.startswith("_"):
            continue
        entry = {"module": module, "name": node.name,
                 "tier": "contract" if node.name in listed else "public-unlisted"}
        if node.name in listed:
            entry["sources"] = listed[node.name].split()
        if isinstance(node, ast.FunctionDef):
            entry["kind"] = "function"
            entry["params"] = _params(node)
        else:
            entry["kind"] = "class"
            methods = {m.name: m for m in node.body if isinstance(m, ast.FunctionDef)}
            if "__init__" in methods:
                entry["params"] = _params(methods["__init__"], skip_self=True)
            else:  # a dataclass: its constructor takes the annotated fields, in order
                entry["params"] = [
                    {"name": f.target.id, "kind": "POSITIONAL_OR_KEYWORD",
                     **({"default": _default(f.value)} if f.value is not None else {})}
                    for f in node.body if isinstance(f, ast.AnnAssign) and isinstance(f.target, ast.Name)
                ]
            entry["methods"] = {name: _params(m, skip_self=True)
                                for name, m in methods.items() if not name.startswith("_")}
        out.append(entry)
    missing = set(listed) - {e["name"] for e in out}
    if missing:
        raise SystemExit(f"{module}: contract symbols not found in source: {sorted(missing)}")
    return out


def main():
    symbols = [s for module, path in MODULES.items() for s in _symbols(module, path)]
    OUT.write_text(json.dumps({"symbols": symbols}, indent=2) + "\n")
    n_contract = sum(s["tier"] == "contract" for s in symbols)
    print(f"wrote {OUT.relative_to(SRC.parent)}: {n_contract} contract symbols, "
          f"{len(symbols) - n_contract} public-unlisted")


if __name__ == "__main__":
    main()
