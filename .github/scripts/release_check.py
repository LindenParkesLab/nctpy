"""Check that a release tag agrees with the repository, and write the release notes.

Used by .github/workflows/release.yml; runnable locally before tagging:

    python .github/scripts/release_check.py v1.1.0 release_notes.md

For tag vX.Y.Z it checks that `__version__` in src/nctpy/__init__.py, `version` in CITATION.cff and
a dated `## X.Y.Z (YYYY-MM-DD)` section in CHANGELOG.md all say X.Y.Z, and that CITATION.cff's
`date-released` matches the changelog date. It writes that changelog section to the notes file.
"""

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def read(path):
    return (REPO / path).read_text(encoding="utf-8")


def main(tag, notes_path):
    errors = []
    version = tag[1:] if tag.startswith("v") else tag
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        errors.append(f"tag {tag!r} is not of the form vX.Y.Z")

    found = re.search(r'^__version__ = "([^"]+)"', read("src/nctpy/__init__.py"), re.M)
    if not found or found.group(1) != version:
        errors.append(f"src/nctpy/__init__.py has __version__ {found and found.group(1)!r}, not {version!r}")

    citation = read("CITATION.cff")
    cff_version = re.search(r"^version: *(\S+)", citation, re.M)
    if not cff_version or cff_version.group(1).strip("\"'") != version:
        errors.append(f"CITATION.cff has version {cff_version and cff_version.group(1)!r}, not {version!r}")

    changelog = read("CHANGELOG.md")
    heading = re.search(rf"^## {re.escape(version)} \((\d{{4}}-\d{{2}}-\d{{2}})\)\s*$", changelog, re.M)
    if not heading:
        errors.append(f"CHANGELOG.md has no dated section '## {version} (YYYY-MM-DD)'")
    else:
        cff_date = re.search(r"^date-released: *(\S+)", citation, re.M)
        if not cff_date or cff_date.group(1).strip("\"'") != heading.group(1):
            errors.append(
                f"CITATION.cff date-released {cff_date and cff_date.group(1)!r} "
                f"does not match the changelog date {heading.group(1)!r}"
            )
        following = re.search(r"^## ", changelog[heading.end() :], re.M)
        section = changelog[heading.end() : heading.end() + following.start() if following else None]
        Path(notes_path).write_text(section.strip() + "\n", encoding="utf-8")

    if errors:
        print("release check failed:\n  " + "\n  ".join(errors))
        return 1
    print(f"release check passed for {version}; notes written to {notes_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
