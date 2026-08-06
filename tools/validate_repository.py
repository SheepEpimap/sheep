#!/usr/bin/env python3
"""Run repository-level static checks using only the Python standard library."""

from __future__ import annotations

import argparse
import ast
import csv
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


TEXT_SUFFIXES = {
    ".md",
    ".py",
    ".r",
    ".sh",
    ".sbatch",
    ".smk",
    ".txt",
    ".tsv",
    ".csv",
    ".yml",
    ".yaml",
    ".toml",
    ".env",
    ".ipynb",
}


@dataclass
class Check:
    name: str
    status: str
    detail: str


def run(command: list[str], cwd: Path) -> tuple[int, str]:
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return completed.returncode, completed.stdout.strip()


def check_python(root: Path) -> Check:
    files = sorted(root.rglob("*.py"))
    errors: list[str] = []
    for path in files:
        try:
            ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
        except (SyntaxError, UnicodeError) as exc:
            errors.append(f"{path.relative_to(root)}: {exc}")
    if errors:
        return Check("Python AST", "FAIL", "; ".join(errors))
    return Check("Python AST", "PASS", f"{len(files)} files parsed")


def check_notebooks(root: Path) -> Check:
    files = sorted(root.rglob("*.ipynb"))
    errors: list[str] = []
    for path in files:
        try:
            payload = json.loads(path.read_text(encoding="utf-8-sig"))
            if not isinstance(payload.get("cells"), list):
                errors.append(f"{path.relative_to(root)}: missing cells list")
        except (json.JSONDecodeError, UnicodeError) as exc:
            errors.append(f"{path.relative_to(root)}: {exc}")
    if errors:
        return Check("Notebook JSON", "FAIL", "; ".join(errors))
    return Check("Notebook JSON", "PASS", f"{len(files)} notebooks parsed")


def check_bash(root: Path) -> Check:
    files = sorted([*root.rglob("*.sh"), *root.rglob("*.sbatch")])
    candidates = [shutil.which("bash")]
    git = shutil.which("git")
    if git:
        git_sh = Path(git).resolve().parents[1] / "usr" / "bin" / "sh.exe"
        if git_sh.is_file():
            candidates.append(str(git_sh))
    candidates.append(shutil.which("sh"))
    bash = None
    for candidate in dict.fromkeys(item for item in candidates if item):
        version_code, version_output = run([candidate, "--version"], root)
        if not version_code and "bash" in version_output.lower():
            bash = candidate
            break
    if bash is None:
        return Check("Bash syntax", "SKIP", f"bash unavailable; {len(files)} files found")
    errors: list[str] = []
    for path in files:
        code, output = run([bash, "-n", str(path)], root)
        if code:
            errors.append(f"{path.relative_to(root)}: {output}")
    if errors:
        return Check("Bash syntax", "FAIL", "; ".join(errors))
    return Check("Bash syntax", "PASS", f"{len(files)} files passed bash -n")


def check_r(root: Path) -> Check:
    files = sorted(root.rglob("*.R"))
    rscript = shutil.which("Rscript")
    if not rscript and sys.platform == "win32":
        candidates = sorted(Path("C:/Program Files/R").glob("R-*/bin/Rscript.exe"), reverse=True)
        if candidates:
            rscript = str(candidates[0])
    if not rscript:
        return Check("R syntax", "SKIP", f"Rscript unavailable; {len(files)} files found")
    errors: list[str] = []
    for path in files:
        expression = f"parse(file={str(path)!r})"
        code, output = run([rscript, "-e", expression], root)
        if code:
            errors.append(f"{path.relative_to(root)}: {output}")
    if errors:
        return Check("R syntax", "FAIL", "; ".join(errors))
    return Check("R syntax", "PASS", f"{len(files)} files parsed")


def check_snakemake(root: Path) -> Check:
    files = sorted(
        path
        for path in root.rglob("*")
        if path.is_file() and (path.name == "Snakefile" or path.suffix == ".smk")
    )
    snakemake = shutil.which("snakemake")
    if not snakemake:
        return Check("Snakemake syntax", "SKIP", f"snakemake unavailable; {len(files)} files found")
    errors: list[str] = []
    for path in files:
        code, output = run([snakemake, "-s", str(path), "--lint"], path.parent)
        if code:
            errors.append(f"{path.relative_to(root)}: {output}")
    if errors:
        return Check("Snakemake syntax", "FAIL", "; ".join(errors))
    return Check("Snakemake syntax", "PASS", f"{len(files)} files passed snakemake --lint")


def check_figure_map(root: Path) -> Check:
    map_path = root / "FIGURE_CODE_MAP.tsv"
    if not map_path.is_file():
        return Check("Figure map", "FAIL", "FIGURE_CODE_MAP.tsv is missing")
    missing: list[str] = []
    invalid_relations: list[str] = []
    allowed_relations = {
        "PANEL_CODE",
        "PANEL_INPUT",
        "PARTIAL_PANEL",
        "VALIDATION",
        "ANALYSIS_SUPPORT",
    }
    rows = 0
    with map_path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows += 1
            script = row.get("script", "")
            relation = row.get("relation", "")
            if not script or not (root / script).is_file():
                missing.append(script or f"row {rows}: empty script")
            if relation not in allowed_relations:
                invalid_relations.append(f"{script}: {relation}")
    if missing:
        return Check("Figure map", "FAIL", f"missing mapped files: {', '.join(missing)}")
    if invalid_relations:
        return Check("Figure map", "FAIL", f"invalid relations: {', '.join(invalid_relations)}")
    return Check("Figure map", "PASS", f"{rows} mapped scripts exist and use approved relationship labels")


def check_markdown_links(root: Path) -> Check:
    pattern = re.compile(r"(?<!!)\[[^\]]+\]\((?:<([^>]+)>|([^)]+))\)")
    missing: list[str] = []
    checked = 0
    for path in sorted(root.rglob("*.md")):
        text = path.read_text(encoding="utf-8-sig")
        for match in pattern.finditer(text):
            raw_target = (match.group(1) or match.group(2)).strip()
            target = raw_target.split(maxsplit=1)[0] if match.group(2) else raw_target
            if not target or target.startswith(("#", "http://", "https://", "mailto:")):
                continue
            checked += 1
            local = target.split("#", 1)[0]
            if local and not (path.parent / local).resolve().exists():
                missing.append(f"{path.relative_to(root)} -> {target}")
    if missing:
        return Check("Markdown links", "FAIL", "; ".join(missing))
    return Check("Markdown links", "PASS", f"{checked} local links resolved")


def check_text_encoding(root: Path) -> Check:
    failures: list[str] = []
    count = 0
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES and path.name not in {"Snakefile", ".gitignore", ".gitattributes"}:
            continue
        count += 1
        try:
            path.read_text(encoding="utf-8-sig")
        except UnicodeError as exc:
            failures.append(f"{path.relative_to(root)}: {exc}")
    if failures:
        return Check("UTF-8 text", "FAIL", "; ".join(failures))
    return Check("UTF-8 text", "PASS", f"{count} text files decoded")


def check_large_files(root: Path) -> Check:
    limit = 100 * 1024 * 1024
    repository_files = [
        path
        for path in root.rglob("*")
        if path.is_file() and ".git" not in path.relative_to(root).parts
    ]
    large = [
        f"{path.relative_to(root)} ({path.stat().st_size} bytes)"
        for path in repository_files
        if path.stat().st_size >= limit
    ]
    if large:
        return Check("GitHub file size", "FAIL", "; ".join(large))
    return Check("GitHub file size", "PASS", f"{len(repository_files)} files; none >= 100 MiB")


def check_secrets(root: Path) -> Check:
    patterns = {
        "AWS access key": re.compile(r"\bAKIA[0-9A-Z]{16}\b"),
        "GitHub token": re.compile(r"\bgh[pousr]_[A-Za-z0-9]{30,}\b"),
        "OpenAI-style key": re.compile(r"\bsk-[A-Za-z0-9_-]{20,}\b"),
        "private key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    }
    hits: list[str] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES and path.name not in {"Snakefile", ".gitignore", ".gitattributes"}:
            continue
        text = path.read_text(encoding="utf-8-sig")
        for label, pattern in patterns.items():
            if pattern.search(text):
                hits.append(f"{path.relative_to(root)}: {label}")
    if hits:
        return Check("Credential patterns", "FAIL", "; ".join(hits))
    return Check("Credential patterns", "PASS", "no high-confidence credential patterns found")


def audit_absolute_paths(root: Path) -> Check:
    patterns = ["/vol2/", "/public/home/", "/storage/public/home/", "/data/home/"]
    counts = dict.fromkeys(patterns, 0)
    files: set[Path] = set()
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES and path.name != "Snakefile":
            continue
        text = path.read_text(encoding="utf-8-sig")
        for pattern in patterns:
            count = text.count(pattern)
            if count:
                counts[pattern] += count
                files.add(path)
    detail = ", ".join(f"{key}={value}" for key, value in counts.items())
    return Check("Absolute path audit", "INFO", f"{detail}; {len(files)} files affected")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", nargs="?", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()

    checks = [
        check_python(root),
        check_notebooks(root),
        check_bash(root),
        check_r(root),
        check_snakemake(root),
        check_figure_map(root),
        check_markdown_links(root),
        check_text_encoding(root),
        check_large_files(root),
        check_secrets(root),
        audit_absolute_paths(root),
    ]

    print("| Check | Status | Detail |")
    print("|---|---|---|")
    for check in checks:
        detail = check.detail.replace("|", "\\|").replace("\n", " ")
        print(f"| {check.name} | {check.status} | {detail} |")

    failures = [check for check in checks if check.status == "FAIL"]
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
