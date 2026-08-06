#!/usr/bin/env python3
"""Rebuild pipelines as one folder per H1 and one Markdown file per H2/H3."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from dataclasses import dataclass
from pathlib import Path


HEADING_RE = re.compile(r"^(#{1,6})\s+(.+?)\s*$")
INVALID_WINDOWS_CHARS = re.compile(r'[<>:"/\\|?*]')


@dataclass(frozen=True)
class OutputSection:
    path: Path
    heading: str
    level: int
    start: int
    end: int
    content: str


def clean_heading(raw: str) -> str:
    value = raw.strip()
    value = re.sub(r"^\*\*(.+)\*\*$", r"\1", value)
    return value.strip()


def safe_name(raw: str) -> str:
    value = clean_heading(raw)
    value = INVALID_WINDOWS_CHARS.sub("_", value)
    value = re.sub(r"_+", "_", value)
    value = value.rstrip(" .")
    return value or "untitled"


def h1_basename(title: str) -> str:
    value = re.sub(r"^\d+\.\s*", "", clean_heading(title))
    return safe_name(value)


def has_code_fence(lines: list[str]) -> bool:
    return any(line.lstrip().startswith("```") for line in lines)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def source_slice(lines: list[str], start_index: int, end_index: int) -> str:
    return "".join(lines[start_index:end_index])


def build_sections(source: Path, destination: Path) -> list[OutputSection]:
    lines = source.read_text(encoding="utf-8-sig").splitlines(keepends=True)
    headings: list[tuple[int, int, str]] = []
    in_fence = False
    fence_marker = ""
    for index, line in enumerate(lines):
        stripped = line.lstrip()
        if stripped.startswith("```") or stripped.startswith("~~~"):
            marker = stripped[:3]
            if not in_fence:
                in_fence = True
                fence_marker = marker
            elif marker == fence_marker:
                in_fence = False
                fence_marker = ""
            continue
        if in_fence:
            continue
        match = HEADING_RE.match(line.rstrip("\r\n"))
        if match:
            headings.append((index, len(match.group(1)), clean_heading(match.group(2))))

    h1_entries = [entry for entry in headings if entry[1] == 1]
    outputs: list[OutputSection] = []

    for h1_number, (h1_index, _, h1_title) in enumerate(h1_entries):
        h1_end = h1_entries[h1_number + 1][0] if h1_number + 1 < len(h1_entries) else len(lines)
        folder = destination / safe_name(h1_title)
        folder.mkdir(parents=True, exist_ok=False)

        local_headings = [entry for entry in headings if h1_index < entry[0] < h1_end]
        h2_entries = [entry for entry in local_headings if entry[1] == 2]
        listed: list[OutputSection] = []

        if h2_entries:
            preamble_end = h2_entries[0][0]
            preamble_lines = lines[h1_index:preamble_end]
            if has_code_fence(preamble_lines):
                relative = Path(safe_name(h1_basename(h1_title)) + ".md")
                listed.append(
                    OutputSection(
                        path=folder / relative,
                        heading=h1_basename(h1_title),
                        level=1,
                        start=h1_index + 1,
                        end=preamble_end,
                        content=source_slice(lines, h1_index, preamble_end),
                    )
                )

            for number, (section_index, level, section_title) in enumerate(h2_entries):
                section_end = h2_entries[number + 1][0] if number + 1 < len(h2_entries) else h1_end
                relative = Path(safe_name(section_title) + ".md")
                listed.append(
                    OutputSection(
                        path=folder / relative,
                        heading=section_title,
                        level=level,
                        start=section_index + 1,
                        end=section_end,
                        content=source_slice(lines, section_index, section_end),
                    )
                )
        else:
            h3_entries = [entry for entry in local_headings if entry[1] == 3]
            if h3_entries:
                preamble_end = h3_entries[0][0]
                relative = Path(safe_name(h1_basename(h1_title)) + ".md")
                listed.append(
                    OutputSection(
                        path=folder / relative,
                        heading=h1_basename(h1_title),
                        level=1,
                        start=h1_index + 1,
                        end=preamble_end,
                        content=source_slice(lines, h1_index, preamble_end),
                    )
                )
                for number, (section_index, level, section_title) in enumerate(h3_entries):
                    section_end = h3_entries[number + 1][0] if number + 1 < len(h3_entries) else h1_end
                    relative = Path(safe_name(section_title) + ".md")
                    listed.append(
                        OutputSection(
                            path=folder / relative,
                            heading=section_title,
                            level=level,
                            start=section_index + 1,
                            end=section_end,
                            content=source_slice(lines, section_index, section_end),
                        )
                    )
            else:
                relative = Path(safe_name(h1_basename(h1_title)) + ".md")
                listed.append(
                    OutputSection(
                        path=folder / relative,
                        heading=h1_basename(h1_title),
                        level=1,
                        start=h1_index + 1,
                        end=h1_end,
                        content=source_slice(lines, h1_index, h1_end),
                    )
                )

        for section in listed:
            provenance = f"<!-- Source: SheepEpimap.md lines {section.start}-{section.end} -->\n\n"
            write_text(section.path, provenance + section.content)
            outputs.append(section)

        intro_lines = lines[h1_index + 1 : h2_entries[0][0]] if h2_entries else []
        intro_text = "".join(intro_lines).strip() if intro_lines and not has_code_fence(intro_lines) else ""
        readme = [
            f"# {h1_title}\n",
            "\n",
            "This folder follows the heading hierarchy of the original `SheepEpimap.md`; code blocks, explanations, and lower-level headings under the same heading are kept in one Markdown file.\n",
        ]
        if intro_text:
            readme.extend(["\n", intro_text, "\n"])
        readme.extend(["\n", "## File order\n", "\n"])
        for section in listed:
            filename = section.path.name
            readme.append(f"- [{section.heading}](<{filename}>)\n")
        write_text(folder / "README.md", "".join(readme))

    return outputs


def write_source_map(root: Path, outputs: list[OutputSection]) -> None:
    path = root / "SOURCE_CODE_MAP.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["path", "source_lines", "heading_level", "heading", "format", "status"])
        for section in outputs:
            writer.writerow(
                [
                    "pipelines/" + section.path.relative_to(root).as_posix(),
                    f"{section.start}-{section.end}",
                    f"H{section.level}",
                    section.heading,
                    "markdown with original fenced code",
                    "requires configuration and data",
                ]
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--repository-root", required=True, type=Path)
    args = parser.parse_args()

    source = args.source.resolve()
    destination = args.destination.resolve()
    root = args.repository_root.resolve()
    if destination.exists():
        raise SystemExit(f"destination already exists: {destination}")
    if not source.is_file():
        raise SystemExit(f"source not found: {source}")
    if destination.parent != root.parent:
        raise SystemExit("staging destination must be a sibling of the repository root")

    destination.mkdir(parents=True)
    try:
        outputs = build_sections(source, destination)
    except Exception:
        shutil.rmtree(destination, ignore_errors=True)
        raise

    write_source_map(destination, outputs)
    print(f"Generated {len(outputs)} grouped code files in {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
