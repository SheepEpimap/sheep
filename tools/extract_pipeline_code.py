#!/usr/bin/env python3
"""Export SheepEpimap Markdown code fences as numbered native code files."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from dataclasses import dataclass
from pathlib import Path


HEADING_RE = re.compile(r"^(#{1,6})\s+(.+?)\s*$")
FENCE_RE = re.compile(r"^\s*(```|~~~)(.*)$")
INVALID_WINDOWS_CHARS = re.compile(r'[<>:"/\\|?*]')


@dataclass(frozen=True)
class CodeBlock:
    number: int
    language_label: str
    h1: str
    h2: str
    h3: str
    start_line: int
    end_line: int
    fence_start_line: int
    fence_end_line: int
    content: str


# These names describe each source block without fragmenting it further.
DESCRIPTION_BY_BLOCK = {
    1: "ATAC-Seq_data_analysis",
    2: "CUT_Tag_data_analysis",
    3: "call_peak",
    4: "SPP",
    5: "FRIP",
    6: "Heatmap_Matrix_all_score_CUT_Tag_ATAC_normalization",
    7: "Heatmap_Matrix_all_score_RNA-Seq_normalization",
    8: "Heatmap_Matrix_all_score_summary_and_correlation",
    9: "Signal_around_TSS",
    10: "PCA_for_each_assay_generate_input",
    11: "PCA_for_each_assay_prepare_table",
    12: "PCA_for_each_assay_plot",
    13: "ChromImpute_ChromHMM",
    14: "TSR_core_workflow",
    15: "TSR_summary",
    16: "TSR_link_to_target_gene",
    17: "Promoter_target_genes",
    18: "Enhancer_target_gene_summary",
    19: "Enhancer_assign_cCRE_IDs",
    20: "GREAT_liftover_regions",
    21: "TSR_GO_enrichment",
    22: "Human_Phenotype_enrichment",
    23: "Mouse_Phenotype_enrichment",
    24: "Mouse_Phenotype_heatmap_plot",
    25: "TSR_MOTIF_HOMER",
    26: "TSR_MOTIF_MEME_FIMO",
    27: "TSE_scores",
    28: "Network_construction",
    29: "Network_statistics",
    30: "Selection_signatures_inputs",
    31: "Selection_signatures_CDS_overlap",
    32: "Top_10_percent_GAT_environment",
    33: "Top_10_percent_fold_enrichment",
    34: "Top_10_percent_FST_plot",
    35: "Calculate_point_FST_neweurope_oldeurope",
    36: "Calculate_point_FST_eur_CEA",
    37: "Calculate_point_FST_merge",
}


FORMAT_BY_BLOCK = {
    1: ("snakemake", ".smk"),
    2: ("snakemake", ".smk"),
    3: ("snakemake", ".smk"),
    6: ("snakemake", ".smk"),
    7: ("snakemake", ".smk"),
    12: ("R", ".R"),
    24: ("R", ".R"),
    34: ("python", ".py"),
}


def clean_heading(raw: str) -> str:
    value = raw.strip()
    value = re.sub(r"^\*\*(.+)\*\*$", r"\1", value)
    return value.strip()


def safe_folder_name(raw: str) -> str:
    value = INVALID_WINDOWS_CHARS.sub("_", clean_heading(raw))
    value = re.sub(r"_+", "_", value).rstrip(" .")
    return value or "untitled"


def safe_file_stem(raw: str) -> str:
    value = INVALID_WINDOWS_CHARS.sub("_", clean_heading(raw))
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_.")
    return value or "untitled"


def is_comment_only(content: str) -> bool:
    meaningful = [
        line.strip()
        for line in content.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    return not meaningful


def parse_blocks(path: Path) -> list[CodeBlock]:
    lines = path.read_text(encoding="utf-8-sig").splitlines(keepends=True)
    blocks: list[CodeBlock] = []
    h1 = ""
    h2 = ""
    h3 = ""
    open_index: int | None = None
    marker = ""
    language = ""

    for index, line in enumerate(lines):
        if open_index is None:
            heading_match = HEADING_RE.match(line.rstrip("\r\n"))
            if heading_match:
                level = len(heading_match.group(1))
                title = clean_heading(heading_match.group(2))
                if level == 1:
                    h1, h2, h3 = title, "", ""
                elif level == 2:
                    h2, h3 = title, ""
                elif level == 3:
                    h3 = title
            fence_match = FENCE_RE.match(line.rstrip("\r\n"))
            if fence_match:
                open_index = index
                marker = fence_match.group(1)
                language = fence_match.group(2).strip()
            continue

        stripped = line.lstrip()
        if stripped.startswith(marker) and not stripped[len(marker) :].strip():
            content = "".join(lines[open_index + 1 : index])
            blocks.append(
                CodeBlock(
                    number=len(blocks) + 1,
                    language_label=language,
                    h1=h1,
                    h2=h2,
                    h3=h3,
                    start_line=open_index + 2,
                    end_line=index,
                    fence_start_line=open_index + 1,
                    fence_end_line=index + 1,
                    content=content,
                )
            )
            open_index = None
            marker = ""
            language = ""

    if open_index is not None:
        raise ValueError(f"Unclosed code fence starting at line {open_index + 1}: {path}")
    return blocks


def heading_for(block: CodeBlock) -> tuple[str, str]:
    if block.h3:
        return "H3", block.h3
    if block.h2:
        return "H2", block.h2
    return "H1", block.h1


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8", newline="\n")


def build_repository(
    source: Path,
    original_source: Path,
    destination: Path,
) -> tuple[int, int]:
    blocks = parse_blocks(source)
    original_blocks = parse_blocks(original_source)
    if len(blocks) != len(original_blocks):
        raise ValueError(
            f"Translated source has {len(blocks)} blocks; original has {len(original_blocks)}"
        )
    for translated, original in zip(blocks, original_blocks):
        translated_lines = (
            translated.start_line,
            translated.end_line,
            translated.fence_start_line,
            translated.fence_end_line,
        )
        original_lines = (
            original.start_line,
            original.end_line,
            original.fence_start_line,
            original.fence_end_line,
        )
        if translated_lines != original_lines:
            raise ValueError(f"Source line structure differs at block {translated.number}")

    export_blocks = [block for block in blocks if not is_comment_only(block.content)]
    skipped_blocks = [block for block in blocks if is_comment_only(block.content)]
    missing_names = [block.number for block in export_blocks if block.number not in DESCRIPTION_BY_BLOCK]
    if missing_names:
        raise ValueError(f"Missing descriptions for source blocks: {missing_names}")

    stage_order: list[str] = []
    stage_blocks: dict[str, list[CodeBlock]] = {}
    for block in export_blocks:
        if block.h1 not in stage_blocks:
            stage_order.append(block.h1)
            stage_blocks[block.h1] = []
        stage_blocks[block.h1].append(block)

    rows: list[dict[str, str]] = []
    root_readme = [
        "# Numbered pipeline code\n",
        "\n",
        "Each directory follows an H1 section in `SheepEpimap.md`. Native code files are numbered in source order within each directory, and spaces in filenames are replaced with underscores. Each Markdown code fence is kept intact in one file.\n",
        "\n",
        "## Directories\n",
        "\n",
    ]

    for h1 in stage_order:
        folder_name = safe_folder_name(h1)
        folder = destination / folder_name
        folder.mkdir(parents=True, exist_ok=False)
        folder_blocks = stage_blocks[h1]
        root_readme.append(f"- [{h1}](<{folder_name}/README.md>) — {len(folder_blocks)} code files\n")
        stage_readme = [
            f"# {h1}\n",
            "\n",
            "Files are numbered in the original source order. Each file contains one complete fenced code block from the English-translated source archive.\n",
            "\n",
            "| Order | File | Source heading | Original lines | Type |\n",
            "|---:|---|---|---:|---|\n",
        ]
        for order, block in enumerate(folder_blocks, start=1):
            language, extension = FORMAT_BY_BLOCK.get(block.number, ("bash", ".sh"))
            stem = safe_file_stem(DESCRIPTION_BY_BLOCK[block.number])
            filename = f"{order:02d}_{stem}{extension}"
            output_path = folder / filename
            write_text(output_path, block.content)
            heading_level, heading = heading_for(block)
            stage_readme.append(
                f"| {order:02d} | [{filename}]({filename}) | {heading.replace('|', '\\|')} | {block.start_line}-{block.end_line} | {language} |\n"
            )
            rows.append(
                {
                    "path": "pipelines/" + output_path.relative_to(destination).as_posix(),
                    "source_block": str(block.number),
                    "source_lines": f"{block.start_line}-{block.end_line}",
                    "source_fence_lines": f"{block.fence_start_line}-{block.fence_end_line}",
                    "heading_level": heading_level,
                    "heading": heading,
                    "language": language,
                    "format": "native code file",
                    "status": "requires configuration and data",
                }
            )
        write_text(folder / "README.md", "".join(stage_readme))

    root_readme.extend(
        [
            "\n",
            "## Extraction note\n",
            "\n",
            f"Exported {len(export_blocks)} non-placeholder code blocks. "
            f"Skipped {len(skipped_blocks)} comment-only placeholder fences "
            f"(source blocks {', '.join(str(block.number) for block in skipped_blocks)}).\n",
        ]
    )
    write_text(destination / "README.md", "".join(root_readme))

    map_path = destination / "SOURCE_CODE_MAP.tsv"
    fieldnames = [
        "path",
        "source_block",
        "source_lines",
        "source_fence_lines",
        "heading_level",
        "heading",
        "language",
        "format",
        "status",
    ]
    with map_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(export_blocks), len(skipped_blocks)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--original-source", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--repository-root", required=True, type=Path)
    args = parser.parse_args()

    source = args.source.resolve()
    original_source = args.original_source.resolve()
    destination = args.destination.resolve()
    root = args.repository_root.resolve()
    if destination.exists():
        raise SystemExit(f"Destination already exists: {destination}")
    if not source.is_file() or not original_source.is_file():
        raise SystemExit("Source or original source is missing")
    if destination.parent != root.parent:
        raise SystemExit("Staging destination must be a sibling of the repository root")

    destination.mkdir(parents=True)
    try:
        exported, skipped = build_repository(source, original_source, destination)
    except Exception:
        shutil.rmtree(destination, ignore_errors=True)
        raise
    print(f"Exported {exported} code blocks; skipped {skipped} placeholders")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
