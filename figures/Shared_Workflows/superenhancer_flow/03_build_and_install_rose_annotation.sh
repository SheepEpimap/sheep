#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=01_config.sh
source "${SCRIPT_DIR}/01_config.sh"

require_cmd awk
require_cmd python3
require_file "${GTF_FILE}"
[[ -d "${ROSE_HOME}" ]] || die "ROSE_HOME does not exist: ${ROSE_HOME}"

local_annotation="${WORK_ROOT}/annotation/${ROSE_ANNOTATION_NAME}"
mkdir -p "$(dirname "${local_annotation}")"

log "Building ROSE UCSC-style annotation from: ${GTF_FILE}"
if [[ "${GTF_FILE}" == *.gz ]]; then
  reader=(gzip -dc "${GTF_FILE}")
else
  reader=(cat "${GTF_FILE}")
fi

{
  printf '#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\n'
  "${reader[@]}" | awk -F '\t' 'BEGIN{OFS="\t"}
    $3=="gene" {
      id="";
      n=split($9,a,";");
      for(i=1;i<=n;i++) {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", a[i]);
        if(a[i] ~ /^gene_id[[:space:]]+/) {
          sub(/^gene_id[[:space:]]+/, "", a[i]);
          gsub(/"/, "", a[i]);
          id=a[i];
          break;
        }
      }
      if(id=="") next;
      chr=$1;
      if(chr !~ /^chr/) chr="chr" chr;
      # Coordinate handling intentionally follows the original workflow to keep
      # historical results comparable. Change txStart/cdsStart to $4-1 only if
      # the local ROSE installation explicitly requires strict 0-based starts.
      print 0,id,chr,$7,$4,$5,$4,$5,".",".",".",".",id
    }'
} > "${local_annotation}"

require_file "${local_annotation}"

# ROSE releases differ in whether genomeDict resolves annotation relative to
# ROSE_HOME or ROSE_HOME/bin. Install the same generated file in both locations.
for annotation_dir in "${ROSE_HOME}/annotation" "${ROSE_HOME}/bin/annotation"; do
  mkdir -p "${annotation_dir}" || die "Cannot create ROSE annotation directory: ${annotation_dir}"
  cp -f "${local_annotation}" "${annotation_dir}/${ROSE_ANNOTATION_NAME}"
done

log "Registering ${ROSE_GENOME_KEY} in ROSE_main.py and ROSE_geneMapper.py"
python3 - "${ROSE_HOME}" "${ROSE_GENOME_KEY}" "${ROSE_ANNOTATION_NAME}" <<'PY'
from pathlib import Path
import shutil
import sys

rose_home = Path(sys.argv[1])
genome_key = sys.argv[2]
annotation_name = sys.argv[3]

for name in ("ROSE_main.py", "ROSE_geneMapper.py"):
    path = rose_home / "bin" / name
    if not path.exists():
        raise SystemExit(f"Missing ROSE script: {path}")
    text = path.read_text()
    if f"'{genome_key}'" in text or f'"{genome_key}"' in text:
        print(f"Already registered: {path}")
        continue
    marker = "genomeDict = {"
    if marker not in text:
        raise SystemExit(f"Could not locate genomeDict in {path}")
    backup = path.with_suffix(path.suffix + ".pre_oviari2.bak")
    if not backup.exists():
        shutil.copy2(path, backup)
    entry = f"\n        '{genome_key}':'%s/annotation/{annotation_name}' % (cwd),"
    path.write_text(text.replace(marker, marker + entry, 1))
    print(f"Patched: {path}")
PY

log "ROSE annotation installed: ${local_annotation}"
