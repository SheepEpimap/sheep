#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
failed=0

for script in "${SCRIPT_DIR}"/*.sh; do
    if bash -n "${script}"; then
        printf 'OK bash -n: %s\n' "$(basename "${script}")"
    else
        failed=1
    fi
done

if grep -RIl $'\r' "${ROOT}/scripts" "${ROOT}/envs" \
    --include='*.sh' --include='*.R' --include='*.md' | grep -q .; then
    echo 'ERROR: CRLF line endings detected'
    failed=1
fi

grep -q 'ATAC_QVALUE="${ATAC_QVALUE:-0.01}"' "${SCRIPT_DIR}/00_config.sh" || {
    echo 'ERROR: ATAC q-value rule is missing'
    failed=1
}

grep -q 'OTHER_PEAK_QVALUE="${OTHER_PEAK_QVALUE:-0.05}"' "${SCRIPT_DIR}/00_config.sh" || {
    echo 'ERROR: non-ATAC q-value rule is missing'
    failed=1
}

if grep -RIn '/vol2/' "${SCRIPT_DIR}" --include='*.sh' | grep -v '/00_config.sh:'; then
    echo 'ERROR: hard-coded /vol2 path found outside 00_config.sh'
    failed=1
fi

for required in README.md CODE_AVAILABILITY.md CONTRIBUTING.md .gitignore .gitattributes; do
    if [[ ! -s "${ROOT}/${required}" ]]; then
        echo "ERROR: missing repository file: ${required}"
        failed=1
    fi
done

if command -v python3 >/dev/null 2>&1; then
    if ! python3 - "${ROOT}/README.md" "${ROOT}/resources/README.md" <<'PY'
import re
import sys
from pathlib import Path

pattern = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff]")
failed = False
for filename in sys.argv[1:]:
    text = Path(filename).read_text(encoding="utf-8")
    if pattern.search(text):
        print(f"ERROR: CJK characters detected in {filename}")
        failed = True
raise SystemExit(1 if failed else 0)
PY
    then
        failed=1
    fi
fi

if (( failed == 0 )); then
    echo 'Static package checks passed.'
fi
exit "${failed}"
