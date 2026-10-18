#!/usr/bin/env bash
# Acceptance check for the SLOW5/BLOW5 writer (internal/genomics/slow5, #15).
#
# Validates BLOW5/SLOW5 files with the reference tools:
#   slow5tools quickcheck, view, stats, index, get, merge
#   blue-crab s2p (BLOW5 -> pod5) and pod5 inspect, if installed
#
# Usage:
#   scripts/acceptance/slow5_validate.sh            # generate sample files, then validate
#   scripts/acceptance/slow5_validate.sh FILE|DIR...  # validate existing .blow5/.slow5 files
#
# Tools (override the paths with environment variables):
#   SLOW5TOOLS  slow5tools, built with zstd support   (required)
#               https://github.com/hasindu2008/slow5tools  (make zstd=1)
#   BLUECRAB    blue-crab, the SLOW5 <-> pod5 converter  (optional)
#               pip install blue-crab   (Python 3.8-3.12; provides `pod5` too)
#   POD5        pod5 command-line tool                   (optional)
#
# Note: slow5tools itself has no s2p command; BLOW5 -> pod5 conversion is
# blue-crab's `s2p`, from the same authors.
#
# Exit status: 0 all checks passed; 1 a check failed; 2 slow5tools missing;
# 3 slow5tools checks passed but the pod5 step was skipped (blue-crab missing).
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
slow5tools=${SLOW5TOOLS:-$(command -v slow5tools || true)}
bluecrab=${BLUECRAB:-$(command -v blue-crab || true)}
pod5tool=${POD5:-$(command -v pod5 || true)}

if [[ -z "$slow5tools" || ! -x "$slow5tools" ]]; then
	echo "slow5_validate: slow5tools not found. Install it (https://github.com/hasindu2008/slow5tools, build with 'make zstd=1')" >&2
	echo "                or set SLOW5TOOLS=/path/to/slow5tools." >&2
	exit 2
fi
if ! "$slow5tools" --version >/dev/null 2>&1; then
	echo "slow5_validate: $slow5tools does not run" >&2
	exit 2
fi

work=$(mktemp -d "${TMPDIR:-/tmp}/slow5_validate.XXXXXX")
trap 'rm -rf "$work"' EXIT

failures=0
pass() { printf 'PASS  %s\n' "$*"; }
fail() { printf 'FAIL  %s\n' "$*"; failures=$((failures + 1)); }
# run DESCRIPTION CMD... : run quietly, report, keep going.
run() {
	local desc=$1
	shift
	if "$@" >"$work/out" 2>"$work/err"; then
		pass "$desc"
	else
		fail "$desc"
		sed 's/^/      /' "$work/err" | tail -n 5
	fi
}

files=()
if [[ $# -eq 0 ]]; then
	sample="$work/sample"
	echo "Generating sample files in $sample"
	(cd "$repo_root" && PHANTOM_SLOW5_SAMPLE_DIR="$sample" go test -tags '!release' -count=1 \
		-run '^TestWriteSampleFiles$' ./internal/genomics/slow5/ >/dev/null)
	set -- "$sample"
fi
for arg in "$@"; do
	if [[ -d "$arg" ]]; then
		while IFS= read -r f; do files+=("$f"); done < <(find "$arg" -type f \( -name '*.blow5' -o -name '*.slow5' \) | sort)
	else
		files+=("$arg")
	fi
done
if [[ ${#files[@]} -eq 0 ]]; then
	echo "slow5_validate: no .blow5 or .slow5 files given" >&2
	exit 1
fi
echo "slow5tools: $("$slow5tools" --version 2>/dev/null | tail -n 1)  ($slow5tools)"
echo "Validating ${#files[@]} file(s)"

for f in "${files[@]}"; do
	name=$(basename "$f")
	copy="$work/files/$name"
	mkdir -p "$work/files"
	cp "$f" "$copy" # index files go next to the copy, not the input
	run "$name: quickcheck" "$slow5tools" quickcheck "$copy"
	run "$name: view (full decode)" "$slow5tools" view "$copy" --to slow5 -o "$work/view.slow5"
	if "$slow5tools" stats "$copy" >"$work/stats" 2>/dev/null; then
		n=$(awk -F'\t' '$1 == "number of records" {print $2}' "$work/stats")
		pass "$name: stats ($n records)"
	else
		fail "$name: stats"
		n=0
	fi
	run "$name: index" "$slow5tools" index "$copy"
	if [[ "${n:-0}" -gt 0 ]]; then
		first=$(awk -F'\t' '!/^[#@]/ {print $1; exit}' "$work/view.slow5")
		last=$(awk -F'\t' '!/^[#@]/ {id = $1} END {print id}' "$work/view.slow5")
		run "$name: get $first $last" "$slow5tools" get "$copy" "$first" "$last" --to slow5 -o "$work/get.slow5"
		got=$(grep -vc '^[#@]' "$work/get.slow5" || true)
		if [[ "$got" == 2 ]]; then
			pass "$name: get returned 2 records"
		else
			fail "$name: get returned $got records, want 2"
		fi
	fi
	# Re-encode to slow5tools' default BLOW5 (zlib + svb-zd) and back to
	# ASCII: the records must be unchanged.
	run "$name: re-encode zlib/svb-zd" "$slow5tools" view "$copy" -c zlib -s svb-zd -o "$work/re.blow5"
	if "$slow5tools" view "$work/re.blow5" --to slow5 -o "$work/re.slow5" 2>/dev/null &&
		cmp -s <(grep -v '^#slow5_version' "$work/view.slow5") <(grep -v '^#slow5_version' "$work/re.slow5"); then
		pass "$name: re-encoded records identical"
	else
		fail "$name: re-encoded records differ"
	fi
done

# Chunk directories must merge into one valid file.
chunkdirs=$(for f in "${files[@]}"; do [[ "$f" == *.blow5 ]] && dirname "$f"; done | sort | uniq -d || true)
for d in $chunkdirs; do
	run "merge $(basename "$d")/" "$slow5tools" merge "$d" -o "$work/merged.blow5"
	[[ -f "$work/merged.blow5" ]] && run "merged file quickcheck" "$slow5tools" quickcheck "$work/merged.blow5"
done

pod5_skipped=0
if [[ -n "$bluecrab" && -x "$bluecrab" ]]; then
	for f in "${files[@]}"; do
		[[ "$f" == *.blow5 ]] || continue
		name=$(basename "$f" .blow5)
		out="$work/pod5/$name.pod5"
		mkdir -p "$work/pod5"
		src=$f
		# The pyslow5 wheels on PyPI are built without zstd, so blue-crab
		# cannot open zstd BLOW5 (byte 9 of the header = 2). Re-encode those
		# to zlib with slow5tools first; the records are unchanged.
		if [[ $(od -An -tu1 -j9 -N1 "$f" | tr -d ' ') == 2 ]]; then
			src="$work/pod5/$name.zlib.blow5"
			"$slow5tools" view "$f" -c zlib -s svb-zd -o "$src" 2>/dev/null
			echo "NOTE  $name: zstd record compression; converting a zlib copy (blue-crab's pyslow5 lacks zstd)"
		fi
		run "$name: blue-crab s2p -> pod5" "$bluecrab" s2p "$src" -o "$out"
		if [[ -n "$pod5tool" && -x "$pod5tool" && -f "$out" ]]; then
			run "$name: pod5 inspect summary" "$pod5tool" inspect summary "$out"
		fi
		# Round trip pod5 -> BLOW5 must reproduce the samples.
		if [[ -f "$out" ]]; then
			run "$name: blue-crab p2s" "$bluecrab" p2s "$out" -o "$work/pod5/$name.back.blow5"
			if "$slow5tools" view "$work/pod5/$name.back.blow5" --to slow5 2>/dev/null | awk -F'\t' '!/^[#@]/ {print $1 "\t" $7 "\t" $8}' | sort >"$work/back.sig" &&
				"$slow5tools" view "$f" --to slow5 2>/dev/null | awk -F'\t' '!/^[#@]/ {print $1 "\t" $7 "\t" $8}' | sort >"$work/orig.sig" &&
				cmp -s "$work/back.sig" "$work/orig.sig"; then
				pass "$name: pod5 round trip preserves every read ID and sample"
			else
				fail "$name: pod5 round trip changed reads or samples"
			fi
		fi
	done
else
	echo "SKIP  pod5 conversion: blue-crab not found (pip install blue-crab, or set BLUECRAB=...)"
	pod5_skipped=1
fi

echo
if [[ $failures -gt 0 ]]; then
	echo "slow5_validate: $failures check(s) FAILED"
	exit 1
fi
if [[ $pod5_skipped -eq 1 ]]; then
	echo "slow5_validate: slow5tools checks passed; pod5 step skipped"
	exit 3
fi
echo "slow5_validate: all checks passed"
