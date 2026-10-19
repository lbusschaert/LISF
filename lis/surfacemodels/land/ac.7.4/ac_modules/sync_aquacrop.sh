#!/usr/bin/env bash
#-----------------------------------------------------------------------
# Copy the AquaCrop Fortran source of a given commit into this
# directory, prefixing every AquaCrop module (ac_X -> ac74_X) and file
# (X.f90 -> ac74_X.f90), so that this version can be compiled in the
# same LIS executable as other AquaCrop versions (e.g., ac.7.2).
#
# The renaming is the only change made to the AquaCrop source; it is
# undone and compared to the commit with --check.
#
# usage:
#   ./sync_aquacrop.sh <AquaCrop repo> [<commit>] [<subdir>]
#   ./sync_aquacrop.sh --check
#
#   <commit> defaults to HEAD, <subdir> to src. The commit is taken
#   from git (git archive), uncommitted changes are never copied.
#   --check re-reads repo, commit and subdir from VERSION.txt.
#
# 19 Oct 2026: Louise Busschaert; initial implementation
#-----------------------------------------------------------------------
set -euo pipefail

prefix=ac74
dst=$(cd "$(dirname "$0")" && pwd)
wrappers=$(dirname "$dst")
# standalone main program, not used in LIS
exclude="aquacrop.F90"

if [[ "${1:-}" == "--check" ]]; then
    mode=check
    repo=$(sed -n 's/^repo: *//p' "$dst/VERSION.txt")
    commit=$(sed -n 's/^commit: *//p' "$dst/VERSION.txt")
    subdir=$(sed -n 's/^subdir: *//p' "$dst/VERSION.txt")
else
    mode=sync
    repo=${1:?usage: $0 <AquaCrop repo> [<commit>] [<subdir>]}
    commit=$(git -C "$repo" rev-parse "${2:-HEAD}^{commit}")
    subdir=${3:-src}
fi
repo=$(cd "$repo" && git rev-parse --show-toplevel)

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
git -C "$repo" archive "$commit" "$subdir" | tar -x -C "$tmp"
src="$tmp/$subdir"

sources=$(cd "$src" && ls *.f90 *.F90 2>/dev/null | grep -vxF "$exclude")
modules=$(cd "$src" && grep -hiE '^\s*module\s+ac_[a-z0-9_]+\s*(!.*)?$' $sources \
          | awk '{print tolower($2)}' | sort -u)

# case-preserving, whole-word: ac_global -> ac74_global, AC_GLOBAL -> AC74_GLOBAL
fwd=() ; rev=()
for m in $modules; do
    fwd+=(-e "s/\b(ac)_(${m#ac_})\b/\1${prefix#ac}_\2/Ig")
    rev+=(-e "s/\b(ac)${prefix#ac}_(${m#ac_})\b/\1_\2/Ig")
done

if [[ $mode == sync ]]; then
    if grep -qiwE "${prefix}_[a-z0-9_]+" $(printf "$src/%s " $sources); then
        echo "ERROR: source already contains ${prefix}_ names, renaming is not reversible" >&2
        exit 1
    fi
    rm -f "$dst"/${prefix}_*.f90 "$dst"/${prefix}_*.F90
    for f in $sources; do
        sed -E "${fwd[@]}" "$src/$f" > "$dst/${prefix}_$f"
    done
    for f in LICENSE AUTHORS.md; do
        for d in "$src" "$tmp"; do
            [[ -f "$d/$f" ]] && { cp "$d/$f" "$dst/$f"; break; }
        done
        [[ -f "$dst/$f" ]] || git -C "$repo" show "$commit:$f" > "$dst/$f" 2>/dev/null || true
    done
    # the LIS wrappers use the same modules
    sed -i -E "${fwd[@]}" "$wrappers"/*.F90
    cat > "$dst/VERSION.txt" <<EOF
AquaCrop source synced with sync_aquacrop.sh, do not edit the ${prefix}_* files by hand.
describe: $(git -C "$repo" describe --tags --always "$commit")
commit: $commit
subdir: $subdir
repo: $repo
modules: $(echo $modules)
EOF
    echo "synced $(echo $sources | wc -w) files from $commit:$subdir"
else
    status=0
    for f in $sources; do
        if ! diff -q <(sed -E "${rev[@]}" "$dst/${prefix}_$f") "$src/$f" >/dev/null; then
            echo "DIFFERS: ${prefix}_$f"; status=1
        fi
    done
    for f in "$dst"/${prefix}_*.[fF]90; do
        b=$(basename "$f"); b=${b#${prefix}_}
        grep -qxF "$b" <<< "$sources" || { echo "EXTRA:   ${prefix}_$b"; status=1; }
    done
    [[ $status == 0 ]] && echo "OK: identical to $commit:$subdir up to the ${prefix}_ renaming"
    exit $status
fi
