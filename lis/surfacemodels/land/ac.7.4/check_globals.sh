#!/usr/bin/env bash
#-----------------------------------------------------------------------
# List the module-level variables of the AquaCrop modules (ac_modules)
# that ac74_main.F90 does not restore per tile, with SetX(...) or X = ...
#
# AquaCrop keeps its state in module variables, which all tiles share,
# so every one that carries state from one day to the next must be
# stored in AC74_struc and restored before AdvanceOneTimeStep. Run it
# after each sync; anything listed must be restored or added to SKIP
# below, with the reason.
#
# usage: ./check_globals.sh
#
# 19 Oct 2026: Louise Busschaert; initial implementation
#-----------------------------------------------------------------------
set -euo pipefail

dir=$(cd "$(dirname "$0")" && pwd)
main="$dir/ac74_main.F90"

# file names, descriptions, paths, file units: no state
NOSTATE='file|filefull|description|pathname|outputname|_iostat|_filename|_exists|^f[a-z0-9]+$|warning_log|handler|^message$'

# reviewed, not per-tile state (name: reason)
SKIP="
cuttings: restored with SetManagement_Cuttings
rootzonewc: restored with the SetRootZoneWC_* of its components
etodataset: 10-daily and monthly climate records, not with daily forcing
raindataset: 10-daily and monthly climate records, not with daily forcing
tmindataset: 10-daily and monthly climate records, not with daily forcing
tmaxdataset: 10-daily and monthly climate records, not with daily forcing
globalirriecw: set from the irrigation file, the same for all tiles
lastirridap: only for the irrigation output file, not written in LIS
externalrun: same for all tiles
repeattoday: standalone run loop only
kspolcold: reset every day
kspolheat: reset every day
ksaeration: reset every day
management_brel: set before each run initialization (AC74_set_run_input)
tmintnxreference365daysrun: rebuilt from the 12 months in main
tmaxtnxreference365daysrun: rebuilt from the 12 months in main
tmintnxreference12monthsrun: set from AC74_struc in main
tmaxtnxreference12monthsrun: set from AC74_struc in main
tmincropreferencerun: only used when a run is initialized
tmaxcropreferencerun: only used when a run is initialized
tmincropreference: only used when a run is initialized
tmaxcropreference: only used when a run is initialized
tmintnxreference365days: only used when a run is initialized
tmaxtnxreference365days: only used when a run is initialized
tnxreferenceyear: set before a run is initialized
tmin: forcing, set every day
tmax: forcing, set every day
rain: forcing, set every day
eto: forcing, set every day
"

# module-level variables: declarations before 'contains', outside types
vars=$(for f in "$dir"/ac_modules/ac74_*.[fF]90; do
    awk -v F="$(basename "$f")" '
        tolower($0) ~ /^\s*contains\s*$/ { exit }
        tolower($0) ~ /^\s*type\s*(,|::|[a-z])/ && tolower($0) !~ /^\s*type\s*\(/ { t = 1 }
        t && tolower($0) ~ /^\s*end\s*type/ { t = 0; next }
        !t && /::/ && tolower($0) !~ /parameter|public|private|interface|procedure|^\s*use/ {
            sub(/!.*/, ""); split($0, a, "::")
            gsub(/\([^)]*\)/, "", a[2]); gsub(/=[^,]*/, "", a[2]); gsub(/ /, "", a[2])
            m = split(a[2], v, ",")
            for (i = 1; i <= m; i++) if (v[i] ~ /^[A-Za-z]/) print F, tolower(v[i])
        }' "$f"
done | sort -u -k2)

status=0
while read -r f v; do
    [[ $v =~ $NOSTATE ]] && continue
    grep -qiE "^${v}:" <<< "$SKIP" && continue
    if ! grep -qiE "\bSet${v}\b|^\s*${v}\s*=" "$main"; then
        echo "not restored: $v ($f)"
        status=1
    fi
done <<< "$vars"
[[ $status == 0 ]] && echo "OK: every module-level state is restored in ac74_main.F90"
exit $status
