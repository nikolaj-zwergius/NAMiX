SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")
python3 "$SCRIPTPATH/NAMIX_bat.py" $*
