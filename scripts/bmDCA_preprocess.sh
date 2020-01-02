#! /bin/bash
set -eu

OPTIONS=i:d:rh
LONGOPTIONS=input:,destination:,reweight,help
PARSED=$(getopt -o ${OPTIONS} --long ${LONGOPTIONS} -n "$0" -- "$@")
eval set -- "$PARSED"

while [ $# -ge 1 ]; do
  case "$1" in
    -i|--input)
      infile="${2}"
      shift 2
      ;;
    -d|--destination)
      destdir="${2}"
      shift 2
      ;;
    -r|--reweight)
      reweight=true
      shift
      ;;
    -h|--help)
      print_usage
      exit
      ;;
    --)
      shift
      break
      ;;
    *)
      echo "ERROR"
      exit 3
      ;;
  esac
done

mkdir -p ${destdir}
if [[ -v reweight ]]; then
  preprocess -i "${infile}" -d "${destdir}" -r
else
  preprocess -i "${infile}" -d "${destdir}"
fi
