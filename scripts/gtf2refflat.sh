
GTF_IN=$1
REFFLAT_OUT=$2

tools/gtfToGenePred \
    -genePredExt \
    $GTF_IN refFlat.tmp.txt

paste <(cut -f 12 refFlat.tmp.txt) \
    <(cut -f 1-10 refFlat.tmp.txt) > $REFFLAT_OUT
