
#
# non-linear registration at full resolution
#

blockmatching -ref ref.mha -flo flo.mha -res res.mha -trsf-type vectorfield -res-trsf vector.trsf -elastic-sigma 2.5 -fluid-sigma 1.0

#
# creation of images of lower resolution
# note that resolutions can be different
#

applyTrsf flo.mha lower-flo.mha -iso 2.0 -resize -res-trsf lower-flo.trsf
applyTrsf ref.mha lower-ref.mha -iso 2.0 -resize -res-trsf lower-ref.trsf

#
# non-linear registration at lower resolution
#

blockmatching -ref lower-ref.mha -flo lower-flo.mha -res lower-res.mha -trsf-type vectorfield -res-trsf lower-vector.trsf -elastic-sigma 2.5 -fluid-sigma 1.0

#
# resample at full resolution the result at the lower resolution
# 
#

invTrsf lower-ref.trsf inv-lower-ref.trsf
applyTrsf lower-res.mha resampled-lower-res.mha -trsf inv-lower-ref.trsf -template ref.mha

#
# compute a transformation at full resolution by composition
#

composeTrsf resampled-lower-vector.trsf -trsfs lower-flo.trsf lower-vector.trsf inv-lower-ref.trsf -template ref.mha
applyTrsf flo.mha lower-trsf-res.mha -trsf resampled-lower-vector.trsf -template ref.mha


exit


for p in *.mha
do
    copy $p `basename $p .mha`.pgm -norma -o 1
    convert `basename $p .mha`.pgm `basename $p .mha`.png
    rm `basename $p .mha`.pgm
done
