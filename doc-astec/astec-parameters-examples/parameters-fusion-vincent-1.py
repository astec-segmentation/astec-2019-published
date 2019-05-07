
# path to the embryo data

PATH_EMBRYO = '/Users/greg/COLLABORATIONS/DIGEM/DATA/TEST/180927-Vincent-Nodal-St8'

# Embryo Name

EN = '180927-Vincent-Nodal-St8'

# time points to be processed
# begin: first time point
# end: last time point
# delta: Delta between two time points (if one does not want
#        to deal with every single time point) (default = 1)

begin = 0
end = 0
delta = 1

# output image format
#
# default_image_suffix is for all images (including auxiliary ones)

default_image_suffix = 'mha'

# raw_ori: image orientation, value is in {'right','left'}
#   if im2 angle - im1 angle < 0 => raw_ori = 'right'
# raw_mirrors: value is in {True, False}
# 	- the standard value of this parameter is False
# 	- in case of axial symmetry between left
# 	  and right cameras, then set to True
# raw_resolution: acquisition voxel size
# raw_delay: increment to to be added to the time values
#   (values in range [begin,end]) when generating the
#   fused image names (default is 0)

raw_ori = 'left'
raw_mirrors = False
raw_resolution = (.195, .195, 1.)
raw_delay = 0

# Isotropic resolution of the final fused image

target_resolution = .3

