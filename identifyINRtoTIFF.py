from definitions import *

from ImageHandling import imread, imsave
from lineage import timeNamed
im=imread(timeNamed(postsegment_files,begin))

import numpy as np

ids=np.unique(im)

print ids

