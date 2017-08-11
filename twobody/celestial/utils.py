# Third-party
from astropy.utils.misc import check_broadcast
import numpy as np

__all__ = ['ArrayProcessor']

class ArrayProcessor(object):

    def __init__(self, *arrs):
        self.arrs = [np.array(arr) for arr in arrs]

    def prepare_arrays(self):
        """Make sure input arrays are all C-contiguous and have same shape."""
        self.max_shape = None
        for arr in self.arrs:
            if self.max_shape is None:
                self.max_shape = arr.shape
            elif arr.shape > self.max_shape:
                self.max_shape = arr.shape

        orig_shapes = []
        arrs_1d = []
        for arr in self.arrs:
            orig_shapes.append(arr.shape)
            arr = np.broadcast_to(arr, self.max_shape).ravel()
            arrs_1d.append(np.ascontiguousarray(arr.astype(np.float64)))

        if not check_broadcast(orig_shapes):
            raise ValueError("Shapes are not broadcastable: {0}"
                             .format(orig_shapes))

        return arrs_1d

    def prepare_result(self, res):
        return res.reshape(self.max_shape)
