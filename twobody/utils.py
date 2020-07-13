# Third-party
from astropy.time import Time
from astropy.utils import check_broadcast
import numpy as np

__all__ = ['ArrayProcessor']


class ArrayProcessor:

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


def format_doc(docstring, *args, **kwargs):
    """
    A modified version of `astropy.utils.decorators.format_doc` that first
    applies ``.format()`` with ``__doc__``, then calls ``.format()`` with the
    arguments.

    Replaces the docstring of the decorated object and then formats it.
    The formatting works like :meth:`str.format` and if the decorated object
    already has a docstring this docstring can be included in the new
    documentation if you use the ``{__doc__}`` placeholder.
    Its primary use is for reusing a *long* docstring in multiple functions
    when it is the same or only slightly different between them.
    Parameters
    ----------
    docstring : str or object or None
        The docstring that will replace the docstring of the decorated
        object. If it is an object like a function or class it will
        take the docstring of this object. If it is a string it will use the
        string itself. One special case is if the string is ``None`` then
        it will use the decorated functions docstring and formats it.
    args :
        passed to :meth:`str.format`.
    kwargs :
        passed to :meth:`str.format`. If the function has a (not empty)
        docstring the original docstring is added to the kwargs with the
        keyword ``'__doc__'``.
    Raises
    ------
    ValueError
        If the ``docstring`` (or interpreted docstring if it was ``None``
        or not a string) is empty.
    IndexError, KeyError
        If a placeholder in the (interpreted) ``docstring`` was not filled. see
        :meth:`str.format` for more information.
    Notes
    -----
    Using this decorator allows, for example Sphinx, to parse the
    correct docstring.
    """
    def set_docstring(obj):
        if docstring is None:
            # None means: use the objects __doc__
            doc = obj.__doc__
            # Delete documentation in this case so we don't end up with
            # awkwardly self-inserted docs.
            obj.__doc__ = None
        elif isinstance(docstring, str):
            # String: use the string that was given
            doc = docstring
        else:
            # Something else: Use the __doc__ of this
            doc = docstring.__doc__

        if not doc:
            # In case the docstring is empty it's probably not what was wanted.
            raise ValueError('docstring must be a string or containing a '
                             'docstring that is not empty.')

        # If the original has a not-empty docstring append it to the format
        # kwargs.
        _doc = obj.__doc__ or ''
        doc = doc.format(__doc__=_doc)
        obj.__doc__ = doc.format(*args, **kwargs)
        return obj
    return set_docstring
