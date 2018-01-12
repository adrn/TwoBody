# Third-party
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import (frame_transform_graph, CoordinateAttribute,
                                 FunctionTransform, AffineTransform, ICRS)
from astropy.coordinates.matrix_utilities import (rotation_matrix,
                                                  matrix_product,
                                                  matrix_transpose)

__all__ = ['ReferencePlaneFrame']

_cache = {}


def _make_cls(framecls):
    """Analogous to
    astropy.coordinates.builtin_frames.skyoffset._make_skyoffset_cls -- see
    notes in there.
    """

    if framecls in _cache:
        return _cache[framecls]

    # the class of a class object is the metaclass
    framemeta = framecls.__class__

    class ReferencePlaneMeta(framemeta):
        """
        This metaclass renames the class to be "ReferencePlane<framecls>" and
        also adjusts the frame specific representation info so that spherical
        names are always "lon" and "lat" (instead of, e.g., "ra" and "dec").
        """

        def __new__(cls, name, bases, members):
            # Only 'origin' is needed here, to set the origin frame properly.
            members['origin'] = CoordinateAttribute(frame=framecls,
                                                    default=None)

            # This has to be done because FrameMeta will set these attributes
            # to the defaults from BaseCoordinateFrame when it creates the base
            # SkyOffsetFrame class initially.
            members['_default_representation'] = framecls._default_representation
            members['_default_differential'] = framecls._default_differential

            newname = name[:-5] if name.endswith('Frame') else name
            newname += framecls.__name__

            return super().__new__(cls, newname, bases, members)

    _ReferencePlaneFramecls = ReferencePlaneMeta('ReferencePlaneFrame',
        (ReferencePlaneFrame, framecls),
        {'__doc__': ReferencePlaneFrame.__doc__})

    @frame_transform_graph.transform(FunctionTransform,
                                     _ReferencePlaneFramecls,
                                     _ReferencePlaneFramecls)
    def referenceplane_to_referenceplane(from_coord, to_frame):
        """Transform between two reference plane frames."""

        if not from_coord.origin.has_data or not to_frame.origin.has_data:
            raise ValueError("Reference plane origin frame object must have "
                             "data, i.e. it must have a position specified.")

        # This transform goes through the parent frames on each side.
        # from_frame -> from_frame.origin -> to_frame.origin -> to_frame
        intermediate_from = from_coord.transform_to(from_coord.origin)
        intermediate_to = intermediate_from.transform_to(to_frame.origin)
        return intermediate_to.transform_to(to_frame)

    @frame_transform_graph.transform(AffineTransform, _ReferencePlaneFramecls,
                                     framecls)
    def referenceplane_to_coord(reference_plane_coord, to_frame):
        """Convert a reference plane coordinate to the original system"""

        if (reference_plane_coord.origin is None or
                not reference_plane_coord.origin.has_data):
            raise ValueError("Reference plane origin frame object must have "
                             "data, i.e. it must have a position specified.")

        origin = reference_plane_coord.origin.spherical
        if origin.distance.unit == u.dimensionless_unscaled:
            dist = 0. * u.pc
        else:
            dist = origin.distance
        sun = coord.CartesianRepresentation([0, 0, 1] * dist)

        mat1 = rotation_matrix(90*u.deg+origin.lat, 'y')
        mat2 = rotation_matrix(-origin.lon, 'z')
        M = matrix_product(mat2, mat1)

        offset = (-sun).transform(M)
        return M, offset

    @frame_transform_graph.transform(AffineTransform, framecls,
                                     _ReferencePlaneFramecls)
    def coord_to_referenceplane(from_coord, to_reference_plane):
        """Convert a sky coordinate to a reference-plane frame."""

        if (to_reference_plane.origin is None or
                not to_reference_plane.origin.has_data):
            raise ValueError("Reference plane origin frame object must have "
                             "data, i.e. it must have a position specified.")

        M, offset = referenceplane_to_coord(to_reference_plane, from_coord)
        M = matrix_transpose(M)
        offset = (-offset).transform(M)
        return M, offset

    _cache[framecls] = _ReferencePlaneFramecls
    return _ReferencePlaneFramecls


class ReferencePlaneFrame(coord.BaseCoordinateFrame):
    """A coordinate frame aligned with the reference plane coordinates of a
    Kepler orbit, centered on the barycenter or reference point of the orbit.

    See :ref:`celestial-reference-plane` for more information.

    ``ReferencePlaneFrame`` objects always have generic component names for
    spherical coordinates of ``lon``/``lat``, *not* the component names for the
    frame of the ``origin``.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other
        keywords).
    origin : `SkyCoord` or low-level coordinate object.
        The coordinate which specifies the origin of this frame. This is
        typically a sky position and a distance to the barycenter or reference
        point of the orbit at a particular epoch.

    Notes
    -----
    ``ReferencePlaneFrame`` is a factory class.  That is, the objects that it
    yields are *not* actually objects of class ``ReferencePlaneFrame``.
    Instead, distinct classes are created on-the-fly for whatever the frame
    class is of ``origin``.
    """

    origin = CoordinateAttribute(default=None, frame=None)

    def __new__(cls, *args, **kwargs):
        # We don't want to call this method if we've already set up
        # an skyoffset frame for this class.
        if not (issubclass(cls, ReferencePlaneFrame) and
                cls is not ReferencePlaneFrame):
            # We get the origin argument, and handle it here. Default is ICRS:
            # the user might want to use arbitrary reference plane coordinates
            # without every transforming them
            origin_frame = kwargs.get('origin', ICRS())

            if hasattr(origin_frame, 'frame'):
                origin_frame = origin_frame.frame
            newcls = _make_cls(origin_frame.__class__)
            return newcls.__new__(newcls, *args, **kwargs)

        # http://stackoverflow.com/questions/19277399/why-does-object-new-work-differently-in-these-three-cases
        # See above for why this is necessary. Basically, because some child
        # may override __new__, we must override it here to never pass
        # arguments to the object.__new__ method.
        if super().__new__ is object.__new__:
            return super().__new__(cls)
        return super().__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.origin is not None and not self.origin.has_data:
            raise ValueError('The origin supplied to ReferencePlaneFrame has '
                             'no data.')
        if self.has_data and hasattr(self.data, 'lon'):
            self.data.lon.wrap_angle = 180*u.deg
