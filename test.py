import coanacatl
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry import Polygon
from shapely.geometry import MultiPoint
from shapely.geometry import MultiLineString
from shapely.geometry import MultiPolygon


features = [
    dict(
        geometry=Point(0, 0),
        properties={
            'string': 'string_value',
            'long': 4294967297L,
            'int': 1,
            'float': 1.0,
            'bool': True,
        },
        id=1
    ),
    dict(
        geometry=LineString([(0, 0), (1, 1)]),
        properties={'baz': 'bat'},
        id=None
    ),
    dict(
        geometry=Point(0, 0).buffer(1),
        properties={'blah': 'blah', 'id': 123},
        id=3
    ),
    dict(
        geometry=MultiPoint([(0, 0), (1, 1)]),
        properties={'foo': 'bar', 'boolean': False},
        id=None
    ),
    dict(
        geometry=MultiLineString([[(0, 0), (1, 0)], [(0, 1), (1, 1)]]),
        properties={'foo': 'bar'},
        id=None
    ),
    dict(
        geometry=Point(0, 0).buffer(0.4).union(Point(1, 1).buffer(0.4)),
        properties={'blah': 'blah'},
        id=4
    ),
]

layers = [dict(
    name='layer',
    features=features,
)]

bounds = (0, 0, 1, 1)
extents = 4096

tile_data = coanacatl.encode(layers, bounds, extents)
print repr(tile_data)
with open('foo.mvt', 'w') as fh:
    fh.write(tile_data)
