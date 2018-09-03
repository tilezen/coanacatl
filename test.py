from unittest import TestCase


class GeometryTest(TestCase):

    def _generate_tile(self, features):
        import coanacatl

        layers = [dict(
            name='layer',
            features=features,
        )]

        bounds = (0, 0, 1, 1)
        extents = 4096

        tile_data = coanacatl.encode(layers, bounds, extents)
        self.assertTrue(tile_data)
        return tile_data

    def test_point(self):
        from shapely.geometry import Point

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
        ]

        self._generate_tile(features)

    def test_linestring(self):
        from shapely.geometry import LineString

        features = [
            dict(
                geometry=LineString([(0, 0), (1, 1)]),
                properties={'baz': 'bat'},
                id=None
            ),
        ]

        self._generate_tile(features)

    def test_polygon(self):
        from shapely.geometry import Point

        features = [
            dict(
                geometry=Point(0, 0).buffer(1),
                properties={'blah': 'blah', 'id': 123},
                id=3
            ),
        ]

        self._generate_tile(features)

    def test_multipoint(self):
        from shapely.geometry import MultiPoint

        features = [
            dict(
                geometry=MultiPoint([(0, 0), (1, 1)]),
                properties={'foo': 'bar', 'boolean': False},
                id=None
            ),
        ]

        self._generate_tile(features)

    def test_multilinestring(self):
        from shapely.geometry import MultiLineString

        features = [
            dict(
                geometry=MultiLineString([[(0, 0), (1, 0)], [(0, 1), (1, 1)]]),
                properties={'foo': 'bar'},
                id=None
            ),
        ]

        self._generate_tile(features)

    def test_multipolygon(self):
        from shapely.geometry import Point

        features = [
            dict(
                geometry=Point(0, 0).buffer(0.4).union(Point(1, 1).buffer(0.4)),
                properties={'blah': 'blah'},
                id=4
            ),
        ]

        self._generate_tile(features)
