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
                properties={},
                id=1
            ),
        ]

        self._generate_tile(features)

    def test_linestring(self):
        from shapely.geometry import LineString

        features = [
            dict(
                geometry=LineString([(0, 0), (1, 1)]),
                properties={},
                id=None
            ),
        ]

        self._generate_tile(features)

    def test_polygon(self):
        from shapely.geometry import Point

        features = [
            dict(
                geometry=Point(0, 0).buffer(1),
                properties={},
                id=3
            ),
        ]

        self._generate_tile(features)

    def test_multipoint(self):
        from shapely.geometry import MultiPoint

        features = [
            dict(
                geometry=MultiPoint([(0, 0), (1, 1)]),
                properties={},
                id=None
            ),
        ]

        self._generate_tile(features)

    def test_multilinestring(self):
        from shapely.geometry import MultiLineString

        features = [
            dict(
                geometry=MultiLineString([[(0, 0), (1, 0)], [(0, 1), (1, 1)]]),
                properties={},
                id=None
            ),
        ]

        self._generate_tile(features)

    def test_multipolygon(self):
        from shapely.geometry import Point

        features = [
            dict(
                geometry=Point(0, 0).buffer(0.4).union(
                    Point(1, 1).buffer(0.4)),
                properties={},
                id=4
            ),
        ]

        self._generate_tile(features)


class PropertyTest(TestCase):

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

    def test_property_types(self):
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

    def test_unicode_property_value(self):
        from shapely.geometry import Point

        features = [
            dict(
                geometry=Point(0, 0),
                properties={
                    'string': unicode('unicode_value'),
                },
                id=1
            ),
        ]

        self._generate_tile(features)

    def test_unicode_property_key(self):
        from shapely.geometry import Point

        features = [
            dict(
                geometry=Point(0, 0),
                properties={
                    unicode('unicode'): 'string_value',
                },
                id=1
            ),
        ]

        self._generate_tile(features)

    def test_unicode_layer_name(self):
        import coanacatl
        from shapely.geometry import Point

        layers = [dict(
            name=unicode('layer'),
            features=[
                dict(
                    geometry=Point(0, 0),
                    properties={
                        'foo': 'bar',
                    },
                    id=1
                ),
            ],
        )]

        bounds = (0, 0, 1, 1)
        extents = 4096

        tile_data = coanacatl.encode(layers, bounds, extents)
        self.assertTrue(tile_data)
        return tile_data


class DiscardEmptyTest(TestCase):

    def test_discard_empty_polygon(self):
        """
        Test that a polygon which becomes invalid (no rings returned from
        Wagyu) is discarded, rather than being output as an empty geometry.
        """

        import coanacatl
        from shapely.wkb import loads as wkb_loads

        shape = wkb_loads(
            '0103000000010000000400000000000000DAA90C4100000050AECA5741001F3A'
            '68F4A90C41647D6020AFCA5741001F3A68F4A90C4195DC9C73AECA5741000000'
            '00DAA90C4100000050AECA5741'.decode('hex'))

        layers = [dict(
            name=unicode('layer'),
            features=[
                dict(
                    geometry=shape,
                    properties={
                        'foo': 'bar',
                    },
                    id=1
                ),
            ],
        )]

        bounds = (156543.03392808512, 6183449.840157587,
                  234814.5508921072, 6261721.357121607)
        extents = 8192

        tile_data = coanacatl.encode(layers, bounds, extents)
        self.assertEqual('', tile_data)

    def test_discard_empty_line(self):
        import coanacatl
        from shapely.geometry import LineString

        shape = LineString([(0, 0), (1.0e-10, 0)])

        layers = [dict(
            name=unicode('layer'),
            features=[
                dict(
                    geometry=shape,
                    properties={
                        'foo': 'bar',
                    },
                    id=1
                ),
            ],
        )]

        bounds = (0, 0, 1, 1)
        extents = 4096

        tile_data = coanacatl.encode(layers, bounds, extents)
        self.assertEqual('', tile_data)
