# Coanacatl

Coanacatl is a Python extension library for writing [Mapbox Vector Tiles](https://www.mapbox.com/vector-tiles/specification/), wrapping Mapbox's [Wagyu](https://github.com/mapbox/wagyu) and [vtzero](https://github.com/mapbox/vtzero) libraries.

## How to use it?

**NOTE: This is pre-alpha proof-of-concept level software. I don't recommend you use this in production code yet!**

```python
tile_data = coanacatl.encode(layers, bounds, extents)
```

Where:

* `layers` is either a dictionary containing `features` and `name` keys, with the `features` being a list of dictionaries containing `geometry`, `properties` and `id` entries. The `geometry` must be a [Shapely](https://shapely.readthedocs.io/en/stable/) geometry object, `properties` must be a dictionary, and `id` must be either `None` or a positive integer.
* `bounds` is a 4-tuple containing the `(min_x, min_y, max_x, max_y)` bounding box of the tile.
* `extents` is a positive integer giving the number of coordinates in the `x` and `y` directions of the tile, typically 4096.

## Building and installing

You will need to install a C++11 build system and the GEOS library, e.g: if you are on Ubuntu or Debian:

```
sudo apt install build-essential libgeos-dev
```

**NOTE: probably other stuff as well! Please [file an issue](https://github.com/tilezen/coanacatl/issues/new) if you find you need additional dependencies.**

Please make sure that the git submodules are initialise and updated:

```
git submodule update --init --recursive
```

Then run the `setup.py` installation procedure:

```
python setup.py install
```

## Current limitations

* Only point, linestring, polygon and multi-versions of those are supported. Linear rings and geometry collections are currently not supported.
* Property dictionary keys must be strings, as per the MVT spec. Property dictionary values can be boolean, integer, floating point or strings.
* There are **no tests**!
* Error checking of return values from the GEOS API is inadequate, and needs shoring up.
* There needs to be a better way to return warnings/errors to the user, perhaps as a list of objects, so that the user can determine if it's enough to fail the tile or just log.

## Contributing

Coanacatl welcomes contributions! If you find a bug, please look through the [existing issues](https://github.com/tilezen/coanacatl/issues) and, if you don't find something related, [please file a new one](https://github.com/tilezen/coanacatl/issues/new). It helps if you can give as much detail as possible when writing an issue, particularly about the versions of software you are using. If you can reduce it to a reproducible test case, that's awesome and helps us a lot when we're trying to fix it!

To submit changes, please fork this project into your own organisation, make your changes on an appropriately-named branch, and submit a Pull Request. Thank you!

## What's with the name?

"Coanacatl" is the Nahuatl word for [snake meat](https://es.wiktionary.org/wiki/coanacatl), which seemed fitting as this is a Python wrapper around [Wagyu](https://github.com/mapbox/wagyu). Plus some other trimmings, but I couldn't think of a way to work [vtzero](https://github.com/mapbox/vtzero) in there as well!
