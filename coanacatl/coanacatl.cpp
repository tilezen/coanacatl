#include <boost/python.hpp>

#include <string>
#include <memory>
#include <set>
#include <iostream>
#include <fstream>

#include <geos_c.h>

#include <mapbox/geometry/wagyu/wagyu.hpp>
#include <vtzero/builder.hpp>

namespace mg = mapbox::geometry;
namespace mgw = mapbox::geometry::wagyu;
namespace bp = boost::python;

typedef int coord;
typedef mg::point<coord> point;
typedef mg::linear_ring<coord> linear_ring;
typedef mg::polygon<coord> polygon;
typedef mg::multi_polygon<coord> multi_polygon;
typedef mgw::wagyu<coord> wagyu;

// set aside this special value of fid to represent None, as we often will set
// IDs to None and keep the real ID in the properties.
const uint64_t FID_NONE = std::numeric_limits<uint64_t>::max();

#if GEOS_CAPI_VERSION_MINOR >= 9
#define INIT_GEOS GEOS_init_r()
#define FINISH_GEOS GEOS_finish_r
#else
#include <stdarg.h>
void _coanacatl_printf(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  vprintf(fmt, ap);
  va_end(ap);
}
#define INIT_GEOS initGEOS_r(_coanacatl_printf, _coanacatl_printf)
#define FINISH_GEOS finishGEOS_r
#endif

namespace {

/**
 * Extract a string from a Python object. The object _must_ be either a str or
 * unicode object, else an exception will be thrown.
 */
std::string extract_utf8_string(bp::object value) {
  PyObject *value_ptr = value.ptr();

  if (PyUnicode_Check(value_ptr)) {
    bp::object encoded = bp::str(value).encode("utf-8");
    std::string v = bp::extract<std::string>(encoded);
    return v;

  } else if (PyString_Check(value_ptr)) {
    std::string v = bp::extract<std::string>(value);
    return v;

  } else {
    std::ostringstream out;
    bp::object repr_py = value.attr("__repr__")();
    std::string repr = bp::extract<std::string>(repr_py);
    out << "Unable to convert Python object of type "
        << value_ptr->ob_type->tp_name << " to string: " << repr;
    throw std::runtime_error(out.str());
  }
}

} // end anonymous namespace

class encoder {
public:
  encoder(bp::tuple bounds, size_t extents)
    : m_extents(extents)
    , m_geos_ctx(INIT_GEOS) {
    if (bp::len(bounds) != 4) {
      throw std::runtime_error("Bounds tuple must have 4 elements.");
    }
    m_minx = bp::extract<double>(bounds[0]);
    m_miny = bp::extract<double>(bounds[1]);
    m_maxx = bp::extract<double>(bounds[2]);
    m_maxy = bp::extract<double>(bounds[3]);
  }

  ~encoder() {
    FINISH_GEOS(m_geos_ctx);
  }

  void encode_layer(bp::object layer);

  void serialize(std::string &buffer) const {
    m_tile_builder.serialize(buffer);
  }

private:
  size_t m_extents;
  GEOSContextHandle_t m_geos_ctx;
  double m_minx, m_miny, m_maxx, m_maxy;
  vtzero::tile_builder m_tile_builder;
  std::set<std::string> m_layer_names;

  void add_id(vtzero::feature_builder &fb, uint64_t fid);
  void add_properties(vtzero::feature_builder &fb, bp::dict props);
  vtzero::point translate(double x, double y);
  void encode_feature(
    const GEOSGeom_t *geometry,
    bp::dict props,
    uint64_t fid,
    vtzero::layer_builder &lb);
  void encode_point(
    const GEOSGeom_t *geometry,
    bp::dict props,
    uint64_t fid,
    vtzero::layer_builder &lb);
  void encode_linestring(
    const GEOSGeom_t *geometry,
    bp::dict props,
    uint64_t fid,
    vtzero::layer_builder &lb);
  void encode_polygon(
    const GEOSGeom_t *geometry,
    bp::dict props,
    uint64_t fid,
    vtzero::layer_builder &lb);
  void encode_multi_point(
    const GEOSGeom_t *geometry,
    bp::dict props,
    uint64_t fid,
    vtzero::layer_builder &lb);
  void encode_multi_linestring(
    const GEOSGeom_t *geometry,
    bp::dict props,
    uint64_t fid,
    vtzero::layer_builder &lb);
  void encode_multi_polygon(
    const GEOSGeom_t *geometry,
    bp::dict props,
    uint64_t fid,
    vtzero::layer_builder &lb);
  polygon geom_to_poly(const GEOSGeom_t *geometry);
  linear_ring ring_to_mbg(const GEOSGeom_t *linear_ring);
  std::vector<vtzero::point> distinct_linestring(const GEOSGeom_t *geometry);
  void encode_wagyu_result(
    wagyu &w,
    bp::dict props,
    uint64_t fid,
    vtzero::layer_builder &lb);
};

void encoder::encode_layer(bp::object layer) {
  std::string layer_name = extract_utf8_string(layer["name"]);

  if (m_layer_names.count(layer_name) > 0) {
    throw std::runtime_error("Duplicate layer names are not allowed.");
  }
  m_layer_names.insert(layer_name);

  vtzero::layer_builder lb(m_tile_builder, layer_name, 2, m_extents);

  bp::object features = layer["features"];
  const size_t num_features = bp::len(features);
  for (size_t i = 0; i < num_features; ++i) {
    bp::object feature = features[i];
    bp::object shape = feature["geometry"];
    bp::dict props = bp::extract<bp::dict>(feature["properties"]);
    bp::object py_fid = feature["id"];
    uint64_t fid = py_fid.is_none() ? FID_NONE : bp::extract<uint64_t>(py_fid);

    // extract the shapely geometry pointer. this is probably very fragile.
    size_t geom_ptr = bp::extract<size_t>(shape.attr("_geom"));
    const GEOSGeom_t *geometry = reinterpret_cast<const GEOSGeom_t*>(geom_ptr);

    encode_feature(geometry, props, fid, lb);
  }
}

void encoder::add_id(vtzero::feature_builder &fb, uint64_t fid) {
  if (fid != FID_NONE) {
    fb.set_id(fid);
  }
}

void encoder::add_properties(vtzero::feature_builder &fb, bp::dict props) {
  bp::list items = props.items();
  const size_t num_items = bp::len(items);
  for (size_t i = 0; i < num_items; ++i) {
    bp::object item = items[i];
    std::string k = extract_utf8_string(item[0]);
    bp::object value = item[1];
    PyObject *value_ptr = value.ptr();

    if (PyBool_Check(value_ptr)) {
      bool v = bp::extract<bool>(value);
      fb.add_property(k, v);

    } else if (PyFloat_Check(value_ptr)) {
      double v = bp::extract<double>(value);
      fb.add_property(k, v);

    } else if (PyLong_Check(value_ptr) || PyInt_Check(value_ptr)) {
      int64_t v = bp::extract<int64_t>(value);
      fb.add_property(k, v);

    } else if (PyUnicode_Check(value_ptr) || PyString_Check(value_ptr)) {
      std::string v = extract_utf8_string(value);
      fb.add_property(k, v);

    } else {
      // we don't understand how to write out this type, so we just ignore
      // it for now.

      // std::ostringstream out;
      // bp::object repr_py = value.attr("__repr__")();
      // std::string repr = bp::extract<std::string>(repr_py);
      // out << "Unable to handle Python object of type "
      //     << value_ptr->ob_type->tp_name << ": " << repr;
      // throw std::runtime_error(out.str());
    }
  }
}

vtzero::point encoder::translate(double x, double y) {
  int quant_x = double(m_extents) * (x - m_minx) / (m_maxx - m_minx);
  // remember: MVT tile coordinates are y-down, "screen" space. mercator
  // coordinates are y-up, "world" space, so we need to flip the y coord.
  int quant_y = double(m_extents) * (m_maxy - y) / (m_maxy - m_miny);
  return vtzero::point(quant_x, quant_y);
}

void encoder::encode_point(
  const GEOSGeom_t *geometry,
  bp::dict props,
  uint64_t fid,
  vtzero::layer_builder &lb) {

  vtzero::point_feature_builder fb{lb};
  add_id(fb, fid);

  double x, y;
  GEOSGeomGetX_r(m_geos_ctx, geometry, &x);
  GEOSGeomGetY_r(m_geos_ctx, geometry, &y);
  vtzero::point p = translate(x, y);
  fb.add_point(p);

  add_properties(fb, props);
}

std::vector<vtzero::point> encoder::distinct_linestring(
  const GEOSGeom_t *geometry) {
  // NOTE: coord_seq is owned by the geometry, so we do not free it.
  const GEOSCoordSequence *coord_seq = GEOSGeom_getCoordSeq_r(m_geos_ctx, geometry);
  if (coord_seq == nullptr) {
    throw std::runtime_error("Error calling GEOSGeom_getCoordSeq_r");
  }
  unsigned int num_points = 0;
  int status = GEOSCoordSeq_getSize_r(m_geos_ctx, coord_seq, &num_points);
  if (status == 0) {
    throw std::runtime_error("Error calling GEOSCoordSeq_getSize_r");
  }
  if (num_points < 2) {
    throw std::runtime_error("Must have at least 2 points in a linestring");
  }

  // only use _distinct_ points. it's not allowed to have two adjacent
  // points with the same coordinates.
  std::vector<vtzero::point> distinct;
  distinct.reserve(num_points);

  double x, y;
  for (unsigned int i = 0; i < num_points; ++i) {
    GEOSCoordSeq_getX_r(m_geos_ctx, coord_seq, i, &x);
    GEOSCoordSeq_getY_r(m_geos_ctx, coord_seq, i, &y);
    vtzero::point p = translate(x, y);
    if (i == 0 || p != distinct[i-1]) {
      distinct.push_back(p);
    }
  }

  return distinct;
}

void encoder::encode_linestring(
  const GEOSGeom_t *geometry,
  bp::dict props,
  uint64_t fid,
  vtzero::layer_builder &lb) {

  auto distinct = distinct_linestring(geometry);

  // ignore a geometry with not enough distinct points.
  // TODO: we should log this!
  if (distinct.size() >= 2) {
    vtzero::linestring_feature_builder fb{lb};
    add_id(fb, fid);

    fb.add_linestring_from_container(distinct);

    add_properties(fb, props);
  }
}

linear_ring encoder::ring_to_mbg(const GEOSGeom_t *geos_ring) {
  // NOTE: coord_seq is owned by the geometry, so we do not free it.
  const GEOSCoordSequence *coord_seq = GEOSGeom_getCoordSeq_r(m_geos_ctx, geos_ring);
  if (coord_seq == nullptr) {
    throw std::runtime_error("Error calling GEOSGeom_getCoordSeq_r");
  }
  unsigned int num_points = 0;
  int status = GEOSCoordSeq_getSize_r(m_geos_ctx, coord_seq, &num_points);
  if (status == 0) {
    throw std::runtime_error("Error calling GEOSCoordSeq_getSize_r");
  }
  if (num_points < 2) {
    throw std::runtime_error("Must have at least 2 points in a linear ring");
  }

  double x, y;
  linear_ring ring;
  ring.reserve(num_points);
  for (unsigned int i = 0; i < num_points; ++i) {
    GEOSCoordSeq_getX_r(m_geos_ctx, coord_seq, i, &x);
    GEOSCoordSeq_getY_r(m_geos_ctx, coord_seq, i, &y);
    vtzero::point p = translate(x, y);
    ring.push_back(point(p.x, p.y));
  }

  return ring;
}

polygon encoder::geom_to_poly(const GEOSGeom_t *geometry) {
  polygon p;

  p.emplace_back(ring_to_mbg(GEOSGetExteriorRing_r(m_geos_ctx, geometry)));

  const int num_int_rings = GEOSGetNumInteriorRings_r(m_geos_ctx, geometry);
  if (num_int_rings < 0) {
    throw std::runtime_error("Error calling GEOSGetNumInteriorRings_r");
  }

  for (int i = 0; i < num_int_rings; ++i) {
    p.emplace_back(ring_to_mbg(GEOSGetInteriorRingN_r(m_geos_ctx, geometry, i)));
  }

  return p;
}

void encoder::encode_polygon(
  const GEOSGeom_t *geometry,
  bp::dict props,
  uint64_t fid,
  vtzero::layer_builder &lb) {

  wagyu w;
  polygon poly = geom_to_poly(geometry);
  w.add_polygon(poly);

  encode_wagyu_result(w, props, fid, lb);
}

void encoder::encode_multi_point(
  const GEOSGeom_t *multi_geometry,
  bp::dict props,
  uint64_t fid,
  vtzero::layer_builder &lb) {

  vtzero::point_feature_builder fb{lb};
  add_id(fb, fid);

  const int num_geoms = GEOSGetNumGeometries_r(m_geos_ctx, multi_geometry);
  if (num_geoms < 0) {
    throw std::runtime_error("Error calling GEOSGetNumGeometries_r");
  }

  fb.add_points(num_geoms);
  for (int i = 0; i < num_geoms; ++i) {
    const GEOSGeom_t *geom = GEOSGetGeometryN_r(m_geos_ctx, multi_geometry, i);
    if (geom == nullptr) {
      throw std::runtime_error("Error calling GEOSGetGeometryN_r");
    }
    double x, y;
    GEOSGeomGetX_r(m_geos_ctx, geom, &x);
    GEOSGeomGetY_r(m_geos_ctx, geom, &y);
    vtzero::point p = translate(x, y);
    fb.set_point(p);
  }

  add_properties(fb, props);
}

void encoder::encode_multi_linestring(
  const GEOSGeom_t *multi_geometry,
  bp::dict props,
  uint64_t fid,
  vtzero::layer_builder &lb) {

  bool output_started = false;
  vtzero::linestring_feature_builder fb{lb};

  const int num_geoms = GEOSGetNumGeometries_r(m_geos_ctx, multi_geometry);
  if (num_geoms < 0) {
    throw std::runtime_error("Error calling GEOSGetNumGeometries_r");
  }

  for (int i = 0; i < num_geoms; ++i) {
    const GEOSGeom_t *geom = GEOSGetGeometryN_r(m_geos_ctx, multi_geometry, i);
    if (geom == nullptr) {
      throw std::runtime_error("Error calling GEOSGetGeometryN_r");
    }
    auto distinct = distinct_linestring(geom);
    // ignore degenerate linestrings without at least 2 distinct points.
    // TODO: log warnings.
    if (distinct.size() >= 2) {
      // lazy output of ID, in case we don't actually end up outputting anything
      // at all, then we can just ignore the whole feature.
      if (!output_started) {
        add_id(fb, fid);
        output_started = true;
      }

      fb.add_linestring_from_container(distinct);
    }
  }

  if (output_started) {
    add_properties(fb, props);
  }
}

void encoder::encode_multi_polygon(
  const GEOSGeom_t *multi_geometry,
  bp::dict props,
  uint64_t fid,
  vtzero::layer_builder &lb) {

  wagyu w;

  const int num_geoms = GEOSGetNumGeometries_r(m_geos_ctx, multi_geometry);
  if (num_geoms < 0) {
    throw std::runtime_error("Error calling GEOSGetNumGeometries_r");
  }

  for (int i = 0; i < num_geoms; ++i) {
    const GEOSGeom_t *geom = GEOSGetGeometryN_r(m_geos_ctx, multi_geometry, i);
    if (geom == nullptr) {
      throw std::runtime_error("Error calling GEOSGetGeometryN_r");
    }
    w.add_polygon(geom_to_poly(geom));
  }

  encode_wagyu_result(w, props, fid, lb);
}

void encoder::encode_wagyu_result(
  wagyu &w,
  bp::dict props,
  uint64_t fid,
  vtzero::layer_builder &lb) {

  multi_polygon result;
  bool ok = w.execute(mgw::clip_type_union, result,
    mgw::fill_type_even_odd,
    mgw::fill_type_even_odd);

  if (!ok || result.empty()) {
    // this is probably because the polygon ended up being degenerate in integer
    // coordinates.

    // TODO: warning?
    return;
  }

  vtzero::polygon_feature_builder fb{lb};
  add_id(fb, fid);

  for (const polygon &p : result) {
    for (const linear_ring &lr : p) {
      fb.add_ring_from_container(lr);
    }
  }

  add_properties(fb, props);
}

void encoder::encode_feature(
  const GEOSGeom_t *geometry,
  bp::dict props,
  uint64_t fid,
  vtzero::layer_builder &lb) {

  int geom_type = GEOSGeomTypeId_r(m_geos_ctx, geometry);

  switch (geom_type) {
  case GEOS_POINT:
    encode_point(geometry, props, fid, lb);
    break;

  case GEOS_LINESTRING:
    encode_linestring(geometry, props, fid, lb);
    break;

  case GEOS_POLYGON:
    encode_polygon(geometry, props, fid, lb);
    break;

  case GEOS_MULTIPOINT:
    encode_multi_point(geometry, props, fid, lb);
    break;

  case GEOS_MULTILINESTRING:
    encode_multi_linestring(geometry, props, fid, lb);
    break;

  case GEOS_MULTIPOLYGON:
    encode_multi_polygon(geometry, props, fid, lb);
    break;

  case GEOS_LINEARRING:
  case GEOS_GEOMETRYCOLLECTION:
  default:
    std::ostringstream out;
    out << "Unable to handle geometry type " << geom_type << ".";
    throw std::runtime_error(out.str());
  }
}

std::string encode(bp::list layers, bp::tuple bounds, size_t extents) {
  encoder enc(bounds, extents);

  const size_t len = bp::len(layers);
  for (size_t i = 0; i < len; ++i) {
    bp::object layer = layers[i];
    enc.encode_layer(layer);
  }

  std::string buffer;
  enc.serialize(buffer);
  return buffer;
}

BOOST_PYTHON_MODULE(coanacatl) {
  using namespace boost::python;

  def("encode", encode);
}
