#include "probe.h"
#include "polygon.h"
#include <vector>
#include <math.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
using namespace std;


#define PI 3.1415926535
#define SWEEPANGLE (26.35 * (PI / 180))
#define MINRANGE         (6.5 / 60. * imgscale)
#define MAXRANGE         (10. / 60. * imgscale)
#define MINRANGE_CADC14         (4.5 / 60. * imgscale)
#define MAXRANGE_CADC14         (7. / 60. * imgscale)

double imgscale;  // mm per degree (we ignore distortion)

enum Mode { ModeGCLEF, ModeM3, ModeDGNF, ModeDGWF, ModeGMACS };

extern Mode mode;


vector<string> split(const string &text, char sep) {
    vector<string> tokens;
    size_t start = 0, end = 0;
    while ((end = text.find(sep, start)) != string::npos) {
        tokens.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    tokens.push_back(text.substr(start));
    return tokens;
}

Point scale(Point v, double m) {
  Point scaled(v.x * m, v.y * m);
  return scaled;
}

Probe::Probe() { }

/**
   Read a polygon from a textfile. Format of the textfile is a single x,y pair per
   line, separated by a single tab.
**/
Polygon load_poly(string filename) {
  ifstream infile(filename);
  vector<string> points;
  string line;
  Polygon poly;

  for ( ; getline(infile, line); ) {
    vector<string> points = split(line, '\t');
    Point current_pt(stod(points[0]), stod(points[1]));

    poly.add_pt(current_pt);
  }

  return poly;
}

// Point calculate_shaft_front(Polygon slider_shaft) {
//   return Point(0, slider_shaft.min_abs_y());
// }
// 
// Point calculate_baffle_ctr(Polygon baffle_tube) {
//   return Point(0, baffle_tube.min_abs_y());
// }

/**
   Probe constructor.

   - First loads the polygons that make up the probe. These polygons are rotated by the
     given angle.
   - The 'Base_' polygons are set so that the probe position can be reset after being moved.
   - Radius us hardcoded in the constructor.
   - rotate_about is the point the probe is translated to before being rotated when it is being
     moved over a star. For the configuration where each probe rotates about its own pivot,
     this attribute should be the origin.
**/
Probe::Probe(double _angle, string slider_body_file,
             string slider_shaft_file, string baffle_tube_file) {
  angle   = _angle;

  Base_SliderBody  = load_poly(slider_body_file).rotate(angle * (PI / 180));
  Base_SliderShaft = load_poly(slider_shaft_file).rotate(angle * (PI / 180));
  Base_BaffleTube  = load_poly(baffle_tube_file).rotate(angle * (PI / 180));

  SliderBody  = Base_SliderBody;
  SliderShaft = Base_SliderShaft;
  BaffleTube  = Base_BaffleTube;

  parts.push_back(SliderBody);
  parts.push_back(BaffleTube);  
  parts.push_back(SliderShaft);

  // SliderShaftFront = calculate_shaft_front(SliderShaft).rotate(angle * (PI / 180));
  // BaffleTubeCtr    = calculate_baffle_ctr(BaffleTube).rotate(angle * (PI / 180));

  SliderShaftFront = Point(0, 69.062).rotate(angle * (PI / 180));
  BaffleTubeCtr = Point(0, 409).rotate(angle * (PI / 180));

  Point origin_vector(0, 1000);

  // The radial position of the rotary axis.
  radius = 1358;

  // BAM: Not sure what the 1.3 is for.  I'm changing radius from 1300 to 1400, so I'll change this also...
  //  axis = scale(origin_vector.rotate(angle * (PI / 180)), 1.300);
  axis = scale(origin_vector.rotate(angle * (PI / 180)), radius / 1000.0);
  center = axis;
  default_star = Star(axis.rotate(angle).x, axis.rotate(angle).y, 20, 20, 0);

  needs_transfer = false;

  Point origin(0, 0);
  rotate_about = origin;
}

Probe::~Probe(void) { }

double vector_length(Point pt) {
  return sqrt( pow(pt.x, 2) + pow(pt.y, 2) );
}

double abs_angle_between_vectors(Point u, Point v) {
  return acos( ((u.x * v.x) + (u.y * v.y)) /
               (vector_length(v) * vector_length(u)) );
}

/**
   The absolute value of the angle between vector can be found by the formula
   given in abs_angle_between_vectors, but in order for the probe to move in the
   correct direction, we need to check how the probe's axis and the star are positioned
   relative to each other.
**/
double angle_between_vectors(Point u, Point v) {
  Point origin(0, 1);
  int mod = -1;

  double u_angle, v_angle;
  u_angle = abs_angle_between_vectors(origin, u);
  v_angle = abs_angle_between_vectors(origin, v);

  int u_quadrant = quadrant(u);
  int v_quadrant = quadrant(v);

  if (u_quadrant == 1) {
    if (v_quadrant == 1 && u_angle < v_angle) {
      mod = 1;
    } else if (v_quadrant == 4) {
      mod = 1;
    }
  } else if (u_quadrant == 2) {
    if (v_quadrant == 2 && v_angle < u_angle) {
      mod = 1;
    }
  } else if (u_quadrant == 3) {
    if (v_quadrant == 3 && v_angle < u_angle) {
      mod = 1;
    }
  } else if (u_quadrant == 4) {
    if (v_quadrant == 4 && u_angle < v_angle) {
      mod = 1;
    } else if (v_quadrant == 3) {
      mod = 1;
    }
  }
  
  return abs_angle_between_vectors(u, v) * mod;
}

Polygon translate_poly(Polygon polygon, Point from, Point to) {
  Polygon translated_poly;
  vector<Point>::iterator it;

  for(it = polygon.points.begin(); it != polygon.points.end(); it++) {
    translated_poly.add_pt(it->translate(from, to));
  }

  return translated_poly;
}

Polygon rotate_poly(Polygon poly, double theta) {
  Polygon rotated_poly;
  vector<Point>::iterator it;
  
  for(it = poly.points.begin(); it != poly.points.end(); ++it) {
    rotated_poly.add_pt(it->rotate(theta));
  }

  return rotated_poly;
}

double distance(Point a, Point b) {
  return sqrt(pow(b.x - a.x, 2) + pow(b.y - a.y, 2));
}

int quadrant(Point pt) {
  if (pt.x >= 0 && pt.y >= 0) { return 1; }
  if (pt.x <= 0 && pt.y >= 0) { return 2; }
  if (pt.x <= 0 && pt.y <= 0) { return 3; }
  return 4;
}

// Angle between probe's axis and the point represented as a vector
// from the origin.
double Probe::angle_to_point(Point pt) {
  Point origin(0, 1);
  int mod = -1;

  double axis_angle, pt_angle;
  axis_angle = abs_angle_between_vectors(origin, axis);
  pt_angle   = abs_angle_between_vectors(origin, pt);

  int pt_quadrant   = quadrant(pt);
  int axis_quadrant = quadrant(axis);

  if (axis_quadrant == 1) {
    if (pt_quadrant == 1 && axis_angle < pt_angle) {
      mod = 1;
    } else if (pt_quadrant == 4) {
      mod = 1;
    }
  } else if (axis_quadrant == 2) {
    if (pt_quadrant == 2 && pt_angle < axis_angle) {
      mod = 1;
    }
  } else if (axis_quadrant == 3) {
    if (pt_quadrant == 3 && pt_angle < axis_angle) {
      mod = 1;
    }
  } else if (axis_quadrant == 4) {
    if (pt_quadrant == 4 && axis_angle < pt_angle) {
      mod = 1;
    } else if (pt_quadrant == 3) {
      mod = 1;
    }
  }

  Point u(pt.x - center.x, pt.y - center.y);

  return abs_angle_between_vectors(u, scale(center, -1)) * mod;
}

/**
mm -> degrees -> arcminutes
(mm / 3600) * 60
**/

/**
   The distance between the star and the front of the probe's slider shaft is 
   given in SAO-AGWS-DOC-0021-AGWSFocalPlaneParameters rev 3
**/
double baffle_separation(Point p) {
  Point origin (0,0);
  Point p_location_arcminutes((p.x / 3600) * 60, (p.y / 3600) * 60);
  double r = distance(p_location_arcminutes, origin);

  if (mode == ModeDGWF) {
      return 189.109 + 1.06287366 * r * r + 1.5212314e-4 * r * r * r * r;
  } else if (mode == ModeGMACS) {
    return 235.289 + 1.35076769 * r *r + 2.53369431e-3 * r * r * r * r;
  } else {
      return  69.062 + 0.83398047 * r * r + 3.112125e-5  * r * r * r * r;
  }
}

vector<Point> get_points(vector<Polygon> polygons) {
  vector<Point> res;
  for (Polygon poly : polygons) {
    for (Point pt : poly.points) {
      res.push_back(pt);
    }
  }

  return res;
}

/**
   Both the slider shaft and the slider body can be positioned correctly using this function.
   The positions of the slider polygons depend on the baffle separation, which is given as a parameter.
   
   - The polygons of the slider shaft and slider body are rotated the angle theta, no matter what. This
     rotation happens outside this function, but it is useful to know that it doesn't need to be calculated
     based on baffle separation.
   - The amount of translation forward or backward along the axis made by the origin (0, 0) and the star
     is determined by how far the slider must be from the star (the baffle separation).
**/
Polygon Probe::move_slider(Polygon slider, double theta, Point pivot, double baffle_sep) {
  Point w = SliderShaftFront.translate(center, rotate_about).rotate(theta).translate(rotate_about, center);
  Point v = Point(pivot.x - w.x, pivot.y - w.y);

  int mod = 1;
  if (distance(w, center) < distance(pivot, center)) {
    mod = -1;
  }

  Point u = scale(v, (v.length() + (mod * baffle_sep)) / v.length());

  Polygon newslider;
  Point origin(0, 0);

  for (int i=0; i<slider.points.size(); i++) {
    newslider.add_pt(slider.points[i].translate(origin, u));
  }

  return newslider;
}

/**
   Translate the baffle tube so that it can be rotated into position over the point pt.
   
   - Create a point, p, that is the given point translated to have the probe center as its origin.
   - Find the angle between the probe's axis and this new point. (this is the angle the probe will
     have to rotate to be positioned over pt.)
   - Find the point on the probe's axis where the p would intersect if rotated around the probe's center.
   - The distance between this intersection and the probe's center is the amount the baffle tube must be
     translated.
**/
void Probe::position_baffle_tube(Point pt) {
  Point origin(0, 0);
  Point c = scale(center, -1);
  Point p = pt.translate(origin, c);
  double theta = angle_between_vectors(p, c);
  Point intersection = p.rotate(theta);

  double distance_from_center = distance(origin, c) - distance(origin, intersection);

  Point unit_axis = axis.normalize();
  for (int i=0; i<BaffleTube.points.size(); i++) {
    BaffleTube.points[i] = BaffleTube.points[i].translate(origin, scale(unit_axis, distance_from_center));
  }

  parts[1] = BaffleTube;
}

void Probe::reset_parts() {
  SliderBody = Base_SliderBody;
  SliderShaft = Base_SliderShaft;
  BaffleTube = Base_BaffleTube;

  parts.clear();

  parts.push_back(SliderBody);
  parts.push_back(BaffleTube);
  parts.push_back(SliderShaft);
}

/**
   Carry out the translation and rotation of all the probe's parts so that it is
   positioned over the given point, pivot.
**/
vector<Polygon> Probe::transform_parts(Point pivot) {
  vector<Polygon> transformed_parts;
  Point origin(0, 0);
  
  double dist_star_from_center   = distance(pivot, origin);
  double dist_slider_from_center = dist_star_from_center + baffle_separation(pivot);

  position_baffle_tube(pivot);

  double theta = angle_to_point(pivot);

  for (Polygon part : parts) {
    Polygon translated_poly = translate_poly(part, center, rotate_about);
    Polygon rotated_poly    = rotate_poly(translated_poly, theta);
    Polygon retranslated    = translate_poly(rotated_poly, rotate_about, center);

    transformed_parts.push_back(retranslated);
  }

  transformed_parts[0] = move_slider(transformed_parts[0], theta, pivot, baffle_separation(pivot));
  transformed_parts[2] = move_slider(transformed_parts[2], theta, pivot, baffle_separation(pivot));

  reset_parts();

  return transformed_parts;
}

/**
   Code derived from explanation given here: http://paulbourke.net/geometry/circlesphere/
   Returns the two points (or none, if there are no intersections) of the circles described
   by the arguments.
**/
vector<Point> circle_intersections(Point P0, double r0, Point P1, double r1) {
  double d = distance(P0, P1);
  vector<Point> intersections;

  if (d < (r0 + r1)) {
    double a = (pow(r0, 2) - pow(r1, 2) + pow(d, 2)) / (2 * d);
    double b = d - a;

    double h = sqrt(pow(r0, 2) - pow(a, 2));

    Point P2(P0.x + a * (P1.x - P0.x) / d,
             P0.y + a * (P1.y - P0.y) / d);

    double i1x = P2.x + h * (P1.y - P0.y) / d;
    double i1y = P2.y - h * (P1.x - P0.x) / d;
    double i2x = P2.x - h * (P1.y - P0.y) / d;
    double i2y = P2.y + h * (P1.x - P0.x) / d;

    Point i1(i1x, i1y);
    Point i2(i2x, i2y);

    intersections.push_back(i1);
    intersections.push_back(i2);
  }

  return intersections;
}

/**
   Returns the amount of angle (in rads) the probe can track the given star.

   - How far a probe can track a star is equal to the angle between the star's
     current position and the position it will have when it leaves the probe's range.
   - The position a star will have when it leaves the probe's range is found by the clockwise
     point in the intersection of two circles: the one created by the probe's center and its
     radius, and the one created by the origin and the radius given by the distance between
     the star and the origin.
**/
double Probe::track_distance(Star s) {
  Point origin(0, 0);
  
  vector<Point> range_extremes = circle_intersections(center, radius, origin, distance(origin, s.point()));

  if (range_extremes.size() < 2) {
    return 0;
  }

  Point star_exit_range;
  if (angle_to_point(range_extremes[0]) < 0) {
    star_exit_range = range_extremes[0];
  } else {
    star_exit_range = range_extremes[1];
  }

  Point u(star_exit_range.x - axis.x, star_exit_range.y - axis.y);
  Point v(s.point().x - axis.x, s.point().y - axis.y);

  return abs_angle_between_vectors(v, u);
}

int partition(vector<Star>& stars, int p, int q, Point center, Star current_star)
{
  Point probe_vector = Point(current_star.x - center.x, current_star.y - center.y);
  Star x = stars[p];
  Point x_vector = Point(x.x - center.x, x.y - center.y);
  int i = p;
  int j;

  for(j=p+1; j<q; j++) {
    Point s_vector = Point(stars[j].x - center.x, stars[j].y - center.y);
    
    if (abs_angle_between_vectors(probe_vector, s_vector) >= abs_angle_between_vectors(probe_vector, x_vector)) {
      i = i+1;
      swap(stars[i], stars[j]);
    }

  }
  
  swap(stars[i], stars[p]);

  return i;
}

// Implementation of quicksort on Star objects. Sorts by magnitude.
void sort_by_transfer_distance(vector<Star>& stars, int p, int q, Point center, Star current_star)
{
    int r;
    if(p < q)
    {
      r = partition(stars, p, q, center, current_star);
      sort_by_transfer_distance(stars, p, r, center, current_star);  
      sort_by_transfer_distance(stars, r+1, q, center, current_star);
    }
}

/**
   Ask a probe to track for a given distance (in degrees).

   - If the probe can follow its base star for the given distance, it does so by adding
     the given distance to its 'distance_tracked' attribute.
   - If the probe cannot follow its base star for the given distance, it looks for a
     star that will be in its range after the starfield has rotated the given distance.
   - If such a star is found, the probe switches its base star and current star.
   - If no backtrack is found, the probe returns a -1 telling the caller it was unable to
     track anything for the given distance.
**/
int Probe::track(double dist) {
  if (! in_range(base_star.rotate(dist))) {

    sort_by_transfer_distance(backward_transfers, 0, backward_transfers.size(), center, current_star);

    bool found = false;
    for (int i=0; i<backward_transfers.size(); i++) {
      Star s = backward_transfers[i];

      if (in_range(s.rotate(dist))) {
        found = true;
        base_star = s;
        current_star  = s.rotate(dist);

        break;
      }
    }

    if (found == false) {
      return -1;
    }

  } else {
    current_star = base_star.rotate(dist);
  }

  distance_tracked = dist;

  return 0;
}

/**
   Return the two points where a star's rotation circle will intersect with the
   probe's tracking circle.
**/
vector<Point> Probe::get_range_extremes(Star s) {
  Point origin(0, 0);
  vector<Point> intersections = circle_intersections(center, radius, origin, distance(origin, s.point()));

  return intersections;
}

/**
   Is the given star in the annulus created by the 6 and 10 arc minute constraints?
**/
bool safe_distance_from_center(Star star) {
  Point star_pt(star.x, star.y), origin(0, 0);
    
    if (mode == ModeGCLEF) {
        imgscale = 3600 * 0.987; //mm per degree
    } else if (mode == ModeM3) {
        imgscale = 3600 * 0.987;
    } else if (mode == ModeDGNF) {
        imgscale = 3600 * 0.987;
    } else if (mode == ModeDGWF) {
        imgscale = 3600 * 1.04938;
    } else if (mode == ModeGMACS) {
        imgscale = 3600 * 65.02345/60.;
    } else {
        cerr << "Unknown mode:" << mode;
        exit(1);
        }
    
  double dist = distance(star_pt, origin); //in mm

    if (mode == ModeGMACS){
        cerr << dist << " CADC14 " << MINRANGE_CADC14 << " " << MAXRANGE_CADC14 << "\n";
      if (MINRANGE_CADC14 < dist && dist < MAXRANGE_CADC14) { //compare in mm
        return true;
      } else {
        return false;
      }
    } else {
        cerr << dist << " " << MINRANGE << " " << MAXRANGE << "\n";
      if (MINRANGE < dist && dist < MAXRANGE) {
        return true;
      } else {
        return false;
      }
    }
}

/**
   Determine whether or not the given star is inside the probe's sweep.
   
   - If the star's rotation circle does not intersect with the probe's sweep,
     the star is not within range.
   - If the star is further from the probe than the probe's radius, the star
     is not within range.
   - If the star is not within the constrained annulus, the star is not within
     range.
   - If the star is not within the +-30 degree probe sweep angle constraint, the
     star is not within range.
**/
bool Probe::in_range(Star s) {
  vector<Point> range_extremes = get_range_extremes(s);
  if (range_extremes.size() < 2) {
    return false;
  }

    cerr << "(" << center.x << ", " << center.y << ")  ("
        << s.point().x << ", " << s.point().y << "), dist = "
        << distance(center, s.point()) << ", radius = " << radius << "\n";
  if (distance(center, s.point()) > radius) { //radius = 1358mm
    return false;
  }

  if (! safe_distance_from_center(s)) {
    return false;
  }

  Point star_enter_range, star_exit_range;
  if (angle_between_vectors(range_extremes[0], axis) < 0) {
    star_exit_range  = range_extremes[0];
    star_enter_range = range_extremes[1];
  } else {
    star_exit_range  = range_extremes[1];
    star_enter_range = range_extremes[0];
  }

  if (abs(angle_to_point(s.point())) < SWEEPANGLE) {
    return true;
  }

  return false;
}

/**
   What amount of rotation (in radians) before a star is within the probe's range?
**/
double Probe::distance_till_inrange(Star s) {
  vector<Point> range_extremes = get_range_extremes(s);

  Point star_enter_range;
  if (angle_between_vectors(range_extremes[0], axis) > 0) {
    star_enter_range = range_extremes[0];
  } else {
    star_enter_range = range_extremes[1];
  }

  return abs_angle_between_vectors(s.point(), star_enter_range);
}

bool Probe::intersects_path_of(Star s) {
  return get_range_extremes(s).size() == 2;
}

/**
   Return a list of all the stars in the starfield that have at least the given
   magnitude.
**/
void Probe::get_backward_transfers(vector<Star> stars, double track_dist, double maglim) {
  for (Star s : stars) {
    if (s.r <= maglim) {
      backward_transfers.push_back(s);
    }
  }
}

bool member(Star s, vector<Star> list) {
  for (Star el : list) {
    if (s.equals(el)) {
      return true;
    }
  }

  return false;
}

/**
   Look for a star that will be in range after the given distance, and update current
   and base star attributes if one is found.
**/
int Probe::backtrack(double dist) {
  sort_by_transfer_distance(backward_transfers, 0, backward_transfers.size(), center, current_star);

  bool found = false;
  for (int i=0; i<backward_transfers.size(); i++) {
    Star s = backward_transfers[i];

    /* Only backtrack to stars we haven't used */
    if (member(s, used_transfers)) {
      continue;
    }

    /* Only backtrack if the star is in range */
    if (in_range(s.rotate(dist))) {
      found = true;
      base_star = s;
      current_star  = base_star.rotate(dist);

      used_transfers.push_back(s);

      break;
    }
  }

  if (found == false) {
    return -1;
  } else {
    return 0;
  }
}
