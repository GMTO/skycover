#include "point.h"
#include "stargroup.h"
#include "probe.h"
#include "polygon.h"
#include <iostream>


double dotprod(Point a, Point b) {
  double sum = 0;
  sum += a.x * b.x;
  sum += a.y * b.y;

  return sum;
}

/**
   The projection of a polygon onto an axis is the minimum and maximum
   dotproduct of all of its points with the axis.
**/
Point project(Polygon poly, Point axis) {
  int i;
  Point projection;

  projection.x = dotprod(poly.points[0], axis);
  projection.y = projection.x;

  vector<Edge> edges = poly.edges();
  for (i=0; i<poly.npoints; i++) {
    projection.x = min(dotprod(edges[i].a, axis), projection.x);
    projection.x = min(dotprod(edges[i].b, axis), projection.x);
    
    projection.y = max(dotprod(edges[i].a, axis), projection.y);
    projection.y = max(dotprod(edges[i].b, axis), projection.y);
  }

  return projection;
}

int between(double x, Point interval) {
  double a = min(interval.x, interval.y);
  double b = max(interval.x, interval.y);

  return a <= x && x <= b;
}

int intersecting(Point u, Point v) {
  int res = 0;

  if (between(u.x, v)) {
    res = 1;
  } else if (between(u.y, v)) {
    res = 1;
  } else if (between(v.x, u)) {
    res = 1;
  } else if (between(v.y, u)) {
    res = 1;
  }

  return res;
}

/**
   The Separating Axis Theorem says that if two convex polygons are not
   intersecting, there is an axis perpendicular to one of the edges of one of the
   polygons, that when both polygons are projected onto this axis, the projections
   do not intersect.
   
   This function returns those axes perpendicular to the edges of a polygon.
**/
vector<Point> SAT_axes(Polygon poly) {
  vector<Edge> edges = poly.edges();
  vector<Point> axes;

  int i;
  for (i=0; i<poly.npoints; i++) {
    axes.push_back(edges[i].normal());
  }

  return axes;
}

/**
   Tests polygon collision by the separating axis theorem.
**/
int colliding(Polygon poly1, Polygon poly2) {
  Point poly1_projection, poly2_projection;
  
  vector<Point> axes = SAT_axes(poly1);
  vector<Point> poly2_axes = SAT_axes(poly2);

  vector<Point>::iterator it;
  for (it=poly2_axes.begin(); it!=poly2_axes.end(); it++) {
      axes.push_back(*it);
  }

  int res = 1;
  vector<Point>::iterator curr_axis;
  for (curr_axis=axes.begin(); curr_axis != axes.end(); ++curr_axis) {
    poly1_projection = project(poly1, *curr_axis);
    poly2_projection = project(poly2, *curr_axis);

    if (! intersecting(poly1_projection, poly2_projection)) {
      res = 0;
      break;
    }
  }

  return res;
}

/**
   This function tests for the intersection of two probes represented as multiple
   polygons. The function tests for collisions of each part of one probe with
   each part of the other.
**/
int colliding_in_parts(vector<Polygon> probe1parts, vector<Polygon> probe2parts) {
  for (Polygon probe1part : probe1parts) {
    for (Polygon probe2part : probe2parts) {
      if (colliding(probe1part, probe2part)) {
        return 1;
      }
    }
  }
  
  return 0;
}

bool probe_is_parked(Probe probe) {
    return false;
}


bool star_is_obscured(Star s, Polygon obscuration) {
  if (obscuration.point_in_poly(s.point())) {
    return true;
  }

  return false;
}

/**
   Determine whether or not the probes would be obscured given their current_star attribute.
**/
bool config_is_obscured_with_current_stars(vector<Probe> probes, Polygon obscuration) {
  for (int i=0; i<probes.size(); i++) {
      if (!probe_is_parked(probes[i]) && star_is_obscured(probes[i].current_star, obscuration)) {
      return true;
    }
  }

  return false;
}

/**
   Given a stargroup, test whether the probes would be obscured when transformed to
   intersect light from the stars in the stargroup.
**/
bool config_is_obscured(StarGroup group, vector<Probe> probes, Polygon obscuration, int max_obscured) {
  int obscured_probes = 0;
  for (int i=0; i<group.stars.size(); i++) {
    if (star_is_obscured(group.stars[i], obscuration)) {
      obscured_probes++;
    }
  }

  if (obscured_probes > max_obscured) {
    return true;
  }

  return false;
}

/**
   Determine whether a series of probes made up of multiple polygons each will have collisions or
   be obscured if transformed to intersect light from the stars in the given stargroup.
   
   - First test collisions between the first and last probes, then test all other pairs (1 and 2,
     2 and 3, 3 and 4).
**/
int has_collisions_in_parts(StarGroup group, vector<Probe> probes, Polygon obscuration, int max_obscured) {
  Probe probe1 = probes[0];
  Probe probe2 = probes[probes.size()-1];

  vector<Polygon> probe1parts = probe1.transform_parts(group.stars[0].point());
  vector<Polygon> probe2parts = probe2.transform_parts(group.stars[group.stars.size()-1].point());

  for (Polygon probe1part : probe1parts) {
    for (Polygon probe2part : probe2parts) {
      if (colliding(probe1part, probe2part)) {
        return 1;
      }
    }
  }

  if (!obscuration.points.empty()) {
    if (config_is_obscured(group, probes, obscuration, max_obscured)) {
      return true;
    }
  }

  for (int i=0; i<probes.size()-1; i++) {
    probe1 = probes[i];
    probe2 = probes[i+1];

    probe1parts = probe1.transform_parts(group.stars[i].point());
    probe2parts = probe2.transform_parts(group.stars[i+1].point());

    for (Polygon probe1part : probe1parts) {
      for (Polygon probe2part : probe2parts) {
        if (colliding(probe1part, probe2part)) {
          return 1;
        }
      }
    }
  }

  return 0;
}

/**
   Perform the same test as has_collisions_in_parts, but this time use the probes' current_star
   attribute.
**/
int has_collisions_with_current_stars(vector<Probe> probes, Polygon obscuration) {
  Probe probe1 = probes[0];
  Probe probe2 = probes[probes.size()-1];

  vector<Polygon> probe1parts = probe1.transform_parts(probes[0].current_star.point());
  vector<Polygon> probe2parts = probe2.transform_parts(probes[probes.size()-1].current_star.point());

  for (Polygon probe1part : probe1parts) {
    for (Polygon probe2part : probe2parts) {
      if (colliding(probe1part, probe2part)) {
        // cout << "collision between first and last" << endl;
        return 1;
      }
    }
  }

  if (!obscuration.points.empty()) {
    if (config_is_obscured_with_current_stars(probes, obscuration)) {
      return true;
    }
  }

  for (int i=0; i<probes.size()-1; i++) {
    probe1 = probes[i];
    probe2 = probes[i+1];

    probe1parts = probe1.transform_parts(probes[i].current_star.point());
    probe2parts = probe2.transform_parts(probes[i+1].current_star.point());

    for (Polygon probe1part : probe1parts) {
      for (Polygon probe2part : probe2parts) {
        if (colliding(probe1part, probe2part)) {
          return 1;
        }
      }
    }
  }

  return 0;
}
