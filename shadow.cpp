#include "point.h"
#include "stargroup.h"
#include "probe.h"
#include "polygon.h"
#include <iostream>

#include <boost/geometry.hpp>
#include <boost/foreach.hpp>

namespace bg = boost::geometry;

/** We use the "boost" package to compute the intersecting area of the probe shadow with the field of view;
    ideally we would use the same polygon type definition everywhere;
    maybe someday we will use boost for everything
**/

#define NSEGS 100

bool shadowing_is_less_than (StarGroup group, vector<Probe> probeshadows, double fovradius, double maxshadow) {

    typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
    typedef bg::model::polygon<point_t> polygon_t; 

    polygon_t fov;

    // The outer edge of the field of view
    for (int i = 0; i<=NSEGS; i++ ){
	double x = fovradius * cos(i * 2 * M_PI / NSEGS);
	double y = -fovradius * sin(i * 2 * M_PI / NSEGS); // Polygon must go clockwise for intersection to work
	bg::append(fov, point_t(x,y));
    }
    //    std::cout << "b" << bg::dsv(fov) << std::endl;

    
    int p  = 0;
    double a = 0;
    int probecount = 0;

    // For each probe
    for (Probe shadow : probeshadows ){
	probecount++;
	vector <Polygon> shadowparts = shadow.transform_parts(group.stars[p++].point());
	// For each section of the probe
	int partcount = 0;
	vector<polygon_t> tempresult;
	polygon_t result;
	vector<polygon_t> shadowpolyall;
	for (Polygon shadowpart : shadowparts) {
	    partcount++;
	    vector<Point>::iterator currpt;

	    polygon_t shadowpoly;
	    // Convert to boost data type
	    for (currpt=shadowpart.points.begin(); currpt!=shadowpart.points.end(); currpt++) {
		bg::append(shadowpoly, point_t(currpt->x, currpt->y));
	    }
	    bg::append(shadowpoly, point_t(shadowpart.points.begin()->x, shadowpart.points.begin()->y));
	    //	    std::cout << "y" << bg::dsv(shadowpoly) << std::endl;
	    shadowpolyall.push_back(shadowpoly);
	}

	// Merge all the three probe parts into a single polygon
	// Do them in this order so that adjacent parts are merged
	bg::union_(shadowpolyall.at(1), shadowpolyall.at(2), tempresult);
	result = tempresult.at(0);
	tempresult.clear();
	bg::union_(shadowpolyall.at(0), result,  tempresult);
	result = tempresult.at(0);
	//	std::cout << "r" << bg::dsv(result) << std::endl;


	// Compute the intersection of the probe with the field of view
	std::deque<polygon_t> output;
	bg::intersection(fov,result, output);
	int i=0;
	BOOST_FOREACH(polygon_t const& p, output)
	    {
		//		std::cout << "g" << bg::dsv(p) << std::endl;
		
		double area = bg::area(p);
		//		std::cout << i++ << ": " << area << std::endl;
		a += area;
	    }
    }

    double fraction = a / (M_PI * fovradius * fovradius);
    //    std::cout << "Frac: " << fraction << std::endl;
    return fraction < maxshadow;
}

