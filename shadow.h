double dotprod(Point a, Point b);
Point project(Polygon poly, Point axis);
int between(double x, Point interval);
int intersecting(Point u, Point v);
vector<Point> SAT_axes(Polygon poly);
int colliding(Polygon poly1, Polygon poly2);
void transform_and_print(StarGroup group, vector<Probe> probes);
int has_collisions_with_current_stars(vector<Probe> probes, Polygon obscuration);
int has_collisions_in_parts(StarGroup group, vector<Probe> probes, Polygon obscuration, int max_obscured);
int colliding_in_parts(vector<Polygon> probe1parts, vector<Polygon> probe2parts);
bool star_is_obscured(Star s, Polygon obscuration);
bool config_is_obscured(StarGroup group, vector<Probe> probes, Polygon obscuration, int max_obscured);

double shadowing(StarGroup group, vector<Probe> probes, double fpradius);
