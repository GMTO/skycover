#include "star.h"
#include "stargroup.h"
#include "probe.h"
#include "prod.h"
#include "probe.h"
#include "collisions.h"
#include "shadow.h"
#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
using namespace std;


#define PI                      3.1415926535
#define MINIMUMTRACK            (60 * (PI / 180))
#define MINIMUMTRACK_DEG        60
#define N_OK_OBSCRD_FOR_PHASING 1
#define N_OK_OBSCRD_FOR_4PROBE  0

double fieldradius;
double maxshadow = 0.15;

bool PRINT = false;
bool SILENT = false;

bool phasing;

enum Mode { ModeGCLEF, ModeM3, ModeDGNF, ModeDGWF, ModeGMACS };
Mode mode;

/**
   Params:  vector of Star vectors.
   Returns: vector of ints denoting the size of each given Star vector.
   
   This function is used to get the sizes of the Star lists passed to
   the CombinationGenerator, which generates combinations of stars the
   probes can move to.
 **/
vector<int> get_list_sizes(vector< vector<Star> > lists) {
    int i;
    vector<int> result;

    for (i=0; i<lists.size(); i++) {
        result.push_back(lists[i].size());
    }

    return result;
}

/**
   Filter a list of stars for those that are in range of the given probe.
**/
vector<Star> in_probe_range(vector<Star> stars, Probe probe, double probe_angle) {
  vector<Star> result;
  int i;

  for (i=0; i<stars.size(); i++) {
    if (safe_distance_from_center(stars[i]) && probe.in_range(stars[i])) {
      stars[i].cablemax = (probe_angle + 90 + stars[i].bear) / 2;
      stars[i].cablemin = min(probe_angle, stars[i].bear);
      result.push_back(stars[i]);
    }
  }

  return result;
}

/**
   Partition function for the quicksort used in sorting stars.
   Sorts by angle from the vector <0, 1>.
**/
int partition(vector<Star>& stars, int p, int q)
{
  Point origin(0, 1);
    Star x = stars[p];
    int i = p;
    int j;

    for(j=p+1; j<q; j++)
    {
      if(angle_between_vectors(origin, stars[j].point()) >= angle_between_vectors(origin, x.point()))
        {
          i = i+1;
          swap(stars[i], stars[j]);
        }

    }

    swap(stars[i], stars[p]);

    return i;
}

/**
   Implementation of quicksort on the Star type.
**/
void star_sort(vector<Star>& stars, int p, int q)
{
    int r;
    if(p < q)
    {
        r = partition(stars, p, q);
        star_sort(stars, p, r);  
        star_sort(stars, r+1, q);
    }
}

/**
   Populate a list of Star lists containing stars that are in range of the probes.
   The first Star list in the returned list contains Stars in range of the first probe,
   and so on.
**/
vector< vector<Star> > get_probe_stars(vector<Star> stars, vector<Probe> probes) {
    int i;
    vector< vector<Star> > result;

    vector<double> probe = { 72, 144, -144, -72 };

    for (i=0; i<probes.size(); i++) {
      result.push_back(in_probe_range(stars, probes[i], probe[i]));
    }

    for (i=0; i<result.size(); i++) {
      star_sort(result[i], 0, result[i].size());
    }

    return result;
}

/**
   Given lists of stars within each probe's reach, and a list of indices, return
   the stars at the given indexes into the star lists.
**/
vector<Star> apply_indices(vector< vector<Star> > probestars, vector<int> indices) {
    int i;
    vector<Star> result;

    for (i=0; i<indices.size(); i++) {
        result.push_back(probestars[i][indices[i]]);
    }

    return result;
}

/**
   Print all the polygons in the system given by each probe positioned over its
   current_star attribute.
**/
void print_with_current_stars(vector<Probe> probes, Polygon obscuration) {
  for (int i=0; i<probes.size(); i++) {
    vector<Polygon> transformed_parts = probes[i].transform_parts(probes[i].current_star.point());
    for (Polygon part : transformed_parts) {
      part.polyprint();
    }
  }

  if (!obscuration.points.empty()) {
    obscuration.polyprint();
  }
}

/**
   Determine whether or not a configuration of probes can be tracked for the full 60 degrees.

   - Set each probe's current_star and base_star to the corresponding Star in the group.
   - Open up a loop that checks for collisions and necessary backtracks in one degree increments
     for the full 60 degrees.

     - At the top of the loop, ask each probe to track for the current number of degrees.
     - Clear each probe's used transfers list because backtracks may change at every degree of
       rotation.
     - If a probe cannot reach any star at the current rotation, return false.

     - Open up a loop that will run until either a configuration is found with no collisions/obscurations,
       or we are unable to find any configuration with no collisions/obscurations.
       
         - Check for collisions between probes 1 and 4. If there is a collision, ask probe 4 to backtrack.
         - If probe 4 is unable to backtrack, meaning there is no other star in its range, set the
           'nobacktrack' flag to true and break.
        
         - Loop through each of the other pairs of probes (1 and 2, 2 and 3, 3 and 4) and perform the same
           checks as those performed for the first and last probes.

         - If there is an obscuration (meaning not DGNF configuration), check if any of the probes are
           obscured in their current location. If so, set 'nobacktrack' and break.
           
         - If 'nobacktrack' was set in any of our loops, break because we hit a spot where backtrack was
           needed but couldn't be found.
           
   - Return whether or not we were able to do a full track.
**/
bool trackable(vector<Probe> probes, StarGroup group, int wfsmag, int gdrmag, Polygon obscuration) {
  for (int k=0; k<probes.size(); k++) {
    probes[k].current_star = group.stars[k];
    probes[k].base_star    = group.stars[k];
    probes[k].backward_transfer_idx = 0;
  }

  bool nobacktrack = false;

  for (int i=0; i<MINIMUMTRACK_DEG; i++) { /* one degree increments */
    for (int j=0; j<probes.size(); j++) {
      // probes[j].track(i * (PI / 180));
      probes[j].used_transfers.clear();
      if (probes[j].track(i * (PI / 180)) == -1) {  /* update the positions of the probes current star */
        return false;
      }
    }


    while (has_collisions_with_current_stars(probes, obscuration)) {  /* check current config as a whole */

	/**
	   if config fails, then break it down and check by pairs, backtracking a probe if necessary 
	**/

      for (int k=0; k<probes.size(); k++) {
        if (colliding_in_parts(probes[k].transform_parts(probes[k].current_star.point()),
                               probes[(k+1)%probes.size()].transform_parts(probes[(k+1)%probes.size()].current_star.point()))) {
          if (probes[k].backtrack(i * (PI / 180)) == -1) {
	    /* Didn't find a star to backtrack to */
            nobacktrack = true;
            break;
          }
        }
      }

      if (!obscuration.points.empty()) {
        for (int l=0; l<probes.size(); l++) {
          if (star_is_obscured(probes[l].current_star, obscuration)) {
            if (probes[l].backtrack(i * (PI / 180)) == -1) {
              nobacktrack = true;
              break;
            }
          }
        }
      }

      if (nobacktrack) { break; }
    }

    /* Add a nobacktrack break check here? */
  }

  if (nobacktrack) {
    return false;
  } else {
    return true;
  }
}

/**
   If any of the probes in probes list cannot track their star for the full 60 degrees,
   fill up their list of possible backward transfers. This list is all stars with magnitude
   at least as bright as the dimmest magnitude between wfsmag and gdrmag.
**/
void populate_backward_transfers(vector<Probe> &probes, StarGroup group, vector<Star> stars, int wfsmag, int gdrmag) {
  for (int i=0; i<probes.size(); i++) {
    double track_dist = probes[i].track_distance(group.stars[i]);
    if (track_dist < MINIMUMTRACK) {
      probes[i].needs_transfer = true;
      probes[i].get_backward_transfers(stars, track_dist, max(wfsmag, gdrmag));
    }
  }
}

/**
   Print all polygons in the system consisting of the probes transformed over their
   corresponding stars in the StarGroup. Also print the obscuration if a non-null polygon
   is given for it.
**/
void transform_and_print(vector<Probe> probes, StarGroup group, Polygon obscuration) {
  for (int i=0; i<probes.size(); i++) {
    vector<Polygon> transformed_parts = probes[i].transform_parts(group.stars[i].point());
    for (Polygon part : transformed_parts) {
      part.polyprint();
    }
  }

  if (!obscuration.points.empty()) {
    obscuration.polyprint();
  }
}

/**
   Perform the same logic as the function 'trackable', but print out the system at
   every rotation in one degree increments.
**/
void track_and_print_probes(vector<Probe> probes, StarGroup group, Polygon obscuration) {
  for (int k=0; k<probes.size(); k++) {
    probes[k].current_star = group.stars[k];
    probes[k].base_star    = group.stars[k];
    probes[k].backward_transfer_idx = 0;
  }

  bool nobacktrack = false;

  for (int i=0; i<MINIMUMTRACK_DEG; i++) {
    for (int j=0; j<probes.size(); j++) {
      probes[j].track(i * (PI / 180));
      probes[j].used_transfers.clear();
    }

    while (has_collisions_with_current_stars(probes, obscuration)) {
      if (colliding_in_parts(probes[0].transform_parts(probes[0].current_star.point()),
                             probes[probes.size()-1].transform_parts(probes[probes.size()-1].current_star.point()))) {
        if (probes[probes.size()-1].backtrack(i * (PI / 180)) == -1) {
          nobacktrack = true;
          break;
        }
      }

      for (int k=0; k<probes.size()-1; k++) {

        if (colliding_in_parts(probes[k].transform_parts(probes[k].current_star.point()),
                               probes[k+1].transform_parts(probes[k+1].current_star.point()))) {

          if (probes[k].backtrack(i * (PI / 180)) == -1) {
            nobacktrack = true;
            break;
          }
        }
      }

      if (!obscuration.points.empty()) {
        for (int l=0; l<probes.size(); l++) {
          if (star_is_obscured(probes[l].current_star, obscuration)) {
            if (probes[l].backtrack(i * (PI / 180)) == -1) {
              nobacktrack = true;
              break;
            }
          }
        }
      }

      if (nobacktrack) { break; }
    }

    print_with_current_stars(probes, obscuration);
  }
}

/**
   Determine whether the given star field contains a valid no-tracking configuration for the
   given wavefront sensor/guide star magnitude pair.
   
   - Create a CombinationGenerator which will hand back index lists giving the next combination
     of stars for the probes to move to.
     
   - Open up a loop that will run until all combinations have been checked, or a valid combination
     has been found.
     
     - At the top of the loop get the next combination of stars.
     - Check whether those stars satisfy the 3-wavefront senor, 1-guide star magnitude constraint.
     - Check for any collisions/obscurations at that configuration.
     - If these tests pass, return true. Print the configuration if the PRINT flag is set.
**/
bool is_valid_pair_notracking(vector<Star> stars, vector< vector<Star> > probestars, vector<Probe> probes,
                              double wfsmag, double gdrmag, Polygon obscuration) {
  CombinationGenerator stargroups(get_list_sizes(probestars));

  while ( !stargroups.done ) {
    StarGroup current_group = StarGroup(apply_indices(probestars, stargroups.next()));
    
    if (current_group.valid(wfsmag, gdrmag)) {

	if ( !has_collisions_in_parts(current_group, probes, obscuration, N_OK_OBSCRD_FOR_4PROBE) && shadowing(current_group, probes, fieldradius) < maxshadow ) {
        if (PRINT) {
          transform_and_print(probes, current_group, obscuration); 
        }
        return true;
      }
    }
  }

  return false;
}

/**
   Determine whether the given star field contains a valid tracking configuration for the given
   wavefront sensor/guide star magnitude pair.
   
   - Create a CombinationGenerator which will hand back index lists giving the next combination
     of stars for the probes to move to.
     
   - Open up a loop that will run until all combinations have been checked, or a valid combination
     has been found.
     
     - At the top of the loop reset each probe's needs_transfer flag to false. If a probe had previously
       needed a backward transfer, it would have done so and now needs to be checked again given its
       new position.
       
     - Ask the combination generator for the next configuration of stars.
     - Setup the probes for testing the next configuration by clearing their backtrack lists and setting
       their current_star to reflect the new group.
       
     - Check if the current group of stars has valid magnitudes (3 wfs and 1 gdr for regular. 3 wfs for
       phasing).
     - Check for any collisions/obscurations in the initial configuration.
     - Populate the backward transfer lists for any of the probes that will need it.
     
     - If a probe needs a transfer but has no stars in its transfer list, the configuration won't work,
       so we go on to the next one.
       
     - Call the trackable() function which traces the configuration through the full 60 degrees of tracking
       and determines whether the configuration is valid.
   
**/
bool is_valid_pair_tracking(vector<Star> stars, vector< vector<Star> > probestars, vector<Probe> probes,
                            double wfsmag, double gdrmag, Polygon obscuration) {

    CombinationGenerator stargroups(get_list_sizes(probestars));
  
    while ( !stargroups.done ) {
      for (Probe p : probes) {
        p.needs_transfer = false;
      }

      StarGroup current_group = StarGroup(apply_indices(probestars, stargroups.next()));

      for (int i=0; i<probes.size(); i++) {
        probes[i].backward_transfers.clear();
        probes[i].current_star = current_group.stars[i];
      }

      if (current_group.valid(wfsmag, gdrmag)) {
        if ( !has_collisions_in_parts(current_group, probes, obscuration, N_OK_OBSCRD_FOR_4PROBE) ) {
          populate_backward_transfers(probes, current_group, stars, wfsmag, gdrmag);

          bool nobacktrack = false;
          for (Probe p : probes) {
            if (p.needs_transfer && p.backward_transfers.size() == 0) {
              nobacktrack = true;
              break;
            }
          }

          if (nobacktrack) {
            continue;
          }

          if (trackable(probes, current_group, wfsmag, gdrmag, obscuration)) {
            if (PRINT) {
              track_and_print_probes(probes, current_group, obscuration); 
            }
            return true;
          }
        }
      }
    }

    return false;
}

/**
   Determines whether or not the given wfsmag has a valid probe configuration in the given starfield.

   - The flag N_OK_OBSCRD_FOR_PHASING is the number of probes that are ok to be obscured. Since
     the phasing calculation only cares about having three stars, this flag is 1.
**/
bool is_valid_phasing_mag(vector< vector<Star> > probestars, vector<Probe> probes, double maglim, Polygon M3) {
  for (int i=0; i<probestars.size(); i++) {
    probestars[i].push_back(probes[i].default_star);
  }

  CombinationGenerator stargroups(get_list_sizes(probestars));

  while ( !stargroups.done ) {
    StarGroup current_group = StarGroup(apply_indices(probestars, stargroups.next()));

    if (current_group.valid_for_phasing(maglim)) {
      if ( !has_collisions_in_parts(current_group, probes, M3, N_OK_OBSCRD_FOR_PHASING) ) {
        return true;
      }
    }
  }

  return false;
}

/**
   Read a starbase star catalogue file and populate a list Star objects. This list represents
   a single starfield.
   
   - Reading starts after line 2 because the first two lines are the header and the dashed line.
**/
vector<Star> load_stars(string filename) {
    vector<Star> stars;
    vector<string> tokens;
    std::ifstream infile(filename);

    int current_line = 1;
    string line; getline(infile, line);
    tokens = split(line, '\t');

    int idx = find(tokens.begin(), tokens.end(), "x")    - tokens.begin();
    int idy = find(tokens.begin(), tokens.end(), "y")    - tokens.begin();
    int idr = find(tokens.begin(), tokens.end(), "R")    - tokens.begin();
    int idb = find(tokens.begin(), tokens.end(), "Bear") - tokens.begin();
    int idj = find(tokens.begin(), tokens.end(), "J")    - tokens.begin();
    int idrj= find(tokens.begin(), tokens.end(), "RJ")   - tokens.begin();
    int siz = tokens.size();

    for ( ; getline(infile, line); ) {
        if (current_line > 2) {
            tokens = split(line, '\t');
            double x = stod(tokens[idx]) * 3600;
            double y = stod(tokens[idy]) * 3600;
            double r = stod(tokens[idr]);
	    double j;
	    if (idj<siz) {
		j = stod(tokens[idj]);
	    } else if (idrj<siz) {
		j = r - stod(tokens[idrj]);
	    } else if (phasing) {
		printf("No column name J or RJ found.  Exiting");
		exit(1);
	    }
            double bear = stod(tokens[idb]);

            stars.push_back(Star(x, y, r, j, bear));
        }
        current_line++;
    }

    return stars;
}

/**
   This function filters the four lists of Stars (one for each probe) for stars that
   are at least as bright as the given magnitude. This is done to narrow down the number
   of calculations that have to be done for a given starfield.
   
   - The 'bin' param is the magnitude.
**/
vector< vector<Star> > probestars_in_bin(vector< vector<Star> > probestars, int bin) {
  vector< vector<Star> > result;
  vector<Star> probeX_stars;

  for ( vector<Star> stars : probestars ) {
    probeX_stars.clear();
    for ( Star s : stars ) {
	if ((phasing && s.j <= bin) || (!phasing && s.r <=bin)) {
	    probeX_stars.push_back(s);
	}
    }
    result.push_back(probeX_stars);
  }

  return result;
}

/**
   Return a list of filenames in the given directory.   

   - This function is used to read the names of star catalogues.
**/
vector<string> files_in_dir(string dirname) {
  ostringstream path;
  vector<string> filenames;
  DIR *pDIR;
  struct dirent *entry;
  if( (pDIR = opendir(dirname.c_str())) != NULL ) {
    while( (entry = readdir(pDIR)) != NULL ) {
      if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 ) {
        path.str(std::string());
        path << dirname << entry->d_name;
        filenames.push_back(path.str());
      }
    }
    closedir(pDIR);
  }

  return filenames;
}

/**
   Write the coordinates of all stars in a list to a file. This function is used to record
   a starfield when a valid config is found.
**/
void write_stars(vector<Star> stars, string filename, int wfsmag, int gdrmag) {
  std::ofstream fout;
  fout.open(filename);

  for (Star s : stars) {
    if (s.r < max(wfsmag, gdrmag)) {
      fout << s.x << " " << s.y << endl;
    }
  }

  fout.close();
}

/**
   This function controls the iteration over multiple files, testing each one and recording
   whether or not the starfield given in the file supports a valid configuration for the given
   phasing magnitude.
   
   - When a valid file is found, the starfield in the file is written to a file in the 'starfiles'
     directory. This info is recorded so that the visual simulation is able to display starfields
     over the probes.
**/
int number_valid_phasing_files(vector<string> starfld_files, vector<Probe> probes,
                               double maglim, int nfiles, Polygon M3) {
  int valid_files = 0;

  for (int i=0; i<nfiles; i++) {
    vector<Star> stars = load_stars(starfld_files[i]);

    vector< vector<Star> > probestars  = get_probe_stars(stars, probes);
    vector< vector<Star> > current_bin = probestars_in_bin(probestars, ceil(maglim));

    if (is_valid_phasing_mag(current_bin, probes, maglim, M3)) {
      valid_files++;

      if (PRINT) {
        ostringstream starfile;
        starfile << "starfiles/starfield" << valid_files << ".cat";
        write_stars(stars, starfile.str(), maglim, maglim);
      }
    } else {
      ostringstream starfile;
      // starfile << "starfiles/starfield" << valid_files << "_invalid.cat";
      // write_stars(stars, starfile.str(), maglim, maglim);
    }
  }
  
  return valid_files;
}

/**
   This function controls the iteration over multiple files, testing each one and recording
   whether or not the starfield given in the file supports a valid configuration for the given
   magnitude pair.

   - The first thing this function does is choose to run the tracking or non-tracking test. The functions
     used for both those tests have the same api, and therefore a function pointer can be switched depending
     on the 'tracking' parameter that is passed in.
     
   - When a valid file is found, the starfield in the file is written to a file in the 'starfiles'
     directory. This info is recorded so that the visual simulation is able to display starfields
     over the probes.
**/
int number_valid_4probe_files(vector<string> starfld_files, vector<Probe> probes, double wfsmag, double gdrmag,
                              int nfiles, Polygon obscuration, bool tracking) {
  int valid_files = 0;

  bool (*is_valid_pair)(vector<Star> stars, vector< vector<Star> > probestars, vector<Probe> probes,
                        double wfsmag, double gdrmag, Polygon obscuration);

  is_valid_pair = is_valid_pair_notracking;
  if (tracking) {
    is_valid_pair = is_valid_pair_tracking;
  }
  
  for (int i=0; i<nfiles; i++) {
    vector<Star> stars = load_stars(starfld_files[i]);

    vector< vector<Star> > probestars  = get_probe_stars(stars, probes);
    vector< vector<Star> > current_bin = probestars_in_bin(probestars, ceil(max(wfsmag, gdrmag)));

    if (is_valid_pair(stars, current_bin, probes, wfsmag, gdrmag, obscuration)) {
      valid_files++;

      if (PRINT) {
        ostringstream starfile;
        starfile << "starfiles/starfield" << valid_files << ".cat";
        write_stars(stars, starfile.str(), wfsmag, gdrmag);
      }
    } else {
      ostringstream starfile;
      // starfile << "starfiles/starfield" << valid_files << "_invalid.cat";
      // write_stars(stars, starfile.str(), wfsmag, gdrmag);
    }
  }
  
  return valid_files;
}

/**
   Load the polygon that represents the m3 obscuration
**/
Polygon get_m3_obscuration() {
  return load_poly("m3_obsc.txt");
}

/**
   Load the polygon that represents the gclef obscuration
**/
Polygon get_gclef_obscuration() {
  return load_poly("gclef_obsc.txt");
}

/**
   Load the polygon that represents the gmacs obscuration
**/
Polygon get_gmacs_obscuration() {
  return load_poly("gmacs_obsc.txt");
}

int main(int argc, char *argv[]) {
  string probe_slider_body_file  = "probe_slider_body.txt";
  string probe_slider_shaft_file = "probe_slider_shaft.txt";
  string probe_baffle_tube_file  = "probe_baffle_tube.txt";

  Probe probe1(0,   probe_slider_body_file, probe_slider_shaft_file, probe_baffle_tube_file);
  Probe probe2(90,  probe_slider_body_file, probe_slider_shaft_file, probe_baffle_tube_file);
  Probe probe3(180, probe_slider_body_file, probe_slider_shaft_file, probe_baffle_tube_file);
  Probe probe4(270, probe_slider_body_file, probe_slider_shaft_file, probe_baffle_tube_file);

  vector<Probe> probes;
  probes.push_back(probe1);
  probes.push_back(probe2);
  probes.push_back(probe3);
  probes.push_back(probe4);

  string shadow_slider_body_file  = "shadow_slider_body.txt";
  string shadow_slider_shaft_file = "shadow_slider_shaft.txt";
  string shadow_baffle_tube_file  = "shadow_baffle_tube.txt";

  Probe shadow1(0,   shadow_slider_body_file, shadow_slider_shaft_file, shadow_baffle_tube_file);
  Probe shadow2(90,  shadow_slider_body_file, shadow_slider_shaft_file, shadow_baffle_tube_file);
  Probe shadow3(180, shadow_slider_body_file, shadow_slider_shaft_file, shadow_baffle_tube_file);
  Probe shadow4(270, shadow_slider_body_file, shadow_slider_shaft_file, shadow_baffle_tube_file);

  vector<Probe> shadows;
  shadows.push_back(shadow1);
  shadows.push_back(shadow2);
  shadows.push_back(shadow3);
  shadows.push_back(shadow4);

  Polygon GCLEF = get_gclef_obscuration();
  Polygon M3    = get_m3_obscuration();
  Polygon GMACS = get_gmacs_obscuration();
  Polygon DGNF;

  /**
  arg format:
      ./agwsvalid <--gclef | --m3 | --dgnf | --dgwf> <--plot | --bool>

  argv[1] -> obscuration type. '--dgnf' for no obscuration
  **/
  
  Polygon obscuration;
  double  x1, y1, x2, y2, x3, y3, x4, y4;
  Star    s1, s2, s3, s4;
  vector <Star> stars;
  char s[100];
  int n;
  int count=1;
  

  double scale;  // mm per degree (we ignore distortion)

  if (argc < 3) {
      cout << "usage: ./agwsvalid <--gclef | --m3 | --dgnf | --dgwf | --gmacs> <--plot | --bool>\nx1 y1 x2 y2 x3 y3 x4 y4\netc" << endl;
      return 0;
  } else {

      if (strcmp(argv[1], "--gclef") == 0) {
	  obscuration = GCLEF;
	  mode = ModeGCLEF;
	  scale = 3600 * 0.98716; 
      } else if (strcmp(argv[1], "--m3") == 0) {
	  obscuration = M3;
	  mode = ModeM3;
	  scale = 3600 * 0.98716;
      } else if (strcmp(argv[1], "--dgnf") == 0) {
	  obscuration = DGNF;
	  mode = ModeDGNF;
	  scale = 3600 * 0.98716;
      } else if (strcmp(argv[1], "--dgwf") == 0) {
	  obscuration = DGNF;
	  mode = ModeDGWF;
	  scale = 3600 * 1.04938;
      } else if (strcmp(argv[1], "--gmacs") == 0) {
        obscuration = GMACS;
        mode = ModeGMACS;
        scale = 3600 * 65.02345/60.; //mm per degree (high order terms ignored)
      } else {
	  cerr << "Unknown mode:" << argv[1];
	  exit(1);
      }
      
      fieldradius = 10 / 60. * scale;//in mm, for radius = 10 arcmin

      if (strcmp(argv[2], "--plot") == 0) {
	  PRINT = 1;
      } else if (strcmp(argv[2], "--boolonly") == 0) {
	  PRINT = 0;
	  SILENT = 1;
      }
      else {
	  PRINT = 0;
      }
      
  }

  while (!feof(stdin)) {
      vector <Star> stars;
      StarGroup g;
      int valid = 1;
      int i;
      double shadowfrac;
      
      if (!fgets(s,100,stdin)) continue;
      if (s[0] == '#') continue;

      n = sscanf(s, "%lf %lf %lf %lf %lf %lf %lf %lf", &x1, &y1, &x2, &y2 , &x3, &y3, &x4, &y4); 

      if (n < 8) break;

      // Convert from degrees to mm (approximately)
      stars.push_back( Star(x1*scale, y1*scale, 0, 0, 0));
      stars.push_back( Star(x2*scale, y2*scale, 0, 0, 0));
      stars.push_back( Star(x3*scale, y3*scale, 0, 0, 0));
      stars.push_back( Star(x4*scale, y4*scale, 0, 0, 0));
      g = StarGroup(stars);
      // g.print();


      // Check that each probe can reach its assigned star
      for (i=0; i<4;i++) {
	  if (!probes[i].in_range(g.star_at(i))) {
	      if ( !SILENT) {
		  cerr <<  "Probe " << i << " cannot reach star " << i << ".\n" ; 
	      }
	      valid = 0;
	  }
      }

      // Check that probes are not obscured
      if (valid) {
	  if ( config_is_obscured(g, probes, obscuration, N_OK_OBSCRD_FOR_4PROBE)) {
	      if ( !SILENT) {
		  cerr << "At least one star is obscured.\n";
	      }
	      valid = 0;
	  }
      }


      // Check that probes do not get too close
      if (valid) {
	  if ( has_collisions_in_parts(g, probes, obscuration, N_OK_OBSCRD_FOR_4PROBE) ) {
	      if ( !SILENT) {
		  cerr << "Probes too close.\n";
	      }
	      valid = 0;
	  }
      }

      // Check in DGWF mode that shadowing is not too large
      shadowfrac = 0;
      if (valid && mode==ModeDGWF) {
	  shadowfrac=shadowing(g, shadows, fieldradius);
          //cerr << "shadow: " << shadowfrac << " " << maxshadow << "\n";
	  if (  shadowfrac > maxshadow ) {
	      if ( !SILENT) {
		  cerr << "Shadowing too large.\n";
	      }
	      valid = 0;
	  }
      }

      if (valid && !SILENT) cerr << "Looks good!\n";
      if (PRINT) {
	  transform_and_print(probes, g, obscuration); 
	  ostringstream starfile;
	  starfile << "starfiles/starfield" << count++ << ".cat";
	  write_stars(stars, starfile.str(), 1,1);
	  fflush(stdout);
      } else {
	  printf("%d %.3f\n", valid, shadowfrac);
	  fflush(stdout);
      }
  }

  return 0;
}
