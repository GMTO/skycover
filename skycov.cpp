#include "star.h"
#include "stargroup.h"
#include "probe.h"
#include "prod.h"
#include "probe.h"
#include "collisions.h"
#include "threadpool.h"
#include <time.h>
#include <regex>
#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <tablelib/table.h>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <numeric>
using namespace std;


#define NUMTHREADS 4
#define PI         3.1415926535
#define MINRANGE   0.1
#define MAXRANGE   0.167

map<string, bool> GlobalValidMags_CurrentFile;
mutex mtx;


vector<Star> apply_indices(vector< vector<Star> > probestars, vector<int> indices) {
    int i;
    vector<Star> result;

    for (i=0; i<indices.size(); i++) {
      if (indices[i] > probestars[i].size()) {
        cout << "bad index" << endl;
        for (int j=0; j<indices.size(); j++) {
          cout << indices[j] << " ";
        }
        cout << endl;
        for (int j=0; j<probestars.size(); j++) {
          cout << probestars[j].size() << " ";
        }
        cout << endl;
      }
      result.push_back(probestars[i][indices[i]]);
    }

    return result;
}

std::vector<std::string> split(const std::string &text, char sep) {
    std::vector<std::string> tokens;
    std::size_t start = 0, end = 0;
    while ((end = text.find(sep, start)) != std::string::npos) {
        tokens.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    tokens.push_back(text.substr(start));
    return tokens;
}

vector<string> dimmer_gdr(string magpair) {
    vector<string> dimmers;
    ostringstream pair;

    vector<string> mags = split(magpair, ':');
    int wfs = stoi(mags[0]);
    int gdr = stoi(mags[1]);

    for (int i=gdr+1; i<=19; i++) {
        pair << wfs << ":" << i;
        dimmers.push_back(pair.str());
        pair.str(std::string());
    }

    return dimmers;
}

string get_wfs(string magpair) {
    vector<string> mags = split(magpair, ':');
    return mags[0];
}

string get_gdr(string magpair) {
    vector<string> mags = split(magpair, ':');
    return mags[1];
}

vector<string> dimmer_wfs(string magpair) {
    vector<string> dimmers;
    ostringstream pair;

    vector<string> mags = split(magpair, ':');
    int wfs = stoi(mags[0]);
    int gdr = stoi(mags[1]);

    for (int i=wfs+1; i<=19; i++) {
        pair << i << ":" << gdr;
        dimmers.push_back(pair.str());
        pair.str(std::string());
    }

    return dimmers;
}

Worker::Worker(CombinationGenerator gen, vector<Probe> _probes, vector< vector<Star> > _probestars) {
  generator = gen;
  probes = _probes;
  probestars = _probestars;
}

void pmap(map<string, bool> m) {
  for (auto const entry : m) {
    cout << entry.first << ": " << entry.second << endl;
  }
}

void Worker::operator()() {
  string magpair;
  StarGroup current_group;

  while (! generator.done) {
    mtx.lock();
    // cout << "getting group" << endl;
    current_group = StarGroup(apply_indices(probestars, generator.next()));
    // cout << "got group" << endl;
    mtx.unlock();

    if ( current_group.valid() ) {
      if ( ! has_collisions(current_group, probes) ) {
        magpair = current_group.magpair();
        valid_mags_map[magpair] = true;

        for ( string dimmer_pair : dimmer_wfs(magpair) ) { valid_mags_map[dimmer_pair] = true; }
        for ( string dimmer_pair : dimmer_gdr(magpair) ) { valid_mags_map[dimmer_pair] = true; }

        // brightest_pair = brightest(magpair, brightest_pair);
      }
    }
  }

  mtx.lock();
  GlobalValidMags_CurrentFile.insert(valid_mags_map.begin(), valid_mags_map.end());
  mtx.unlock();

  // cout << "Done" << endl;
}

void ThreadPool::joinall() {
  for (int i=0; i<workers.size(); i++) {
    // cout << "joining workers[" << i << "]" << endl;
    workers[i].join();
  }
}

ThreadPool::ThreadPool(vector<CombinationGenerator> generators, vector<Probe> probes, vector< vector<Star> > probestars)
{
  for ( CombinationGenerator gen : generators ) {
    workers.push_back(thread(Worker(gen, probes, probestars)));
  }
}

ThreadPool::~ThreadPool() { }

vector<int> get_list_sizes(vector< vector<Star> > lists) {
    int i;
    vector<int> result;

    for (i=0; i<lists.size(); i++) {
        result.push_back(lists[i].size());
    }

    return result;
}

bool safe_distance_from_center(Star star) {
  Point star_pt(star.x, star.y), origin(0, 0);
  
  double dist = distance(star_pt, origin);
  
  if (MINRANGE < dist && dist < MAXRANGE) {
    return true;
  } else {
    return false;
  }
}

vector<Star> in_probe_range(vector<Star> stars, Probe probe) {
  Point star_pt;
  vector<Star> result;
  int i;

  for (i=0; i<stars.size(); i++) {
    star_pt = Point(stars[i].x, stars[i].y);
    if (safe_distance_from_center(stars[i]) && (probe.can_cover(stars[i]))) {
      result.push_back(stars[i]);
    }
  }

  return result;
}

vector< vector<Star> > get_probe_stars(vector<Star> stars, vector<Probe> probes) {
    int i;
    vector< vector<Star> > result;

    for (i=0; i<probes.size(); i++) {
        result.push_back(in_probe_range(stars, probes[i]));
    }

    return result;
}

int mul(int x, int y) { return x * y; }

void get_valid_mags(vector< vector<Star> > probestars, vector<Probe> probes) {
  int i;
  string current_pair;
  StarGroup current_group;
  map<string, bool> valid_mags_map;

  vector<int> list_sizes(get_list_sizes(probestars));
  int numcombos = accumulate(list_sizes.begin(), list_sizes.end(), 1, mul);
  int blocksize = numcombos / NUMTHREADS;

  vector<CombinationGenerator> generators;
  for (i=0; i<NUMTHREADS; i++) {
    generators.push_back(CombinationGenerator(list_sizes, i*blocksize, blocksize));
  }

  time_t starttime, endtime;
  time(&starttime);

  ThreadPool threads(generators, probes, probestars);
  threads.joinall();
  
  time(&endtime);
  cerr << "finding valid mags took " << difftime(endtime, starttime) << " seconds" << endl;
}

vector<Star> load_stars(string filename) {
    vector<Star> stars;
    vector<string> tokens;
    std::ifstream infile(filename);

    int current_line = 1;
    string line; getline(infile, line);
    for ( ; getline(infile, line); ) {
        if (current_line > 2) {
            tokens = split(line, '\t');
            double x = stod(tokens[0]);
            double y = stod(tokens[1]);
            double r = stod(tokens[18]);

            stars.push_back(Star(x, y, r));
        }
        current_line++;
    }

    return stars;
}

void valid_mags_in_starfield(vector<Star> stars, vector<Probe> probes) {
    get_valid_mags(get_probe_stars(stars, probes), probes);
}

vector<string> files_in_dir(string dirname, char *regexp) {
  ostringstream path;
  vector<string> filenames;
  DIR *pDIR;
  struct dirent *entry;
  if( (pDIR = opendir(dirname.c_str())) != NULL ) {
    while( (entry = readdir(pDIR)) != NULL ) {
      if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 ) {
        if (regex_match(entry->d_name, regex(regexp))) {
          path.str(std::string());
          path << dirname << entry->d_name;
          filenames.push_back(path.str());
        }
      }
    }
    closedir(pDIR);
  }

  return filenames;
}

void dogrid_file(map<string, int> ValidMags, double nfiles) {
    cout << "grid\tconfig\twfsmag\tgdrmag\tprob\n";
    cout << "----\t------\t------\t------\t----\n";

    string wfs, gdr;
    double prob;
    for (auto const entry : ValidMags) {
        wfs = get_wfs(entry.first);
        gdr = get_gdr(entry.first);

        prob = entry.second / (nfiles * 1.0);

        cout << "grid\tdgnf\t" << wfs << "\t" << gdr << "\t" << prob << "\n";
    }
}

int main(int argc, char *argv[]) {
    Point probe1_ctr(0.25, 0);
    Point probe2_ctr(0, 0.25);
    Point probe3_ctr(-0.25, 0);
    Point probe4_ctr(0, -0.25);

    double range_width = 0.4;

    Probe probe1(27, 94.5, 72);
    Probe probe2(121.5, 166.5, 144);
    Probe probe3(-27, -94.5, -72);
    Probe probe4(-121.5, -166.5, -144);

    Point origin(0, 0);

    vector<Probe> probes;
    probes.push_back(probe1);
    probes.push_back(probe2);
    probes.push_back(probe3);
    probes.push_back(probe4);

    vector<Star>  stars;

    map<string, int> ValidMagnitudes;

    vector<string> starfield_files = files_in_dir("Bes2/", argv[1]);
    vector<string>::iterator curr_path;
    for (curr_path=starfield_files.begin(); curr_path!=starfield_files.end(); curr_path++) {
      // if (count > 0) { break; }

        cerr << "Processing file " << *curr_path << endl;
        stars = load_stars(*curr_path);

        valid_mags_in_starfield(stars, probes);

        for (auto const entry : GlobalValidMags_CurrentFile) {
            if (ValidMagnitudes.count(entry.first) == 0) {
                ValidMagnitudes[entry.first] = 1;
            } else {
                ValidMagnitudes[entry.first]++;
            }
        }
        GlobalValidMags_CurrentFile.clear();
    }

    dogrid_file(ValidMagnitudes, starfield_files.size());

    return 0;
}
