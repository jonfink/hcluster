#ifndef _HCLUSTER_H_
#define _HCLUSTER_H_

#include <math.h>

#include<map>
#include<list>
#include<vector>
using namespace std;

struct cPoint_t
{
  int index;               // ordered index of this point
  vector<double> point;    // position of the point
  double data;             // measurement
  double min_dist_to_set;  // R_index
  int parent_index;        // parent index
  list<int> children;      // children of this node
};

void hcluster_points(list<cPoint_t*> &input_points, map<int, cPoint_t*> &output_points);
void compute_n_clusters(int cluster_count, map<int, cPoint_t*> &output_points,
                        vector<list<cPoint_t*> > &clusters);
#endif
