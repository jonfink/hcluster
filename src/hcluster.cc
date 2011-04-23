#include<map>
#include<list>
#include<vector>
using namespace std;

#include "hcluster.h"

double norm(vector<double> &a, vector<double> &b)
{
  if(a.size() != b.size()) {
    fprintf(stderr, "norm: vectors must be the same size!\n");
    return -1;
  }

  double val = 0;
  for(unsigned int i = 0; i < a.size(); ++i) {
    val += (a[i] - b[i])*(a[i] - b[i]);
  }

  return sqrt(val);
}

void hcluster_points(list<cPoint_t*> &input_points, map<int, cPoint_t*> &output_points)
{
  int current_index = 1;
  double R=-1;
  double beta = 2.0;
  double alpha = 1.0;
  vector<list<cPoint_t*> > gLevels;
  vector<double> gLevelsMin, gLevelsMax;

  // move point from input to output and label as index 1
  output_points.insert(make_pair(current_index, input_points.front()));
  input_points.pop_front();
  output_points[current_index]->index = current_index;
  //output_points[current_index]->max_dist_to_set = 0.0;
  output_points[current_index]->min_dist_to_set = 0.0;
  output_points[current_index]->parent_index = -1;

  // Introduce L0 (will only hold point 1)
  gLevels.push_back(list<cPoint_t*>());
  gLevels.back().push_back(output_points[current_index]);
  gLevelsMin.push_back(0);
  gLevelsMax.push_back(0);

  while(input_points.size() > 0) {
    list<cPoint_t*>::iterator next_point_iter;
    double max_distance = 0;
    double dist;
    for(list<cPoint_t*>::iterator iter = input_points.begin();
        iter != input_points.end();
        ++iter) {
      // check if latest distance to output set is a new min
      dist = norm((*iter)->point, output_points[current_index]->point);

      if(dist < (*iter)->min_dist_to_set) {
        (*iter)->min_dist_to_set = dist;
      }

      // this could be our next point!
      if((*iter)->min_dist_to_set > max_distance) {
        next_point_iter = iter;
        max_distance = (*iter)->min_dist_to_set;
      }
    }

    // increment current_index and add our next point
    current_index++;
    output_points.insert(make_pair(current_index, *next_point_iter));
    output_points[current_index]->index = current_index;
    input_points.erase(next_point_iter);

    if(current_index == 2) {
      R = alpha*output_points[current_index]->min_dist_to_set;
    }

    // place this point in a granularity level
    int current_level = -1;
    for(unsigned int i = 1; i < gLevels.size(); ++i) {
      if(output_points[current_index]->min_dist_to_set > gLevelsMin[i] &&
         output_points[current_index]->min_dist_to_set <= gLevelsMax[i] ) {
        gLevels.back().push_back(output_points[current_index]);
        current_level = i;
      }
    }
    while(current_level < 0) {
      // create a new level
      gLevels.push_back(list<cPoint_t*>());
      if(gLevels.size() >= 2) {
        gLevelsMin.push_back(R/pow(beta, (double)(gLevels.size()-1.0)));
        gLevelsMax.push_back(R/pow(beta, (double)(gLevels.size()-2.0)));
      }
      else {
        gLevelsMin.push_back(R/beta);
        gLevelsMax.push_back(R);
      }

      if(output_points[current_index]->min_dist_to_set > gLevelsMin.back() &&
         output_points[current_index]->min_dist_to_set <= gLevelsMax.back() ) {
        gLevels.back().push_back(output_points[current_index]);
        current_level = gLevels.size()-1;
      }
    }

    // find parent (closest point in lower granularity level)
    double min_dist = 1000;
    int min_dist_index = 0;
    for(int i = 0; i < current_level; ++i) {
      for(list<cPoint_t*>::iterator iter = gLevels[i].begin();
          iter != gLevels[i].end(); ++iter) {
        dist = norm((*iter)->point, output_points[current_index]->point);
        if(dist < min_dist) {
          min_dist = dist;
          min_dist_index = (*iter)->index;
        }
      }
    }
    output_points[current_index]->parent_index = min_dist_index;

    // add this point to parent's child list
    output_points[output_points[current_index]->parent_index]->children.push_back(current_index);
  }

}

void add_to_cluster(cPoint_t* root_node,
                    map<int,cPoint_t*> &tree_data, int cluster_count,
                    list<cPoint_t*> &cluster)
{
  cluster.push_back(root_node);
  for(list<int>::iterator child_iter = root_node->children.begin();
      child_iter != root_node->children.end();
      ++child_iter)
    {
      if(*child_iter > cluster_count) {
        add_to_cluster(tree_data[(*child_iter)], tree_data, cluster_count, cluster);
      }
    }
}

void compute_n_clusters(int cluster_count, map<int, cPoint_t*> &output_points,
                        vector<list<cPoint_t*> > &clusters)
{
  clusters.resize(cluster_count);
  for(int i=1; i <= cluster_count; ++i) {
    add_to_cluster(output_points[i], output_points, cluster_count, clusters[i-1]);
  }
}
