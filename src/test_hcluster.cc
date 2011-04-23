#include "hcluster.h"

#include<map>
#include<list>
#include<vector>
using namespace std;

int main(int argc, char **argv)
{
  list<cPoint_t*> input_points;
  map<int, cPoint_t*> output_points;

  // Add some random points
  for(int i=0; i < 9; ++i) {
    input_points.push_back(new cPoint_t());
    input_points.back()->point.resize(2);
    input_points.back()->point[0] = (rand() % 1000)/1000.0;
    input_points.back()->point[1] = (rand() % 1000)/1000.0;
    //input_points.back()->max_dist_to_set = 0.0;
    input_points.back()->min_dist_to_set = 1000;
  }

  // print input points
  printf("Inputs:\n");
  printf("-------\n");
  for(list<cPoint_t*>::iterator iter = input_points.begin();
      iter != input_points.end();
      ++iter) {
    printf("  %2.2f %2.2f\n", (*iter)->point[0], (*iter)->point[1]);
  }
  printf("\n");

  hcluster_points(input_points, output_points);

  // print output points
  printf("Outputs:\n");
  printf("--------\n");
  for(map<int, cPoint_t*>::iterator iter = output_points.begin();
      iter != output_points.end();
      ++iter) {
    printf("  [%d] %2.2f %2.2f, R=%2.2f, parent: %d, ",
           iter->first,
           iter->second->point[0], iter->second->point[1],
           iter->second->min_dist_to_set,
           iter->second->parent_index);
    printf("children: ");
    for(list<int>::iterator citer = iter->second->children.begin();
        citer != iter->second->children.end();
        ++citer) {
      printf("%d ", *citer);
    }
    printf("\n");
  }
  printf("\n");

  vector<list<cPoint_t*> > clusters;

  compute_n_clusters(4, output_points, clusters);

  for(unsigned int i=0; i < clusters.size(); ++i) {
    printf("Cluster %d (%d): ", i, (int)clusters[i].size());
    for(list<cPoint_t*>::iterator iter = clusters[i].begin();
        iter != clusters[i].end(); ++iter) {
      printf("%d ", (*iter)->index);
    }
    printf("\n");
  }

  return 0;
}
