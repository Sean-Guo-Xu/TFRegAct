#include <Rcpp.h>
#include <deque>
#include <utility>
#include <vector>

// [[Rcpp::export]]
bool recursive_bayes_ball_cpp(
    Rcpp::List parents,
    Rcpp::List children,
    int exposure_index,
    int outcome_index,
    Rcpp::IntegerVector blocked_nodes) {
  const int n = parents.size();
  const int exposure = exposure_index - 1;
  const int outcome = outcome_index - 1;
  std::vector<unsigned char> blocked(n, 0);
  std::vector<unsigned char> ancestor_of_blocked(n, 0);
  std::deque<int> ancestor_queue;

  for (int value : blocked_nodes) {
    const int node = value - 1;
    if (node >= 0 && node < n && !blocked[node]) {
      blocked[node] = 1;
      ancestor_of_blocked[node] = 1;
      ancestor_queue.push_back(node);
    }
  }
  if (blocked[exposure] || blocked[outcome]) {
    return false;
  }

  while (!ancestor_queue.empty()) {
    const int current = ancestor_queue.front();
    ancestor_queue.pop_front();
    Rcpp::IntegerVector current_parents = parents[current];
    for (int value : current_parents) {
      const int parent = value - 1;
      if (!ancestor_of_blocked[parent]) {
        ancestor_of_blocked[parent] = 1;
        ancestor_queue.push_back(parent);
      }
    }
  }

  std::vector<unsigned char> scheduled_up(n, 0);
  std::vector<unsigned char> scheduled_down(n, 0);
  std::deque<std::pair<int, int> > queue;
  queue.push_back(std::make_pair(exposure, 0));
  scheduled_up[exposure] = 1;

  while (!queue.empty()) {
    const int current = queue.front().first;
    const int direction = queue.front().second;
    queue.pop_front();

    if (current == outcome && !blocked[current]) {
      return false;
    }

    if (direction == 0 && !blocked[current]) {
      Rcpp::IntegerVector current_parents = parents[current];
      for (int value : current_parents) {
        const int node = value - 1;
        if (!scheduled_up[node]) {
          scheduled_up[node] = 1;
          queue.push_back(std::make_pair(node, 0));
        }
      }
      Rcpp::IntegerVector current_children = children[current];
      for (int value : current_children) {
        const int node = value - 1;
        if (!scheduled_down[node]) {
          scheduled_down[node] = 1;
          queue.push_back(std::make_pair(node, 1));
        }
      }
    } else if (direction == 1) {
      if (!blocked[current]) {
        Rcpp::IntegerVector current_children = children[current];
        for (int value : current_children) {
          const int node = value - 1;
          if (!scheduled_down[node]) {
            scheduled_down[node] = 1;
            queue.push_back(std::make_pair(node, 1));
          }
        }
      }
      if (ancestor_of_blocked[current]) {
        Rcpp::IntegerVector current_parents = parents[current];
        for (int value : current_parents) {
          const int node = value - 1;
          if (!scheduled_up[node]) {
            scheduled_up[node] = 1;
            queue.push_back(std::make_pair(node, 0));
          }
        }
      }
    }
  }
  return true;
}
