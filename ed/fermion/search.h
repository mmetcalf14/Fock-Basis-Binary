#ifndef __SEARCH_H__
#define __SEARCH_H__
#include <algorithm>

template<class RandomIt, class T>
inline RandomIt binary_locate(RandomIt first, RandomIt last, const T& val) {
  if(val == *first) return first;
  auto d = std::distance(first, last);
  if(d==1) return first;
  auto center = (first + (d/2));
  if(val < *center) return binary_locate(first, center, val);
  return binary_locate(center, last, val);
}

#endif//__SEARCH_H__
