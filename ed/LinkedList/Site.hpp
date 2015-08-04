#ifndef __Site_H__
#define __Site_H__
#include <vector>

template<typename Tnum = double, class T = int>
class Site
{
public:
  Site();
  Site(const T& item, Site<Tnum, T>* ptrnext = NULL);
  virtual ~Site();
  T data;
  void LinkTo(Site<Tnum, T>* p, Tnum J);
  inline int TotalSites()const{return NumSites;};
  inline int NumNeighbors()const{return NumLinks;};
  inline std::vector< Site<Tnum, T>* > getNeighbors()const{return Neighbor;};
  inline bool VerifySite()const{return Neighbor.size() == Jval.size();};
private:
  static int NumSites;
  int NumLinks = 0;
  std::vector< Site<Tnum, T>* > Neighbor;
  std::vector< Tnum > Jval;
};

#endif//__Site_H__
