#ifndef __SYSTEM_HPP__
#define __SYSTEM_HPP__
#include <vector>

class System {
private:
  int TotalSites = 0;
  std::vector< Node<int>* > SiteList;
public:
  System ();
  virtual ~System ();
  inline void AddSite()
  {
    Node<int>* p = new Node<int>(TotalSites);
    SiteList.push_back(p);
    TotalSites++;
  };
  inline void AddLink( const int lhs, const int rhs )
  {
    SiteList.at(lhs)->LinkTo( SiteList.at(rhs) );
    SiteList.at(rhs)->LinkTo( SiteList.at(lhs) );
  }
};

#endif//__SYSTEM_HPP__
