#ifndef __SITE_H__
#define __SITE_H__

class Site {
private:
  static NumOfSites;
protected:
  std::string Name;
  size_t NnumOrbits;
public:
  Site( const std::string &name, const size_t &orbits );
}

#endif//__SITE_H__
