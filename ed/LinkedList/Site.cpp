#include "LinkedList/Site.hpp"

template<typename Tnum, class T>
int Site<Tnum, T>::NumSites = 0;

template<typename Tnum, class T>
Site<Tnum, T>::Site()
{
  this->data = NumSites;
  NumSites++;
}

template<typename Tnum, class T>
Site<Tnum, T>::~Site(){}

//  This constructor is just to set next pointer of a Site and the data contained.
template<typename Tnum, class T>
Site<Tnum, T>::Site(const T& item, Site<Tnum, T>* p)
{
  NumSites++;
  this->data = item;
  if ( p != NULL ){
    this->NumLinks++;
    this->Neighbor.push_back(p);
  }
}

template<typename Tnum, class T>
void Site<Tnum, T>::LinkTo(Site<Tnum, T>* p, Tnum J)
{
  this->NumLinks++;
  this->Neighbor.push_back(p);
  this->Jval.push_back(J);
}

// //  This methods inserts a Site just after the Site that the method belongs to
// //  TODO: Consider a better implementation
// template<typename Tnum, class T>
// void Site<Tnum, T>::InsertAfter(Site<Tnum, T> *p)
// {
//   // not to lose the rest of the list, we ought to link the rest of the
//   // list to the Site<Tnum, T>* p first
//   p->next = this->next;
//   // now we should link the previous Site to Site<Tnum, T> *p , i.e the Site that we are
//   //inserting after,
//   this->next = p;
// }
//
// // Deletes the Site from the list and returns the deleted Site
// template<typename Tnum, class T>
// Site<Tnum, T>* Site<Tnum, T>::DeleteAfter()
// {
//   // store the next Site in a temporary Site
//   Site<Tnum, T>* tempSite = next;
//   // check if there is a next Site
//   if(next != NULL)
//     next = next->next;
//
//   return tempSite;
// }
//
// template<typename Tnum, class T>
// Site<Tnum, T>* GetSite(const T& item, Site<Tnum, T>* nextptr = NULL)
// {
//   Site<Tnum, T>* newSite; // Local ptr for new Site
//   newSite = new Site<Tnum, T>(item,nextptr);
//   if ( newSite == NULL)
//   {
//     cerr << "Memory allocation failed." << endl;
//     exit(1);
//   }
//   return newSite;
// }

template class Site<double, int>;
template class Site<double, char>;
