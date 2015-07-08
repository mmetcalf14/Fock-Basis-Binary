#ifndef __NODE_H__
#define __NODE_H__

template<class T>
class Node
{
public:
  Node();
  Node(const T& item, Node<T>* ptrnext = NULL);
  T data;
  LinkTo(Node<T>* p);
  inline int NumNeighbors()const{return NumLinks;};
private:
  int NumLinks = 0;
  std::vector< Node<T>* > Neighbor;
};

#endif//__NODE_H__
