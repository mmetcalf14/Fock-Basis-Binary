#include "Node.h"

template<class T>
Node<T>::Node(){}

//  This constructor is just to set next pointer of a node and the data contained.
template<class T>
Node<T>::Node(const T& item, Node<T>* p)
{
  this->data = item;
  if ( p != NULL ){
    this->NumLinks++;
    this->Neighbor.push_back(p);
  }
}

template<class T>
void Node<T>::LinkTo(Node<T>* p)
{
  this->NumLinks++;
  this->Neighbor.push_back(p);
}

// //  This methods inserts a node just after the node that the method belongs to
// //  TODO: Consider a better implementation
// template<class T>
// void Node<T>::InsertAfter(Node<T> *p)
// {
//   // not to lose the rest of the list, we ought to link the rest of the
//   // list to the Node<T>* p first
//   p->next = this->next;
//   // now we should link the previous Node to Node<T> *p , i.e the Node that we are
//   //inserting after,
//   this->next = p;
// }
//
// // Deletes the node from the list and returns the deleted node
// template<class T>
// Node<T>* Node<T>::DeleteAfter()
// {
//   // store the next Node in a temporary Node
//   Node<T>* tempNode = next;
//   // check if there is a next node
//   if(next != NULL)
//     next = next->next;
//
//   return tempNode;
// }
//
// template<class T>
// Node<T>* GetNode(const T& item, Node<T>* nextptr = NULL)
// {
//   Node<T>* newnode; // Local ptr for new node
//   newnode = new Node<T>(item,nextptr);
//   if ( newnode == NULL)
//   {
//     cerr << "Memory allocation failed." << endl;
//     exit(1);
//   }
//   return newnode;
// }
