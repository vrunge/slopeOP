// MIT License
// Copyright (c) 2019 Arnaud Liehrmann

#include "ListPoint.h"
#include "math.h"
using namespace std;
#include <iostream>
#include "stdlib.h"

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

ListPoint::ListPoint(){}

ListPoint::ListPoint(unsigned int size)
{
  resizeConst = size;
  lengthList = 0;
  head = (Point *) calloc(size+1, sizeof(Point));
  for(unsigned int i = 0; i < size; i++)
  {
    head[i].nxt = head+i+1;
  }

  head[size].nxt = NULL;
  tail = head + size;
  lastPosition = head -> nxt;
  currentPosition = head;
  headVec.push_back(head);
}

//####### destructor #######////####### destructor #######////####### destructor #######//
//####### destructor #######////####### destructor #######////####### destructor #######//


ListPoint::~ListPoint()
{
  for (unsigned int i = 0; i< headVec.size(); i++)
  {
    free(headVec[i]);
  }
}

//####### accessor #######////####### accessor #######////####### accessor #######//
//####### accessor #######////####### accessor #######////####### accessor #######//

unsigned int ListPoint::getState()
  {return(currentPosition -> state);}


unsigned int ListPoint::getTime()
  {return(currentPosition -> time);}



//####### move #######////####### move #######////####### move #######//
//####### move #######////####### move #######////####### move #######//


void ListPoint::move()
{
  currentPosition = currentPosition->nxt;
}

//####### addPoint #######////####### addPoint #######////####### addPoint #######//
//####### addPoint #######////####### addPoint #######////####### addPoint #######//

void ListPoint::addPoint(unsigned int s, unsigned int t)
{
  lastPosition->state = s;
  lastPosition->time = t;
  if (lastPosition->nxt == NULL){
    Point * newHead = (Point *) calloc(resizeConst, sizeof(Point));
    for (unsigned int i =0; i<resizeConst-1; i++)
    {
      newHead[i].nxt = newHead+i+1;
    }
    newHead[resizeConst-1].nxt = NULL;
    tail = newHead + resizeConst -1;
    lastPosition -> nxt = newHead;
    headVec.push_back(newHead);
  }
  lastPosition = lastPosition->nxt;
  lengthList++;
}


//####### deleteNxtPointAndMove #######////####### deleteNxtPointAndMove #######////####### deleteNxtPointAndMove #######//
//####### deleteNxtPointAndMove #######////####### deleteNxtPointAndMove #######////####### deleteNxtPointAndMove #######//


void ListPoint::deleteNxtPointAndMove()
{
  tail -> nxt = currentPosition -> nxt;
  currentPosition -> nxt = currentPosition->nxt->nxt;
  tail->nxt->nxt = NULL;
  tail = tail->nxt;
  lengthList--;
}


//####### initializeCurrentPosition #######////####### initializeCurrentPosition #######////####### initializeCurrentPosition #######//
//####### initializeCurrentPosition #######////####### initializeCurrentPosition #######////####### initializeCurrentPosition #######//

void ListPoint::initializeCurrentPosition()
{
  currentPosition = head;
}


//####### getLength #######////####### getLength #######////####### getLength #######//
//####### getLength #######////####### getLength #######////####### getLength #######//

unsigned int ListPoint::getLength()
{
  return lengthList;
}
