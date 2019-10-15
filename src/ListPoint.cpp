// MIT License
// Copyright (c) 2019 Vincent Runge

#include "ListPoint.h"
#include "math.h"
using namespace std;
#include <iostream>
#include "stdlib.h"

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

ListPoint::ListPoint(){}

ListPoint::ListPoint(unsigned int maxSize)
{
  listLength = maxSize;
  nbAdd = 0;

  arrayPoint = new Point[maxSize + 1]; ///with a first unused Point element
  currentPosition = arrayPoint;

  ///We lin all points
  for(unsigned int i = 0; i < maxSize; i++)
  {
    arrayPoint[i].time = 0;
    arrayPoint[i].state = 1000000;
    arrayPoint[i].nxt = &arrayPoint[i+1];
  }
}

//####### destructor #######////####### destructor #######////####### destructor #######//
//####### destructor #######////####### destructor #######////####### destructor #######//


ListPoint::~ListPoint()
{
  delete(arrayPoint);
  arrayPoint = NULL;
  currentPosition = NULL;
}

//####### accessor #######////####### accessor #######////####### accessor #######//
//####### accessor #######////####### accessor #######////####### accessor #######//

unsigned int ListPoint::getState()
  {return(currentPosition -> nxt -> state);}

unsigned int ListPoint::getTime()
  {return(currentPosition -> nxt -> time);}



//####### addPoint #######////####### addPoint #######////####### addPoint #######//
//####### addPoint #######////####### addPoint #######////####### addPoint #######//

void ListPoint::addPoint(unsigned int s, unsigned int t)
{
  arrayPoint[nbAdd + 1].state = s;
  arrayPoint[nbAdd + 1].time = t;
  nbAdd = nbAdd + 1;
}

//####### CanMove #######////####### CanMove #######////####### CanMove #######//
//####### CanMove #######////####### CanMove #######////####### CanMove #######//

bool ListPoint::canMove()
{
  if(currentPosition -> nxt -> nxt -> state != 1000000)
  {
    return(true);
  }
  return(false);
}

//####### move #######////####### move #######////####### move #######//
//####### move #######////####### move #######////####### move #######//

void ListPoint::move()
{
  currentPosition = currentPosition -> nxt;
}

//####### deletePoint #######////####### deletePoint #######////####### deletePoint #######//
//####### deletePoint #######////####### deletePoint #######////####### deletePoint #######//

void ListPoint::deletePoint()
{
  currentPosition -> nxt = currentPosition -> nxt -> nxt;
}

//####### initializeCurrentPosition #######////####### initializeCurrentPosition #######////####### initializeCurrentPosition #######//
//####### initializeCurrentPosition #######////####### initializeCurrentPosition #######////####### initializeCurrentPosition #######//

void ListPoint::initializeCurrentPosition()
{
  currentPosition = arrayPoint;
}
