// MIT License
// Copyright (c) 2019 Vincent Runge

#ifndef LISTPOINT_H
#define LISTPOINT_H

#include "Point.h"
#include <math.h>

class ListPoint ///DANGER: THERE IS "NO" EMPTY LIST. ALWAYS ADD SOME ELEMENTS
{
public:
  ListPoint(); //create an empty list
  unsigned int getState();
  unsigned int getTime();

  void initializeCurrentPosition(); //currentPosition = firstPoint
  void addPoint(Point* P); //at the beginning after (0,0)
  void deleteNxtPoint();
  bool move(); //DANGER we can move only in a non-empty list

  void show();

private:
  Point* firstPoint;
  Point* currentPosition;
};

#endif // LISTPOINT_H
