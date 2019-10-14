// MIT License
// Copyright (c) 2019 Vincent Runge

#ifndef LISTPOINT_H
#define LISTPOINT_H

#include "Point.h"
#include <vector>
#include <math.h>

class ListPoint
{
public:
  ListPoint();
  ListPoint(unsigned int maxSize);
  ~ListPoint();

  unsigned int getState();
  unsigned int getTime();

  void addPoint(unsigned int s, unsigned int t);
  bool move();
  void deletePoint();
  void initializeCurrentPosition();

private:
  Point* arrayPoint;
  Point* currentPosition;
  unsigned int nbAdd;

  unsigned int listLength;

};

#endif // LISTPOINT_H
