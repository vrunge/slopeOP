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
  ListPoint(unsigned int size);
  ~ListPoint();

  unsigned int getState();
  unsigned int getTime();

  void move();
  void addPoint(unsigned int a, unsigned int b);
  void deleteNxtPointAndMove();
  void initializeCurrentPosition();
  unsigned int getLength();

private:
  Point * head;
  Point * currentPosition;
  Point * lastPosition;
  Point * tail;
  unsigned int lengthList;
  unsigned int resizeConst;
  std::vector<Point *> headVec;
};

#endif // LISTPOINT_H
