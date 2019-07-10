// MIT License
// Copyright (c) 2019 Vincent Runge

#ifndef POINT_H
#define POINT_H

#include <math.h>

class Point
{
  public:
    Point(unsigned int s, unsigned int t);
    unsigned int state;
    unsigned int time;
    Point* nxt;
};

#endif // POINT_H
