// MIT License
// Copyright (c) 2019 Vincent Runge

#include "ListPoint.h"
#include<iostream>

#include "math.h"

//####### constructor #######////####### constructor #######////####### constructor #######//
//####### constructor #######////####### constructor #######////####### constructor #######//

ListPoint::ListPoint()
{
  firstPoint = new Point(0,0);
  currentPosition = firstPoint;
  //the first Point = an empty point : we consider point->nxt in functions
}


//####### accessor #######////####### accessor #######////####### accessor #######//
//####### accessor #######////####### accessor #######////####### accessor #######//

unsigned int ListPoint::getState()
  {return(currentPosition -> nxt -> state);}


unsigned int ListPoint::getTime()
  {return(currentPosition -> nxt -> time);}

//####### initializeCurrentPosition #######////####### initializeCurrentPosition #######////####### initializeCurrentPosition #######//
//####### initializeCurrentPosition #######////####### initializeCurrentPosition #######////####### initializeCurrentPosition #######//

void ListPoint::initializeCurrentPosition()
  {currentPosition = firstPoint;}

//####### addPoint #######////####### addPoint #######////####### addPoint #######//
//####### addPoint #######////####### addPoint #######////####### addPoint #######//

void ListPoint::addPoint(Point* P)
{
  P -> nxt = firstPoint -> nxt;
  firstPoint -> nxt = P;
}


//####### deleteNxtPoint #######////####### deleteNxtPoint #######////####### deleteNxtPoint #######//
//####### deleteNxtPoint #######////####### deleteNxtPoint #######////####### deleteNxtPoint #######//

void ListPoint::deleteNxtPoint()
{
  Point* pointToDelete = currentPosition -> nxt;
  currentPosition -> nxt = currentPosition -> nxt -> nxt;
  delete(pointToDelete);
}


//####### move #######////####### move #######////####### move #######//
//####### move #######////####### move #######////####### move #######//

bool ListPoint::move()
{
  bool res = true;
  if(currentPosition -> nxt -> nxt != NULL){currentPosition = currentPosition -> nxt;}
    else{res = false;}

  return(res);
}


//####### show #######////####### show #######////####### show #######//
//####### show #######////####### show #######////####### show #######//

void ListPoint::show()
{
  //std::cout << "State " << getState() << "  Time " << getTime() << std::endl;
  // while(move() == true){std::cout << "State " << getState() << "  Time " << getTime() << std::endl;}
}




