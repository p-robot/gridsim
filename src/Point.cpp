#include "Point.h"

Point::Point()
{

}

Point::Point(double x, double y) :
    x(x), y(y)
{
    //ctor
}

Point::~Point()
{
    //dtor
}

void Point::set_x(double val)
{
    x = val;
}

void Point::set_y(double val)
{
    y = val;
}

bool Point::operator==(Point &other)
{
    return(this->get_x() == other.get_x() &&
           this->get_y() == other.get_y());
}

bool Point::operator!=(Point &other)
{
    return !(*this == other);
}
