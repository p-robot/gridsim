#ifndef POINT_H
#define POINT_H

#include <sstream>

/**
A point on the x-y plane.
*/

class Point
{
    public:
        /**
        Create using x and y coordinates.
        */
        Point();
        Point(double x, double y);
        ~Point();
        double get_x(); //Inlined
        double get_y(); //Inlined
        void set_x(double val);
        void set_y(double val);
        bool operator==(Point &other);
        bool operator!=(Point &other);

    private:
        double x;
        double y;
};

inline std::ostream& operator<<(std::ostream& stream, Point o)
{
    stream << "[" << o.get_x() << "; " << o.get_y() << "]";
    return stream;
}

inline double Point::get_x()
{
    return x;
}

inline double Point::get_y()
{
    return y;
}

#endif // POINT_H
