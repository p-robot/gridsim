#ifndef COMMON_F_H
#define COMMON_F_H

#define ping() std::cout << "*** Ping line " << __LINE__ << ", " << __FILE__ << " ***" << std::endl

#include <typeinfo>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <numeric>

#include "Point.h"

typedef std::chrono::high_resolution_clock::time_point TimePoint;
typedef double (*k_fun_ptr)(double); ///Typedef for a pointer to a kernel function.
typedef double (*d_fun_ptr)(Point, Point); ///Typedef for a pointer to a distance function.

namespace common
{
    TimePoint get_wall_time();
    unsigned long long int delta_duration(TimePoint t1);
    double sq_distance(Point p1, Point p2);
    double eucl_distance(Point p1, Point p2);
    double USDOS2_kernel(double d);
    double USDOS2_kernel_square(double d);
    double USDOS2_kernel_half_gamma(double d);
    double USDOS2_kernel_square_half_gamma(double d);
    double USDOS2_kernel_double_gamma(double d);
    double USDOS2_kernel_square_double_gamma(double d);
    double brand_kernel_a3(double d);
    double brand_kernel_a4(double d);
    double brand_kernel_a5(double d);
    double hayama_kernel(double d);
    double USDOS2_kernel_v2(double d);
    double gaussian_kernel(double d);
    double unif_rand();
    int rand_int(int lo, int hi);
    size_t draw_binom(size_t n, double prob);
    unsigned int get_n_lines(std::ifstream& f);
    void skipBOM(std::ifstream &in);
    std::vector<std::string> split(const std::string&,
                                   char, std::vector<std::string>&);
	std::vector<std::string> split(const std::string&, char);
	std::vector<int> split_stoi(const std::string &s, char delim);
	std::vector<double> split_stod(const std::string &s, char delim);
	double oneMinusExp(double x);

    template<typename T>
    double average(const std::vector<T>& elements)
    {
        double sum = 0.0;
        for(T element : elements)
        {
            sum += double(element);
        }
        return sum / double(elements.size());
    }

    template<typename T>
    double variance_from(const std::vector<T>& elements, double val)
    {
        std::vector<double> sq_diff;
        sq_diff.reserve(elements.size());
        for(T element : elements)
        {
            double diff = double(element) - val;
            sq_diff.push_back(diff * diff);
        }
        return average(sq_diff);
    }

    template<typename T>
    T stringToNum(const std::string& text)
    {
        std::istringstream ss(text);
        T result;
        if (! (ss>>result))
        {
            std::cout << "Failed to convert string " << text << " to "
                      << typeid(T).name() << std::endl;
            exit(EXIT_FAILURE);
        }
        return result;
    }

    ///> Chooses multiple random elements from a vector based on the Fisher-Yates shuffling algorithm.
    /// num_random selected random values are copied to output,
    /// selected values are swapped to the end of the vector so they're not selected again
    /// Used in Grid_manager to select binomial-success farms
    /// Used in Shipping_manager to select random premises in counties
    /// \param[in]	elements	Vector of values from which to choose
    ///	\param[in]	num_random	Number of values to choose
    /// \param[out] output		Vector of randomly chosen values
    template<typename T>
    void random_unique(std::vector<T> elements, size_t num_random, std::vector<T>& output)
        // elements not referenced (&) because we're rearranging it
    {
        std::vector<T> result;
        result.reserve(num_random);
        size_t endIndex = elements.size();
        // endIndex separates non-selected values (elements [0, endIndex-1]) from selected values
        for (size_t i = 1; i<= num_random; i++){
            // choose random number between 0 and 1
            double rUnif = unif_rand();
            // scale up to endIndex, so r is an index in [0, endIndex)
            if (rUnif == 1 ){rUnif=0.999;} // avoids assigning actual endIndex value (out of range)
            size_t r = (size_t)floor(rUnif*double(endIndex));
            // copy value to result
            result.emplace_back(elements[r]);
            // swap r and endIndex-1 (last element of non-selected values)
            std::swap(elements[r], elements[endIndex-1]);
            endIndex--;
        }
     result.swap(output);
    }
}
#endif // COMMON_F_H
