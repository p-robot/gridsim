#define _USE_MATH_DEFINES

#include "common_functions.h"
#include <math.h>
#include <random>
#include <chrono>
#include <functional>

TimePoint common::get_wall_time()
{
    return std::chrono::high_resolution_clock::now();
}

unsigned long long int common::delta_duration(TimePoint t1)
{
    auto d = std::chrono::duration_cast<std::chrono::milliseconds>(common::get_wall_time() - t1).count();
    return (unsigned long long int)(d);
}

//Euclidean distance squared
double common::sq_distance(Point p1, Point p2)
{
    double delta_x = p1.get_x() - p2.get_x();
    double delta_y = p1.get_y() - p2.get_y();
    return (delta_x*delta_x) + (delta_y*delta_y);
}

double common::eucl_distance(Point p1, Point p2)
{
    return std::sqrt(common::sq_distance(p1, p2));
}

//This is the Buhnerkempe kernel.
double common::USDOS2_kernel(double d)
{
    //Original
    double k1 = 0.08912676;
    double k2 = 1600;
    double k3 = 4.6;
    return k1 / (1 + std::pow((d/k2),k3));
}

double common::USDOS2_kernel_square(double d)
{
    //Original
    //Same as USDOS2_kernel but accepts distance^2 instead of distance.
    //This means that sqrt operation can be avoided when calc. distance,
    //which potentially could be a little bit faster. This was not
    //used in the study.
    double k1 = 0.08912676;
    double k2 = 1600.0;
    double k3 = 4.6;
    return k1 / (1 + std::pow(d, (k3*0.5)) / std::pow(k2, k3));

}

double common::USDOS2_kernel_half_gamma(double d)
{
    double k1 = 0.01813338;
    double k2 = 1600;
    double k3 = 2.3;
    return k1 / (1 + std::pow((d/k2),k3));
}

double common::USDOS2_kernel_square_half_gamma(double d)
{
    //Same as kernel() but accepts distance^2 instead of distance.
    //This means that sqrt operation can be avoided when calc. distance,
    //which potentially could be a little bit faster. This was not
    //used in the study.
    double k1 = 0.01813338;
    double k2 = 1600.0;
    double k3 = 4.6;
    return k1 / (1 + std::pow(d, (k3*0.5)) / std::pow(k2, k3));

}

double common::USDOS2_kernel_double_gamma(double d)
{
    double k1 = 0.11489682;
    double k2 = 1600;
    double k3 = 9.2;
    return k1 / (1 + std::pow((d/k2),k3));
}

double common::USDOS2_kernel_square_double_gamma(double d)
{
    //Same as USDOS2_kernel but accepts distance^2 instead of distance.
    //This means that sqrt operation can be avoided when calc. distance,
    //which potentially could be a little bit faster. This was not
    //used in the study.
    double k1 = 0.11489682;
    double k2 = 1600.0;
    double k3 = 9.2;
    return k1 / (1 + std::pow(d, (k3*0.5)) / std::pow(k2, k3));

}

double common::gaussian_kernel(double d)
{
    //Bivariate normal but mu_x=mu_y=0 and sigma_x = sigma_y and rho=0,
    //d is x^2 + y^2 which is the squared eucl distance.
    double L_1_norm = 0.2*(1000.0*1000.0);
    double length_scale = 3.0*1000.0;
    double fact_1 = 1.0 / (2*M_PI*length_scale*length_scale);
    return L_1_norm * fact_1 * std::exp(-(d / (2*length_scale*length_scale)));
}

double common::brand_kernel_a3(double d)
{
    d = d / 1000.0; //Transform to km
    double beta = 0.115;
    double L3 = 10.0;
    double N3 = 1.0 / (2*M_PI);
    return beta * N3 * L3 * std::pow(L3*L3 + d * d, -1.5);
}

double common::brand_kernel_a4(double d)
{
    d = d / 1000.0; //Transform to km
    double beta = 0.115;
    double L4 = 20.0 / (0.5*M_PI);
    double N4 = L4 / M_PI;
    return beta * N4 * L4 * std::pow(L4*L4 + d * d, -2.0);
}

double common::brand_kernel_a5(double d)
{
    d = d / 1000.0; //Transform to km
    double beta = 0.115;
    double L5 = 20.0;
    double N5 = (3.0 * L5 * L5) / (2 * M_PI);
    return beta * N5 * L5 * std::pow(L5*L5 + d*d, -2.5);
}

double common::hayama_kernel(double d)
{
    d = d / 1000.0; //Transform to km
    double Hayama_r0 = 0.58;
    double Hayama_h0 = 0.00074;
    double Hayama_alpha = 2.47;
    return Hayama_h0 * std::pow(1 + (d/Hayama_r0), -Hayama_alpha);
}

double common::USDOS2_kernel_v2(double d)
{
    double USDOS_k1 = 1.293833e-08;
    double USDOS_k2 = 2116.798;
    double USDOS_k3 = 2.38;
    return USDOS_k1 / (1 + std::pow(d/USDOS_k2, USDOS_k3));
}

//Returns number of lines in a file.
unsigned int common::get_n_lines(std::ifstream& f)
{
    unsigned int n_lines = 0;
    std::string tempstring;
    f.seekg(0, f.beg);
    while(std::getline(f, tempstring))
        n_lines++;
    f.clear();
    f.seekg(0, f.beg);
    return n_lines;
}

double common::unif_rand()
{
	static unsigned long long int seed = std::chrono::system_clock::now().time_since_epoch().count();
	static std::mt19937 generator(seed); //Mersenne Twister pseudo-random number generator. Generally considered research-grade.
    static std::uniform_real_distribution<double> unif_dist(0.0, 1.0);
	return unif_dist(generator);
}

int common::rand_int(int lo, int hi)
{
	std::uniform_int_distribution<int> unif_dist(lo, hi);
	static unsigned long long int seed = std::chrono::system_clock::now().time_since_epoch().count();
	static std::mt19937 generator(seed); //Mersenne Twister pseudo-random number generator. Generally considered research-grade.
	return unif_dist(generator);
}

size_t common::draw_binom(size_t n, double prob)
// draw from a binomial distribution based on N farms and prob (calc with focalInf & gridKern)
{
	std::binomial_distribution<size_t> binom_dist(n, prob);
	static unsigned long long int seed = std::chrono::system_clock::now().time_since_epoch().count();
	static std::mt19937 generator(seed); //Mersenne Twister pseudo-random number generator. Generally considered research-grade.
	return binom_dist(generator);
}

// Skips the Byte Order Mark (BOM) that defines UTF-8 in some text files.
// Credit to user 'Contango' at stackoverflow.com
void common::skipBOM(std::ifstream &in)
{
    char test[3] = {0};
    in.read(test, 3);
    if ((unsigned char)test[0] == 0xEF &&
        (unsigned char)test[1] == 0xBB &&
        (unsigned char)test[2] == 0xBF)
    {
        return;
    }
    in.seekg(0);
}

//Split a string into a vector based on arbitrary delimiter.
std::vector<std::string> common::split(const std::string &s,
                                       char delim, std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
	{
		if (!item.empty())
		{
			elems.emplace_back(item);
		}
    }
    return elems;
}

std::vector<std::string> common::split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

std::vector<int> common::split_stoi(const std::string &s, char delim)
{
    std::vector<std::string> elems = common::split(s, delim);
    std::vector<int> elems_int;
    for(std::string this_s : elems)
        elems_int.push_back(std::stoi(this_s));
    return elems_int;
}

std::vector<double> common::split_stod(const std::string &s, char delim)
{
    std::vector<std::string> elems = common::split(s, delim);
    std::vector<double> elems_int;
    for(std::string this_s : elems)
        elems_int.push_back(std::stod(this_s));
    return elems_int;
}

// Based on algorithm described at http://www.johndcook.com/blog/cpp_expm1/
double common::oneMinusExp(double x)
{
	if (x == 0.0)
    {
		return 0.0;
	}
	else if (std::abs(x) < 1e-5)
    {
		return -(x + 0.5*x*x);
	}
	else
    {
		return -(exp(x) - 1.0);
	}
}
