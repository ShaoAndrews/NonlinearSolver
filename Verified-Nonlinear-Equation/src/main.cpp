#include <kv/allsol.hpp>
#include <getopt.h>
#include <string>
#include <fstream>
#include <Eigen/Dense>
#include <iostream>
using namespace std;
static Eigen::MatrixXd f1_mat;
static Eigen::MatrixXd f2_mat;
namespace ub = boost::numeric::ublas;
struct Func
{
    template <class T>
    ub::vector<T> operator()(const ub::vector<T> &x)
    {
        ub::vector<T> y(2);
        // y(0) = ((((((2)*(pow(x(0),3)))+(((2)*(x(0)))*(x(1))))-((22)*(x(0))))+(pow(x(1),2)))+(10))+(pi);
        // y(1) = ((((pow(x(0),2))+(((2)*(x(0)))*(sin(x(1)))))+((2)*(pow(x(1),3))))-((14)*(x(1))))+c0;

        //初始化y(0)
        y(0) = 0;
        y(1) = 0;
        for (int i = 0; i < f1_mat.rows(); i++)
        {
            for (int j = 0; j < f1_mat.cols(); j++)
            {
                y(0) += pow(x(0), i) * pow(x(1), j) * f1_mat(i, j);
            }
        }
        for (int i = 0; i < f2_mat.rows(); i++)
        {
            for (int j = 0; j < f2_mat.cols(); j++)
            {
                y(1) += pow(x(0), i) * pow(x(1), j) * f2_mat(i, j);
            }
        }
        return y;
    }
};

using namespace Eigen;
template <typename T>
T load_csv(const std::string &path)
{
    std::ifstream in;
    in.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(in, line))
    {
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
        {
            double val = std::stod(cell);
            values.push_back(val);
        }
        ++rows;
    }
    return Eigen::Map<const Eigen::Matrix<
        typename T::Scalar,
        T::RowsAtCompileTime,
        T::ColsAtCompileTime,
        Eigen::RowMajor>>(values.data(), rows, values.size() / rows);
}

void usage()
{
    std::cout << "Options:" << endl;
    std::cout << " -a,  --f1=f1.csv           A csv file of polynomial about f1" << endl;
    std::cout << " -b,  --f2=f2.csv           A csv file of polynomial about f2" << endl;
    std::cout << " -o,  --out=OUT             The path of output file" << endl;
    std::cout << "                            default value is './curves@f1@f2.csv'" << endl;
    std::cout << "                            The fist line records the time, and " << endl;
    std::cout << "                            the second line records the number of roots" << endl;
    std::cout << " -t,  --tolerance           Set default precision 1e-6" << endl;
    std::cout << " -h,  --help(no value)      Print the message and exit" << endl;
}

int main(int argc, char *argv[])
{
    string f1_path;
    string f2_path;
    double precision = 1e-6;
    string outpath = "./result.json";
    static struct option long_options[] = {
        {"f1", required_argument, 0, 'a'},
        {"f2", required_argument, 0, 'b'},
        {"out", optional_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {"tolerance", optional_argument, 0, 't'},
        {0, 0, 0, 0}};
    int getopt_ret, option_index;
    bool status = true;
    while (1)
    {
        getopt_ret = getopt_long(argc, argv, "a:b:o::ht::", long_options, &option_index);
        if (getopt_ret == -1)
            break;
        switch (getopt_ret)
        {
        case 0:
            status = false;
            break;
        case 'a':
            f1_path = optarg;
            break;
        case 'b':
            f2_path = optarg;
            break;
        case 'o':
            outpath = optarg;
            break;
        case 'h':
            usage();
            status = false;
            break;
        case 't':
            precision = atof(optarg);
            break;
        case '?':
            status = false;
            break;
        default:
            status = false;
            break;
        }
    }
    if (status)
    {
        f1_mat = load_csv<Eigen::MatrixXd>(f1_path);
        f2_mat = load_csv<Eigen::MatrixXd>(f2_path);
    }

    ub::vector<kv::interval<double>> x;
    std::cout.precision(17);
    x.resize(2);
    x(0) = kv::interval<double>(-1, 1);
    x(1) = kv::interval<double>(0, 1);
    allsol(Func(), x);
}
