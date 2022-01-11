#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <Eigen/Dense>


namespace FEM {

using namespace std;

struct Source {
    Source() = default;

    Source(nlohmann::json j) 
    {
        if (j.at("magnitude") == nullptr ||
            j.at("type") == nullptr ||
            j.at("shape") == nullptr) {
            throw runtime_error("Invalid json object in Source ctor.");
        }

        magnitude = j.at("magnitude").get<double>();
        type = j.at("type").get<string>();
        shape = j.at("shape").get<string>();
    }

    double magnitude;
    string type;
    string shape;
};

nlohmann::json readInputData() 
{
    ifstream jfile("./testsData/data.json");
    cout << "Current path: " << filesystem::current_path().string() << endl;
    if (!jfile.is_open()) {
        throw runtime_error("Can not open.");
    }

    nlohmann::json j;
    jfile >> j;
    return j;
}

Eigen::MatrixXd initialiseContinousMatrix(const int maxnodes)
{
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(maxnodes,maxnodes);
    return res;
}

Eigen::MatrixXd initialiseDiscontinousMatrix(const int maxnodes)
{
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(2*maxnodes-2, 2*maxnodes-2);
    return res;
}

Eigen::VectorXd initialiseRightHandSideTerm(const int maxnodes)
{
    Eigen::VectorXd res = Eigen::VectorXd::Zero(maxnodes);
    return res;
}

Eigen::MatrixXd buildConnectionMatrix(const int maxnodes)
{
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(maxnodes, 2 * maxnodes - 2);
    for (int i = 0; i < maxnodes-1; i++) 
    {
        for (int j = 0; j < 2 * maxnodes - 3; j++) 
        {
            if (i == j / 2 + 1) { res(i, j) = 1.0; }
            else { res(i, j) = 0.0; }
        }
    }
    return res;
}


//double FetchVoltageValue(double x_value){
//    return;
//}

}