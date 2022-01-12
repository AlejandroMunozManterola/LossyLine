#pragma once

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>

#include <nlohmann/json.hpp>
#include <Eigen/Dense>

namespace FEM {

typedef std::vector<double> Coordinates;

struct InputData {
    InputData() = default;
    InputData(nlohmann::json j);

    double voltage = 0.0;
    double resistivity = 0.0;
    double conductivity = 0.0;
    Coordinates coordinates{};
    int nodes = 0;
};

nlohmann::json readInputData(std::string filename);
double calculateElementLength(InputData values);
Eigen::VectorXd buildVoltageVector(InputData values);

class ConnectionMatrix : public Eigen::MatrixXd {
public:
    ConnectionMatrix() = default;
    ConnectionMatrix(int maxNodes);
};

class DiscontinousMatrix : public Eigen::MatrixXd {
public:
    DiscontinousMatrix() = default;
    DiscontinousMatrix(int maxNodes, InputData values);
};

class ContinousMatrix : public Eigen::MatrixXd {
public:
    ContinousMatrix() = default;
    ContinousMatrix(int maxNodes, InputData values);
};

};


//Eigen::VectorXd buildRightHandSideTerm(const int maxNodes, Eigen::MatrixXd continousMatrix, Eigen::VectorXd voltageVector) 
//{
//    
//    auto j = readInputData(); InputData values(j);
//    Eigen::VectorXd res = Eigen::VectorXd::Zero(maxNodes);
//
//    for (int i = 0; i < values.nodes - 2; i++) {
//        res(i) = -1.0 * continousMatrix.coeff(i, values.nodes-1) * voltageVector.coeff(values.nodes-1);
//    }
//
//    return res;
//}
