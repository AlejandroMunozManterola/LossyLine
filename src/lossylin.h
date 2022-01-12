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

//Eigen::MatrixXd buildConnectionMatrix(const int maxNodes)
//{
//    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(maxNodes, 2 * maxNodes - 2);
//    for (int i = 0; i < maxNodes - 1; i++) {
//        for (int j = 0; j < 2 * maxNodes - 3; j++) {
//            if (i == j / 2 + 1) {
//                res(i, j) = 1.0;
//            }
//        }
//    }
//    return res;
//}

//
//Eigen::MatrixXd buildDiscontinousMatrix(int maxNodes)
//{
//    auto j = readInputData(); InputData values(j);
//    double elementLength = calculateElementLength();
//
//    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(2 * maxNodes - 2, 2 * maxNodes - 2);
//
//    for (int i = 0; i < 2 * values.nodes - 4; i = i + 2) {
//        int j = i + 1;
//        res(i, i) = (1.0 / (values.resistivity * elementLength)) + values.conductivity * elementLength / 3.0;
//        res(i, j) = -(1.0 / (values.resistivity * elementLength)) + values.conductivity * elementLength / 6.0; 
//        res(j, i) = -(1.0 / (values.resistivity * elementLength)) + values.conductivity * elementLength / 6.0; 
//        res(j, j) = (1.0 / (values.resistivity * elementLength)) + values.conductivity * elementLength / 3.0;
//    }
//    return res;
//}
//
//
//
//Eigen::MatrixXd buildContinousMatrix(const int maxNodes, Eigen::MatrixXd discontinousMatrix, Eigen::MatrixXd connection)
//{
//    auto j = readInputData(); InputData values(j);
//    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(maxNodes, maxNodes);
//
//    for (int k = 0; k < values.nodes - 1; k++) {
//        for (int l = 0; l < values.nodes - 1; l++) {
//            double sum = 0;
//            for (int j = 0; j < 2 * values.nodes - 3; j++) {
//                for (int i = 0; i < 2 * values.nodes - 3; i++) {
//                    sum = sum + connection.coeff(k, i) * discontinousMatrix.coeff(i, j) * connection.coeff(l, j);
//                }
//            }
//            res(k, l) = sum;
//        }
//    }
//    return res;
//}   

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


//double FetchVoltageValue(double x_value){
//    return;
//}
