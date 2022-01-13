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


class ConnectionMatrix : public Eigen::MatrixXd {
public:
    ConnectionMatrix() = default;
    ConnectionMatrix(int nodes);
};

class DiscontinousMatrix : public Eigen::MatrixXd {
public:
    DiscontinousMatrix() = default;
    DiscontinousMatrix(InputData values);
};

class ContinousMatrix : public Eigen::MatrixXd {
public:
    ContinousMatrix() = default;
    ContinousMatrix(int n, ConnectionMatrix connection, DiscontinousMatrix discontinous);
};


nlohmann::json readInputData(std::string filename);

double calculateElementLength(int n, Coordinates coordinates);

Eigen::VectorXd buildPositionVector(int n, Coordinates coordinates);
Eigen::VectorXd buildVoltageVector(int n, double voltage);
Eigen::VectorXd calculateForwardElimination(ContinousMatrix continous, Eigen::MatrixXd rightHandSide);
Eigen::VectorXd calculateBackSubstitution(ContinousMatrix continous, Eigen::VectorXd voltageVector);

Eigen::MatrixXd calculateRightHandSide(Eigen::VectorXd voltageVector, ContinousMatrix continous);
Eigen::MatrixXd calculateProblemSolution(InputData values, Eigen::VectorXd voltageVector, Eigen::VectorXd positionVector);
Eigen::MatrixXd calculateCholeskyDecomposition(ContinousMatrix continous);


};