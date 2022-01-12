#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <Eigen/Dense>


namespace FEM {

using namespace std;

typedef vector<double> Coordinates;

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

struct Values {
    Values() = default;

    Values(nlohmann::json j)
    {
        if (j.at("voltage") == nullptr ||
            j.at("resistivity") == nullptr ||
            j.at("conductivity") == nullptr ||
            j.at("coordinates") == nullptr ||
            j.at("nodes") == nullptr) {
            throw runtime_error("Invalid json object in Values ctor.");
        }
        voltage = j.at("voltage").get<double>();
        resistivity = j.at("resistivity").get<double>();
        conductivity = j.at("conductivity").get<double>();
        coordinates = j.at("coordinates").get<vector<double>>();
        nodes = j.at("nodes").get<int>();
    }
    double voltage;
    double resistivity;
    double conductivity;
    vector<double> coordinates;
    int nodes;
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

Eigen::MatrixXd buildConnectionMatrix(const int maxNodes)
{
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(maxNodes, 2 * maxNodes - 2);
    for (int i = 0; i < maxNodes-1; i++) {
        for (int j = 0; j < 2 * maxNodes - 3; j++) {
            if (i == j / 2 + 1) { 
                res(i, j) = 1.0; 
            }
        }
    }
    return res;
}


double calculateElementLength() 
{
    auto j = readInputData();
    vector<double> jCoords = j.at("coordinates"); auto jNodes = j.at("nodes");
    double res = (jCoords[1] - jCoords[0]) / (jNodes.get<int>() - 1);   
    return res;
}

Eigen::VectorXd buildVoltageVector() 
{
    auto j = readInputData(); Values values(j);
    Eigen::VectorXd res = Eigen::VectorXd::Zero(values.nodes);
    res(values.nodes) = values.voltage;
    return res;

}

Eigen::MatrixXd buildDiscontinousMatrix(const int maxNodes)
{
    auto j = readInputData();
    auto jResist = j.at("resistivity"); auto jConduct = j.at("conductivity"); auto jNodes = j.at("nodes");

    double elementLength = calculateElementLength();

    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(2 * maxNodes - 2, 2 * maxNodes - 2);

    for (int i = 0; i < 2 * jNodes.get<int>() - 4; i = i + 2) {
        int j = i + 1;
        res(i, i) = (1.0 / (jResist.get<double>() * elementLength)) + jConduct.get<double>() * elementLength / 3.0;
        res(i, j) = -(1.0 / (jResist.get<double>() * elementLength)) + jConduct.get<double>() * elementLength / 6.0; 
        res(j, i) = -(1.0 / (jResist.get<double>() * elementLength)) + jConduct.get<double>() * elementLength / 6.0; 
        res(j, j) = (1.0 / (jResist.get<double>() * elementLength)) + jConduct.get<double>() * elementLength / 3.0;
    }
    return res;
}

Eigen::MatrixXd buildContinousMatrix(const int maxNodes, Eigen::MatrixXd discontinousMatrix, Eigen::MatrixXd connectionMatrix)
{
    auto j = readInputData();
    auto jNodes = j.at("nodes");
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(maxNodes, maxNodes);

    for (int k = 0; k < jNodes.get<int>() - 1; k++) {
        for (int l = 0; l < jNodes.get<int>() - 1; l++) {
            double sum = 0;
            for (int j = 0; j < 2 * jNodes.get<int>() - 3; j++) {
                for (int i = 0; i < 2 * jNodes.get<int>() - 3; i++) {
                    sum = sum + connectionMatrix.coeff(k, i) * discontinousMatrix.coeff(i, j) * connectionMatrix.coeff(l, j);
                }
            }
            res(k, l) = sum;
        }
    }
    return res;
}   

Eigen::VectorXd buildRightHandSideTerm(const int maxNodes, Eigen::MatrixXd continousMatrix, Eigen::VectorXd voltageVector) 
{
    
    auto j = readInputData();
    auto jNodes = j.at("nodes");
    Eigen::VectorXd res = Eigen::VectorXd::Zero(maxNodes);

    for (int i = 0; i < jNodes.get<int>() - 2; i++) {
        res(i) = -1.0 * continousMatrix.coeff(i, jNodes.get<int>()-1) * voltageVector.coeff(jNodes.get<int>()-1);
    }

    return res;
}


//double FetchVoltageValue(double x_value){
//    return;
//}

}