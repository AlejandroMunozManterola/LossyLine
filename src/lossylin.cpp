#include "lossylin.h"

namespace FEM {

InputData::InputData(nlohmann::json j)
{
    if (j.at("voltage") == nullptr ||
        j.at("resistivity") == nullptr ||
        j.at("conductivity") == nullptr ||
        j.at("coordinates") == nullptr ||
        j.at("nodes") == nullptr) {
        throw std::runtime_error("Invalid json object in InputData ctor.");
    }

    voltage = j.at("voltage").get<double>();
    resistivity = j.at("resistivity").get<double>();
    conductivity = j.at("conductivity").get<double>();
    coordinates = j.at("coordinates").get<Coordinates>();
    nodes = j.at("nodes").get<int>();

}

nlohmann::json readInputData(std::string filename)
{
    std::ifstream jfile(filename);

    if (!jfile.is_open()) {
        throw std::runtime_error("Can not open.");
    }

    nlohmann::json j;
    jfile >> j;
    return j;
}

double calculateElementLength(InputData values)
{
    return (values.coordinates[1] - values.coordinates[0]) / (values.nodes - 1);
}

Eigen::VectorXd buildVoltageVector(InputData values)
{
    Eigen::VectorXd res = Eigen::VectorXd::Zero(values.nodes);
    res(values.nodes) = values.voltage;
    return res;
}


ConnectionMatrix::ConnectionMatrix(int maxNodes)
{
    if (maxNodes < 1) {
        throw std::runtime_error("Invalid size for building ConnectionMatrix");
    }

    this->resize(maxNodes, 2 * maxNodes - 2);
    this->setZero();
    for (int i = 0; i < maxNodes - 1; i++) {
        for (int j = 0; j < 2 * maxNodes - 3; j++) {
            if (i == j / 2 + 1) {
                (*this)(i, j) = 1.0;
            }
        }
    }
}


DiscontinousMatrix::DiscontinousMatrix(int maxNodes, InputData values)
{
    double elementLength = calculateElementLength(values);

    if (maxNodes < 1) {
        throw std::runtime_error("Invalid size for building DiscontinousMatrix");
    }
    if (values.resistivity == 0) {
        throw std::runtime_error("Invalid resistivity value (cannot be zero)");
    }
    if (elementLength == 0) {
        throw std::runtime_error("Invalid element length value (cannot be zero)");
    }

    this->resize(2 * maxNodes - 2, 2 * maxNodes - 2);
    this->setZero();

    for (int i = 0; i < 2 * values.nodes - 4; i = i + 2) {
        int j = i + 1;
        (*this)(i, i) = (1.0 / (values.resistivity * elementLength)) + values.conductivity * elementLength / 3.0;
        (*this)(i, j) = -(1.0 / (values.resistivity * elementLength)) + values.conductivity * elementLength / 6.0;
        (*this)(j, i) = -(1.0 / (values.resistivity * elementLength)) + values.conductivity * elementLength / 6.0;
        (*this)(j, j) = (1.0 / (values.resistivity * elementLength)) + values.conductivity * elementLength / 3.0;
    }
}

ContinousMatrix::ContinousMatrix(int maxNodes, InputData values)
{

    if (maxNodes < 1) {
        throw std::runtime_error("Invalid size for building ContinousMatrix");
    }

    ConnectionMatrix connection;
    DiscontinousMatrix discontinous;

    this->resize(maxNodes, maxNodes);
    this->setZero();

    for (int k = 0; k < values.nodes - 1; k++) {
        for (int l = 0; l < values.nodes - 1; l++) {
            double sum = 0;
            for (int j = 0; j < 2 * values.nodes - 3; j++) {
                for (int i = 0; i < 2 * values.nodes - 3; i++) {
                    sum = sum + connection.coeff(k, i) * discontinous.coeff(i, j) * connection.coeff(l, j);
                }
            }
            (*this)(k, l) = sum;
        }
    }
}


}
