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

double calculateElementLength(int n, Coordinates coordinates)
{
    return (coordinates[1] - coordinates[0]) / (n - 1);
}

Eigen::VectorXd buildVoltageVector(int n, double voltage)
{
    Eigen::VectorXd res = Eigen::VectorXd::Zero(n);
    res(n-1) = voltage;
    return res;
}

Eigen::MatrixXd transposeVoltageToRightSide(Eigen::VectorXd voltageVector, ContinousMatrix continous)
{
    Eigen::VectorXd res = Eigen::VectorXd::Zero(voltageVector.rows());    
    for (int i = 0; i < voltageVector.rows(); i++) {
        res(i) = -1.0 * continous.coeff(i, voltageVector.rows() - 1) * voltageVector(voltageVector.rows() - 1);
    }
    
    return res;
}

ConnectionMatrix::ConnectionMatrix(int n)
{
    if (n < 2) {
        throw std::runtime_error("Invalid size for building ConnectionMatrix");
    }

    this->resize(n, 2 * n - 2);
    this->setZero();
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= 2 * n - 2; j++) {
            if (i == j / 2 + 1) {
                (*this)(i-1, j-1) = 1.0;
            }
        }
    }
}

DiscontinousMatrix::DiscontinousMatrix(InputData values)
{
    double elementLength = calculateElementLength(values.nodes, values.coordinates);

    if (values.nodes < 2) {
        throw std::runtime_error("Invalid size for building DiscontinousMatrix");
    }
    if (values.resistivity == 0) {
        throw std::runtime_error("Invalid resistivity value (cannot be zero)");
    }
    if (elementLength == 0) {
        throw std::runtime_error("Invalid element length value (cannot be zero)");
    }

    this->resize(2 * values.nodes - 2, 2 * values.nodes - 2);
    this->setZero();

    for (int i = 0; i < 2 * values.nodes - 3; i = i + 2) {
        int j = i + 1;
        (*this)(i, i) = (1.0 / (values.resistivity * elementLength)) + values.conductivity * elementLength / 3.0;
        (*this)(i, j) = -(1.0 / (values.resistivity * elementLength)) + values.conductivity * elementLength / 6.0;
        (*this)(j, i) = -(1.0 / (values.resistivity * elementLength)) + values.conductivity * elementLength / 6.0;
        (*this)(j, j) = (1.0 / (values.resistivity * elementLength)) + values.conductivity * elementLength / 3.0;
    }
}


ContinousMatrix::ContinousMatrix(int n, ConnectionMatrix connection, DiscontinousMatrix discontinous)
{

    if (n < 1) {
        throw std::runtime_error("Invalid size for building ContinousMatrix");
    }

    this->resize(n, n);
    this->setZero();

    for (int k = 0; k < n; k++) {
        for (int l = 0; l < n; l++) {
            double sum = 0;
            for (int j = 0; j < 2 * n - 2; j++) {
                for (int i = 0; i < 2 * n - 2; i++) {
                    sum = sum + connection.coeff(k, i) * discontinous.coeff(i, j) * connection.coeff(l, j);
                }
            }
            (*this)(k, l) = sum;
        }
    }
}


}
