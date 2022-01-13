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

Eigen::VectorXd buildPositionVector(int n, Coordinates coordinates) 
{
    Eigen::VectorXd res = Eigen::VectorXd::Zero(n);
    return res.setLinSpaced(n, coordinates[0], coordinates[1]);
}

Eigen::VectorXd buildVoltageVector(int n, double voltage)
{
    Eigen::VectorXd res = Eigen::VectorXd::Zero(n);
    res(n-1) = voltage;
    return res;
}

Eigen::MatrixXd calculateRightHandSide(Eigen::VectorXd voltageVector, ContinousMatrix continous)
{
    Eigen::VectorXd res = Eigen::VectorXd::Zero(voltageVector.rows());    
    for (int i = 0; i < voltageVector.rows(); i++) {
        res(i) = -1.0 * continous.coeff(i, voltageVector.rows() - 1) * voltageVector(voltageVector.rows() - 1);
    }
    return res;
}

Eigen::MatrixXd calculateCholeskyDecomposition(ContinousMatrix continous) 
{
    return continous.llt().matrixL();
}

Eigen::VectorXd calculateForwardElimination(ContinousMatrix continous, Eigen::MatrixXd rightHandSide) 
{
    Eigen::VectorXd res = Eigen::VectorXd::Zero(continous.rows());
    res(0) = rightHandSide.coeff(0) / continous.coeff(0, 0);
    for (int i = 1; i < res.rows(); i++) {
        res(i) = rightHandSide.coeff(i);
        for (int j = 0; j < i - 1; j++) {
            res(i) = res.coeff(i) - continous.coeff(i, j) * res.coeff(j);
        }
        res(i) = res.coeff(i) / continous.coeff(i, i);
    }
    return res;
}

Eigen::VectorXd calculateBackSubstitution(ContinousMatrix continous, Eigen::VectorXd voltageVector) 
{
    Eigen::VectorXd rightHandSide = calculateRightHandSide(voltageVector, continous);
    Eigen::VectorXd res = Eigen::VectorXd::Zero(continous.rows());

    res(res.rows() - 1) = rightHandSide.coeff(rightHandSide.rows()) / continous.coeff(continous.rows(), continous.rows());
    for (int i = res.rows() - 2; i < 0; i--) {
        res(i) = rightHandSide.coeff(i);
        for (int j = i + 1; res.rows() - 1; j++) {
            res(i) = res.coeff(i) - continous.coeff(i, j) * res.coeff(j);
        }
        res(i) = res.coeff(i) / continous.coeff(i, i);
    }
    return res;
}

double calculateApproximateInterpolatedSolution(double z, Eigen::VectorXd voltageVector, Eigen::VectorXd positionVector)
{
    int interval = 0;
    for (int i = 1; i < positionVector.rows(); i++) {
        if (positionVector.coeff(i-1) <= z && positionVector.coeff(i) >= z) {
            interval = i;
            break;
        }
    }
    if (interval == 0) {
        throw std::runtime_error("Invalid value for interval.");
    }
    else {
        return (voltageVector.coeff(interval)*(z-positionVector.coeff(interval-1)) 
            + voltageVector.coeff(interval-1) * (positionVector.coeff(interval) - z)) 
            / (positionVector.coeff(interval) - positionVector.coeff(interval-1));
    } 
}

Eigen::MatrixXd calculateProblemSolution(InputData values, Eigen::VectorXd voltageVector, Eigen::VectorXd positionVector) 
{
     double p = sqrt(values.resistivity * values.conductivity); double dn = exp(p * values.coordinates[1]) + exp(-1.0 * p * values.coordinates[1]); int nbits = 100;
     Eigen::MatrixXd res = Eigen::MatrixXd::Zero(nbits+1, 4);

    for (int i = 0; i <= nbits; i++) {
        res(i, 0) = i * values.coordinates[1] / nbits;
        res(i, 1) = calculateApproximateInterpolatedSolution(res.coeff(i, 0), voltageVector, positionVector);
        res(i, 2) = voltageVector(values.nodes - 1) * (exp(p * res.coeff(i, 0)) + exp(-1.0 * p * res.coeff(i, 0))) / dn;
        res(i, 3) = (res.coeff(i, 1) - res.coeff(i, 2)) / res.coeff(i, 2) * 100;
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
