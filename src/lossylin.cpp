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



}
