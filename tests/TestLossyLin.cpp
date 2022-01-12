#include <gtest/gtest.h>

#include "lossylin.h"

using namespace std;

typedef vector<double> Coordinates;

const int maxNodes = 30;

TEST(lossylin, readInputData)
{
    auto j = FEM::readInputData();
    
    auto jVoltage = j.at("voltage");
    ASSERT_TRUE(jVoltage != nullptr);
    EXPECT_EQ(2.0, jVoltage.get<double>());
    
    auto jCoords = j.at("coordinates");
    ASSERT_TRUE(jCoords != nullptr);
    Coordinates coords = jCoords.get<Coordinates>();
    ASSERT_EQ(2, coords.size());
    EXPECT_EQ(0.0, coords[0]);
    EXPECT_EQ(2.0, coords[1]);
}

TEST(lossylin, readSource)
{
    auto j = FEM::readInputData();
    ASSERT_TRUE(j != nullptr);

    auto jSource = j.at("source");
    ASSERT_TRUE(jSource != nullptr);

    FEM::Source source;
    ASSERT_NO_THROW(source = FEM::Source(jSource));

    FEM::Values values(j);
    ASSERT_NO_THROW(values = FEM::Values(j));

    EXPECT_EQ(values.voltage, 2.0);

    //funcionQueUsaFuente(source);
}

TEST(lossylin, buildMatrixesAndRightHandSideTerm)
{
    Eigen::MatrixXd connectionMatrix = FEM::buildConnectionMatrix(maxNodes);
    ASSERT_TRUE(&connectionMatrix != nullptr);
    EXPECT_EQ(maxNodes, connectionMatrix.rows());
    EXPECT_EQ(2 * maxNodes - 2, connectionMatrix.cols());
    EXPECT_EQ(0.0, connectionMatrix.coeff(1, 15));
    EXPECT_EQ(1.0, connectionMatrix.coeff(1, 0));

    Eigen::MatrixXd discontinousMatrix = FEM::buildDiscontinousMatrix(maxNodes);
    ASSERT_TRUE(&discontinousMatrix != nullptr);
    EXPECT_NEAR(2.083, discontinousMatrix.coeff(0, 0), 0.005);
    EXPECT_NEAR(-1.958, discontinousMatrix.coeff(0, 1), 0.005);

}




//TEST(lossylin, DISABLED_FetchVoltageValue_Test)
//{
//    EXPECT_NEAR(FetchVoltageValue(0.0),0.265802,0.05);
//    EXPECT_NEAR(FetchVoltageValue(0.2),0.271136,0.05);
//    EXPECT_NEAR(FetchVoltageValue(1.90),0.908441,0.05);
//    EXPECT_NEAR(FetchVoltageValue(2.0),1.0,0.05);
//}
