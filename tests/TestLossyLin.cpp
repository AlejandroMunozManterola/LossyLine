#include <gtest/gtest.h>

#include "lossylin.h"

using namespace std;

typedef vector<double> Coordinates;

const int maxnodes = 30;
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

    auto jSource = j.at("source");
    ASSERT_TRUE(jSource != nullptr);

    FEM::Source source;
    ASSERT_NO_THROW(source = FEM::Source(jSource));

    //funcionQueUsaFuente(source);
}

TEST(lossylin, initialiseMatrixesAndRightHandSideTerm)
{
    Eigen::MatrixXd continousMatrix = FEM::initialiseContinousMatrix(maxnodes);
    ASSERT_TRUE(&continousMatrix != nullptr);
    EXPECT_EQ(maxnodes, continousMatrix.rows());
    EXPECT_EQ(0.0, continousMatrix.coeff(0,0));
    EXPECT_EQ(0.0, continousMatrix.coeff(continousMatrix.rows()-1,continousMatrix.cols()-1));

    Eigen::MatrixXd discontinousMatrix = FEM::initialiseDiscontinousMatrix(maxnodes);
    ASSERT_TRUE(&discontinousMatrix != nullptr);
    EXPECT_EQ(2*maxnodes-2, discontinousMatrix.rows());
    EXPECT_EQ(0.0, discontinousMatrix.coeff(0, 0));
    EXPECT_EQ(0.0, discontinousMatrix.coeff(discontinousMatrix.rows()-1, discontinousMatrix.cols()-1));

    Eigen::VectorXd rightHandSideTerm = FEM::initialiseRightHandSideTerm(maxnodes);
    ASSERT_TRUE(&rightHandSideTerm != nullptr);
    EXPECT_EQ(maxnodes, rightHandSideTerm.size());
    EXPECT_EQ(0.0, rightHandSideTerm.coeff(0));
    EXPECT_EQ(0.0, rightHandSideTerm.coeff(rightHandSideTerm.rows()-1));

}

TEST(lossylin, buildMatrixesAndRightHandSideTerm)
{
    Eigen::MatrixXd connectionMatrix = FEM::buildConnectionMatrix(maxnodes);
    ASSERT_TRUE(&connectionMatrix != nullptr);
    EXPECT_EQ(maxnodes, connectionMatrix.rows());
    EXPECT_EQ(2 * maxnodes - 2, connectionMatrix.cols());
    EXPECT_EQ(0.0, connectionMatrix.coeff(1, 15));
    EXPECT_EQ(1.0, connectionMatrix.coeff(1, 0));
}




//TEST(lossylin, DISABLED_FetchVoltageValue_Test)
//{
//    EXPECT_NEAR(FetchVoltageValue(0.0),0.265802,0.05);
//    EXPECT_NEAR(FetchVoltageValue(0.2),0.271136,0.05);
//    EXPECT_NEAR(FetchVoltageValue(1.90),0.908441,0.05);
//    EXPECT_NEAR(FetchVoltageValue(2.0),1.0,0.05);
//}
