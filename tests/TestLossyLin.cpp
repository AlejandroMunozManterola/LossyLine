#include <gtest/gtest.h>

#include "lossylin.h"

using namespace std;
using namespace FEM;

typedef vector<double> Coordinates;

TEST(lossylin, readInputData)
{
    ASSERT_NO_THROW(InputData(readInputData("./testsData/data.json")));
    
    InputData values(readInputData("./testsData/data.json"));
    
    EXPECT_EQ(2.0, values.voltage);
    
    ASSERT_EQ(2, values.coordinates.size());
    EXPECT_EQ(0.0, values.coordinates[0]);
    EXPECT_EQ(2.0, values.coordinates[1]);
    
    EXPECT_EQ(values.voltage, 2.0);
}

TEST(lossylin, ConnectionMatrix_ctor)
{

    for (auto n = 2; n < 10; n++) {
        ASSERT_NO_THROW(ConnectionMatrix conn(n));
    }

    InputData values(readInputData("./testsData/data.json"));
    ConnectionMatrix connection(values.nodes);

    EXPECT_EQ(values.nodes, connection.rows());
    EXPECT_EQ(2 * values.nodes - 2, connection.cols());
    EXPECT_EQ(0.0, connection.coeff(0, connection.rows()));
    EXPECT_EQ(1.0, connection.coeff(2, 3));
}

TEST(lossylin, DiscontinuousMatrix_ctor)
{
    InputData values(readInputData("./testsData/data.json"));
    
    ASSERT_NO_THROW(DiscontinousMatrix disc(values));

    DiscontinousMatrix discontinous(values);

    EXPECT_NEAR(2.083, discontinous.coeff(0, 0), 0.005);
    EXPECT_NEAR(-1.958, discontinous.coeff(0, 1), 0.005);
    EXPECT_NEAR(2.083, discontinous.coeff(discontinous.rows() - 1, discontinous.cols() - 1),0.05);

}

TEST(lossylin, ContinousMatrix_ctor)
{
    InputData values(readInputData("./testsData/data.json"));
    ConnectionMatrix connection(values.nodes);
    DiscontinousMatrix discontinous(values);

    ASSERT_NO_THROW(ContinousMatrix disc(values.nodes, connection, discontinous));

    ContinousMatrix continous(values.nodes, connection, discontinous);

    EXPECT_NEAR(2.083, continous.coeff(0, 0),0.05);
    EXPECT_NEAR(-1.953, continous.coeff(0, 1),0.05);
    EXPECT_NEAR(4.166, continous.coeff(1, 1), 0.05);
    EXPECT_NEAR(2.083, continous.coeff(continous.rows()-1,continous.cols()-1), 0.05);

}

TEST(lossylin, transposeVoltageToRightSide)
{
    InputData values(readInputData("./testsData/data.json"));
    ConnectionMatrix connection(values.nodes);
    DiscontinousMatrix discontinous(values);
    ContinousMatrix continous(values.nodes, connection, discontinous);
    Eigen::VectorXd voltageVector = buildVoltageVector(values.nodes, values.voltage);
    Eigen::VectorXd rightHandSide = transposeVoltageToRightSide(voltageVector, continous);

    ASSERT_EQ(voltageVector.rows(), continous.rows());
    EXPECT_EQ(0.0, rightHandSide.coeff(0));
    EXPECT_NEAR(-1.0 * values.voltage * -1.953, rightHandSide.coeff(rightHandSide.rows() - 2), 0.05);
    EXPECT_NEAR(-1.0 * values.voltage * 2.083, rightHandSide.coeff(rightHandSide.rows() - 1), 0.05);
}

//TEST(lossylin, DISABLED_FetchVoltageValue_Test)
//{
//    EXPECT_NEAR(FetchVoltageValue(0.0),0.265802,0.05);
//    EXPECT_NEAR(FetchVoltageValue(0.2),0.271136,0.05);
//    EXPECT_NEAR(FetchVoltageValue(1.90),0.908441,0.05);
//    EXPECT_NEAR(FetchVoltageValue(2.0),1.0,0.05);
//}
