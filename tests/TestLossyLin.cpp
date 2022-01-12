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
    for (auto n = 0; n < 10; n++) {
        ASSERT_NO_THROW(ConnectionMatrix(n));
    }

    auto maxNodes = 5;
    ConnectionMatrix connection(maxNodes);

    EXPECT_EQ(maxNodes, connection.rows());
    EXPECT_EQ(2 * maxNodes - 2, connection.cols());
    EXPECT_EQ(0.0, connection.coeff(1, 5));
    EXPECT_EQ(1.0, connection.coeff(1, 0));
}

TEST(lossylin, DiscontinuousMatrix_ctor)
{
    InputData values(readInputData("./testsData/data.json"));
    
    //for (auto n = 0; n < 10; n++) {
    //    ASSERT_NO_THROW(DiscontinousMatrix(n, values));
    //}
    
    auto maxNodes = 5;
    DiscontinousMatrix discontinous(maxNodes, values);

    EXPECT_NEAR(2.083, discontinous.coeff(0, 0), 0.005);
    EXPECT_NEAR(-1.958, discontinous.coeff(0, 1), 0.005);

}

//TEST(lossylin, DISABLED_FetchVoltageValue_Test)
//{
//    EXPECT_NEAR(FetchVoltageValue(0.0),0.265802,0.05);
//    EXPECT_NEAR(FetchVoltageValue(0.2),0.271136,0.05);
//    EXPECT_NEAR(FetchVoltageValue(1.90),0.908441,0.05);
//    EXPECT_NEAR(FetchVoltageValue(2.0),1.0,0.05);
//}
