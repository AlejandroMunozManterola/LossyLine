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

TEST(lossylin, positionVectorFunction) {

    InputData values(readInputData("./testsData/data.json"));

    Eigen::VectorXd positionVector = buildPositionVector(values.nodes, values.coordinates);

    EXPECT_EQ(0.0, positionVector.coeff(0));
    EXPECT_EQ(0.5, positionVector.coeff(1));
    EXPECT_EQ(1.0, positionVector.coeff(2));
    EXPECT_EQ(1.5, positionVector.coeff(3));
    EXPECT_EQ(2.0, positionVector.coeff(positionVector.rows()-1));

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

TEST(lossylin, calculateRightHandSide)
{
    InputData values(readInputData("./testsData/data.json"));
    ConnectionMatrix connection(values.nodes);
    DiscontinousMatrix discontinous(values);
    ContinousMatrix continous(values.nodes, connection, discontinous);
    Eigen::VectorXd rightHandSide = calculateRightHandSide(buildVoltageVector(values.nodes, values.voltage), continous);

    ASSERT_EQ(rightHandSide.rows(), continous.rows());
    EXPECT_EQ(0.0, rightHandSide.coeff(0));
    EXPECT_NEAR(-1.0 * values.voltage * -1.953, rightHandSide.coeff(rightHandSide.rows() - 2), 0.05);
    EXPECT_NEAR(-1.0 * values.voltage * 2.083, rightHandSide.coeff(rightHandSide.rows() - 1), 0.05);
}

TEST(lossylin, calculateProblemSolution)
{
    InputData values(readInputData("./testsData/data.json"));
    Eigen::VectorXd voltageVector = buildVoltageVector(values.nodes, values.voltage);
    Eigen::VectorXd positionVector = buildPositionVector(values.nodes, values.coordinates);
    Eigen::MatrixXd solution = calculateProblemSolution(values, voltageVector, positionVector);
 
    EXPECT_EQ(0.0, solution.coeff(0, 0));
    EXPECT_EQ(2.0, solution.coeff(solution.rows() - 1, 0));
    EXPECT_NEAR(0.262344, solution.coeff(0, 1), 0.05);
    EXPECT_NEAR(1.0, solution.coeff(solution.rows() - 1, 1), 0.05);
    EXPECT_NEAR(0.265802, solution.coeff(0, 2), 0.05);
    EXPECT_NEAR(1.0, solution.coeff(solution.rows()-1, 2), 0.05);
}
