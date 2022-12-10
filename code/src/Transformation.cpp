#include <igl/opengl/glfw/Viewer.h>
#include <igl/min_quad_with_fixed.h>

#ifndef HALFEDGE_DS_HEADER
#define HALFEDGE_DS_HEADER
#include "HalfedgeDS.cpp"
#endif

Eigen::MatrixXd static laplacianEditing(const Eigen::MatrixXd &vertices,
                                        const HalfedgeDS &halfEdges,
                                        const Eigen::VectorXi &constraintsIds,
                                        const Eigen::MatrixXd &constraints) // use laplac Bletrami matrix to find neighbors
{
    const int nbVertices = vertices.rows();
    std::cout << "nbVertices: " << std::endl
              << nbVertices << std::endl
              << std::endl;
    Eigen::MatrixXd deltas(nbVertices, 3);
    std::vector<Eigen::Triplet<double>> differentialCoordinatesTriplets;
    differentialCoordinatesTriplets.reserve(6 * nbVertices);
    std::cout << "differentialCoordinatesTriplets: " << std::endl
              << differentialCoordinatesTriplets.size() << std::endl
              << std::endl;
    for (size_t vertexId = 0; vertexId < nbVertices; vertexId++)
    {
        std::list<int> neighborsIds;
        int halfEdge0 = halfEdges.getNext(halfEdges.getEdge(vertexId));
        double valency = 1;
        int neighborId = halfEdges.getTarget(halfEdge0);
        neighborsIds.push_back(neighborId);
        Eigen::RowVector3d sumNeighbors = vertices.row(neighborId);
        int halfEdge = halfEdges.getNext(halfEdges.getOpposite(halfEdge0));
        while (halfEdge != halfEdge0)
        {
            valency++;
            neighborId = halfEdges.getTarget(halfEdge);
            neighborsIds.push_back(neighborId);
            sumNeighbors += vertices.row(neighborId);
            halfEdge = halfEdges.getNext(halfEdges.getOpposite(halfEdge));
        }
        std::cout << "valency: " << std::endl
                  << valency << std::endl
                  << std::endl;
        for (const auto &neighborId : neighborsIds)
        {
            differentialCoordinatesTriplets.push_back(Eigen::Triplet<double>(neighborId, vertexId, -1 / valency));
        }
        deltas.row(vertexId) = vertices.row(vertexId) - sumNeighbors / valency;
    }
    Eigen::SparseMatrix<double> differentialCoordinates(nbVertices, nbVertices);
    differentialCoordinates.setFromTriplets(differentialCoordinatesTriplets.begin(), differentialCoordinatesTriplets.end());
    Eigen::SparseMatrix<double> identity(nbVertices, nbVertices);
    identity.setIdentity();
    differentialCoordinates += identity;
    // TEEEEEST
    differentialCoordinates = differentialCoordinates.transpose();
    // ~TEEEEEST
    std::cout << "differentialCoordinates: " << std::endl
              << MatrixXd(differentialCoordinates) << std::endl
              << std::endl;

    const Eigen::SparseMatrix<double> transposedDifferentialCoordinates = differentialCoordinates.transpose();
    std::cout << "transposedDifferentialCoordinates: " << std::endl
              << MatrixXd(transposedDifferentialCoordinates) << std::endl
              << std::endl;
    const Eigen::SparseMatrix<double> A = 2 * transposedDifferentialCoordinates * differentialCoordinates;
    std::cout << "A: " << std::endl
              << MatrixXd(A) << std::endl
              << std::endl;
    const Eigen::MatrixXd BBase = -2 * transposedDifferentialCoordinates * differentialCoordinates;
    std::cout << "BBase: " << std::endl
              << BBase << std::endl
              << std::endl;
    Eigen::SparseMatrix<double> Aeq(nbVertices, nbVertices);
    Aeq.setZero();
    std::cout << "Aeq: " << std::endl
              << MatrixXd(Aeq) << std::endl
              << std::endl;
    const Eigen::VectorXd Beq = Eigen::VectorXd::Zero(nbVertices);
    std::cout << "Beq: " << std::endl
              << Beq << std::endl
              << std::endl;

    Eigen::VectorXd B = BBase * vertices.col(0);
    std::cout << "B: " << std::endl
              << B << std::endl
              << std::endl;
    Eigen::VectorXd newVerticesX(nbVertices);
    std::cout << "newVerticesX: " << std::endl
              << newVerticesX << std::endl
              << std::endl;
    igl::min_quad_with_fixed_data<double> mqwf;
    std::cout << "OK1" << std::endl
              << std::endl
              << std::endl;
    igl::min_quad_with_fixed_precompute(A, constraintsIds, Aeq, true, mqwf);
    std::cout << "OK2" << std::endl
              << std::endl
              << std::endl;
    igl::min_quad_with_fixed_solve(mqwf, B, constraints.col(0), Beq, newVerticesX);
    std::cout << "newVerticesX: " << std::endl
              << newVerticesX << std::endl
              << std::endl;

    B = BBase * vertices.col(1);
    std::cout << "B: " << std::endl
              << B << std::endl
              << std::endl;
    Eigen::VectorXd newVerticesY(nbVertices);
    std::cout << "newVerticesY: " << std::endl
              << newVerticesY << std::endl
              << std::endl;
    igl::min_quad_with_fixed_precompute(A, constraintsIds, Aeq, true, mqwf);
    std::cout << "OK" << std::endl
              << std::endl
              << std::endl;
    igl::min_quad_with_fixed_solve(mqwf, B, constraints.col(1), Beq, newVerticesY);
    std::cout << "newVerticesY: " << std::endl
              << newVerticesY << std::endl
              << std::endl;

    B = BBase * vertices.col(2);
    std::cout << "B: " << std::endl
              << B << std::endl
              << std::endl;
    Eigen::VectorXd newVerticesZ(nbVertices);
    std::cout << "newVerticesZ: " << std::endl
              << newVerticesZ << std::endl
              << std::endl;
    igl::min_quad_with_fixed_precompute(A, constraintsIds, Aeq, true, mqwf);
    std::cout << "OK" << std::endl
              << std::endl
              << std::endl;
    igl::min_quad_with_fixed_solve(mqwf, B, constraints.col(2), Beq, newVerticesZ);
    std::cout << "newVerticesZ: " << std::endl
              << newVerticesZ << std::endl
              << std::endl;

    Eigen::MatrixXd newVertices(nbVertices, 3);
    newVertices << newVerticesX, newVerticesY, newVerticesZ;
    std::cout << "newVertices: " << std::endl
              << newVertices << std::endl
              << std::endl;
    return newVertices;
}