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
    Eigen::MatrixXd deltas(nbVertices, 3);
    std::vector<Eigen::Triplet<double>> differentialCoordinatesTriplets;
    differentialCoordinatesTriplets.reserve(6 * nbVertices);
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
        for (const auto &neighborId : neighborsIds)
        {
            differentialCoordinatesTriplets.push_back(Eigen::Triplet<double>(vertexId, neighborId, -1 / valency));
        }
        deltas.row(vertexId) = vertices.row(vertexId) - sumNeighbors / valency;
    }
    Eigen::SparseMatrix<double> differentialCoordinates(nbVertices, nbVertices);
    differentialCoordinates.setFromTriplets(differentialCoordinatesTriplets.begin(), differentialCoordinatesTriplets.end());
    Eigen::SparseMatrix<double> identity(nbVertices, nbVertices);
    identity.setIdentity();
    differentialCoordinates += identity;

    const Eigen::SparseMatrix<double> transposedDifferentialCoordinates = differentialCoordinates.transpose();
    const Eigen::SparseMatrix<double> A = 2 * transposedDifferentialCoordinates * differentialCoordinates;
    const Eigen::MatrixXd BBase = -2 * transposedDifferentialCoordinates * differentialCoordinates;
    Eigen::SparseMatrix<double> Aeq(nbVertices, nbVertices);
    Aeq.setZero();
    const Eigen::VectorXd Beq = Eigen::VectorXd::Zero(nbVertices);

    Eigen::VectorXd B = BBase * vertices.col(0);
    Eigen::VectorXd newVerticesX(nbVertices);
    igl::min_quad_with_fixed_data<double> mqwf;
    igl::min_quad_with_fixed_precompute(A, constraintsIds, Aeq, true, mqwf);
    igl::min_quad_with_fixed_solve(mqwf, B, constraints.col(0), Beq, newVerticesX);

    B = BBase * vertices.col(1);
    Eigen::VectorXd newVerticesY(nbVertices);
    igl::min_quad_with_fixed_precompute(A, constraintsIds, Aeq, true, mqwf);
    igl::min_quad_with_fixed_solve(mqwf, B, constraints.col(1), Beq, newVerticesY);

    B = BBase * vertices.col(2);
    Eigen::VectorXd newVerticesZ(nbVertices);
    igl::min_quad_with_fixed_precompute(A, constraintsIds, Aeq, true, mqwf);
    igl::min_quad_with_fixed_solve(mqwf, B, constraints.col(2), Beq, newVerticesZ);

    Eigen::MatrixXd newVertices(nbVertices, 3);
    newVertices << newVerticesX, newVerticesY, newVerticesZ;
    return newVertices;
}