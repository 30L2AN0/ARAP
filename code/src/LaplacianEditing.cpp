#include <igl/opengl/glfw/Viewer.h>
#include <igl/min_quad_with_fixed.h>

// #ifndef HALFEDGE_DS_HEADER
// #define HALFEDGE_DS_HEADER
// #include "HalfedgeDS.cpp"
// #include "HalfedgeBuilder.cpp"
// #endif

#include <map>

Eigen::MatrixXd static laplacianEditing(const Eigen::MatrixXd &vertices,
                                        // const Eigen::MatrixXi &faces,            // uncomment if run with he
                                        const Eigen::VectorXi &constraintsIds,
                                        const Eigen::MatrixXd &constraints,
                                        const Eigen::SparseMatrix<double> &LB)
{
    // std::cout << "started laplacianEditing" << std::endl;

    const int nbVertices = vertices.rows();
    Eigen::MatrixXd deltas(nbVertices, 3);
    std::vector<Eigen::Triplet<double>> differentialCoordinatesTriplets;
    differentialCoordinatesTriplets.reserve(6 * nbVertices);


    /* run with half edges */

    // HalfedgeBuilder *builder = new HalfedgeBuilder();
    // HalfedgeDS halfEdges = (builder->createMeshWithFaces(vertices.rows(), faces));

    // for (size_t vertexId = 0; vertexId < nbVertices; vertexId++)
    // {
    //     std::list<int> neighborsIds;
    //     int halfEdge0 = halfEdges.getNext(halfEdges.getEdge(vertexId));
    //     double valency = 1;
    //     int neighborId = halfEdges.getTarget(halfEdge0);
    //     neighborsIds.push_back(neighborId);
    //     Eigen::RowVector3d sumNeighbors = vertices.row(neighborId);
    //     int halfEdge = halfEdges.getNext(halfEdges.getOpposite(halfEdge0));
    //     while (halfEdge != halfEdge0)
    //     {
    //         valency++;
    //         neighborId = halfEdges.getTarget(halfEdge);
    //         neighborsIds.push_back(neighborId);
    //         sumNeighbors += vertices.row(neighborId);
    //         halfEdge = halfEdges.getNext(halfEdges.getOpposite(halfEdge));
    //     }
    //     for (const auto &neighborId : neighborsIds)
    //     {
    //         differentialCoordinatesTriplets.push_back(Eigen::Triplet<double>(vertexId, neighborId, -1 / valency));
    //     }
    //     deltas.row(vertexId) = vertices.row(vertexId) - sumNeighbors / valency;
    // }

    // delete builder;


    /* run without half edges */

    for (int vertexId = 0; vertexId < LB.outerSize(); ++vertexId)
    {
        double valency = 0;
        std::list<int> neighborsIds;
        Eigen::RowVector3d sumNeighbors = Eigen::RowVector3d::Zero();

        for (Eigen::SparseMatrix<double>::InnerIterator it(LB, vertexId); it; ++it)
        {
            int neighborId = it.row();
            if (vertexId != neighborId && it.value() != 0)
            {
                ++valency;
                neighborsIds.push_back(neighborId);
                sumNeighbors += vertices.row(neighborId);
            }
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

    Eigen::VectorXd B;
    Eigen::MatrixXd newVertices(nbVertices, 3);
    for (int i = 0; i < vertices.cols(); ++i)
    {
        B = BBase * vertices.col(i);
        Eigen::VectorXd newVerticesColumn(nbVertices);
        igl::min_quad_with_fixed_data<double> mqwf;
        igl::min_quad_with_fixed_precompute(A, constraintsIds, Aeq, true, mqwf);
        igl::min_quad_with_fixed_solve(mqwf, B, constraints.col(i), Beq, newVerticesColumn);
        newVertices.col(i) = newVerticesColumn;
    }

    return newVertices;
}

Eigen::MatrixXd static laplacianEditingFromMap(const Eigen::MatrixXd &vertices,
                                               const Eigen::MatrixXi &faces,
                                               const std::map<int, Eigen::Vector3d> &constraints_map,
                                               const Eigen::SparseMatrix<double> &LB)
{
    size_t nbConstraints = constraints_map.size();
    Eigen::VectorXi constraintsIds(nbConstraints);
    Eigen::MatrixXd constraintsValues(nbConstraints, 3);

    size_t i = 0;
    for (const auto &[point, value] : constraints_map)
    {
        constraintsIds(i) = point;
        constraintsValues.row(i) = value;
        i++;
    }

    // std::cout << "Calling original laplacianEditing" << std::endl;

    // return laplacianEditing(vertices, faces, constraintsIds, constraintsValues, LB);     // choose if run with he
    return laplacianEditing(vertices, constraintsIds, constraintsValues, LB);     // choose if run without he
}