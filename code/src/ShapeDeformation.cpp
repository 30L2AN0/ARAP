
#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>

#include <map>

#include "LaplacianEditing.cpp"

using namespace Eigen;

class MeshDeformation
{

public:
    /**
     * Create the deformation object
     **/
    MeshDeformation(const MatrixXd &V_original, const MatrixXi &F_original, MatrixXd &V_guess)
        : V0(V_original), F(F_original), V1(V_guess)
    {
        std::cout << "As rigit as possible deformation class was created" << std::endl;
    }

    /**
     * Initialize the data structures
     **/
    void initialize(bool compute_initial_guess)
    {
        std::cout << "Starting arap deformation initialization..." << std::endl;
        nVertices = V0.rows();

        // Compute the descrete Laplace-Beltrami operator
        igl::cotmatrix(V0, F, LB);
        // std::cout << "laplacian:\n" << MatrixXd(LB) << std::endl << std::endl;
        std::cout << "Initialized Laplace-Beltrami" << std::endl;

        Rs.resize(nVertices, Matrix3d::Identity());

        if (compute_initial_guess)
        {
            initializeConstraints();
            initializeEquation();
            initializeFirstGuess();
        }

        initialized = true;
    }

    /**
     * Check if initialized
     **/
    bool isInitialized()
    {
        return initialized;
    }

    void setConstraints(const VectorXi &constraintsIds,
                        const MatrixXd &new_constraints)
    {
        constraints.clear();
        for (int i = 0; i < constraintsIds.size(); ++i)
        {
            constraints[constraintsIds[i]] = new_constraints.row(i);
        }

        initializeEquation();
    }

    void performOneIteration()
    {
        estimateRotations();
        estimatePositions();
    }

    /**
     * Return the initial vertices
     **/
    MatrixXd getVertexCoordinates()
    {
        return V0;
    }

    /**
     * Return the new vertices
     **/
    MatrixXd getEstimatedVertexCoordinates()
    {
        return V1;
    }

    /**
     * Return the faces
     **/
    MatrixXi getFaces()
    {
        return F;
    }

private:
    void initializeConstraints()
    {
        constraints.clear();
        MatrixXd change = MatrixXd(1, 3);

        /* Constraints for the bar2 */

        // change << 75., -90., 0.;
        // // double cs = 1 / sqrt(2) * 100;
        // // change << cs, -cs, 0.;
        // // int i = 0;
        // // while (V0(i, 1) == 0) {
        // //     constraints[i] = V0.row(i);
        // //     ++i;
        // // }

        // // while (V0(i, 1) > 150) {
        // //     constraints[i] = V0.row(i) + change;
        // //     ++i;
        // // }

        // for (int i = 0; i < nVertices; ++i) {
        //     if (V1(i, 1) > 140) {
        //         constraints[i] = V0.row(i) + change;
        //     }

        //     if (V1(i, 1) < 8) {
        //         constraints[i] = V0.row(i);
        //     }
        // }

        /* Constraints for the square */

        // change << 0., -0.5, -0.75;
        // // change << 0., -0.75, 0.;
        // for (int i = 0; i < nVertices; ++i) {
        //     if (V1(i, 1) > 0.9) {
        //         constraints[i] = V0.row(i) + change;
        //     }

        //     if (V1(i, 1) < 0.1) {
        //         constraints[i] = V0.row(i);
        //     }
        // }

        /* Constraints for the square with spikes */

        change << 0.5, 0.75, 0.;
        for (int i = 0; i < nVertices; ++i)
        {
            if (V1(i, 0) < 40.)
            {
                constraints[i] = V0.row(i) + change;
            }

            if (V1(i, 0) > 40.85)
            {
                constraints[i] = V0.row(i);
            }
        }
    }

    void initializeEquation()
    {
        setTheLeftSide();
        sparse_solver.analyzePattern(leftSide);
        sparse_solver.factorize(leftSide);
        rightSide.resize(nVertices, 3);
        std::cout << "Initialized solver" << std::endl;
    }

    void initializeFirstGuess()
    {
        V1 = laplacianEditingFromMap(V0, F, constraints, LB);
    }

    void setTheLeftSide()
    {
        leftSide = LB;
        for (const auto &[point, new_value] : constraints)
        {
            leftSide.row(point) *= 0;
            leftSide.coeffRef(point, point) = 1;
        }
    }

    void computeTheRightSide()
    {
        std::cout << "Computing the right side..." << std::endl;
        rightSide.setZero();

        for (int i = 0; i < LB.outerSize(); ++i)
        {

            for (SparseMatrix<double>::InnerIterator it(LB, i); it; ++it)
            {
                int j = it.row();

                // skip zeros
                if (it.value() != 0 && i != j)
                {
                    rightSide.row(i) += it.value() / 2. * (V0.row(j) - V0.row(i)) * (Rs[i] + Rs[j]).transpose();
                }
            }
        }

        for (const auto &[point, new_value] : constraints)
        {
            rightSide.row(point) = new_value;
        }

        std::cout << "Computed the right side" << std::endl;
    }

    void estimateRotations()
    {
        std::cout << "Computing the estimation rotations..." << std::endl;
        for (int i = 0; i < LB.outerSize(); ++i)
        {
            Matrix3d C = Matrix3d::Zero();
            Matrix3d id = Matrix3d::Identity();

            for (SparseMatrix<double>::InnerIterator it(LB, i); it; ++it)
            {
                int j = it.row();

                // skip zeros
                if (it.value() != 0 && j != i)
                {
                    C += (V0.row(j) - V0.row(i)).transpose() * (V1.row(j) - V1.row(i)) * it.value();
                }
            }

            JacobiSVD<MatrixXd> svd(C, ComputeFullU | ComputeFullV);
            id(2, 2) = (svd.matrixU() * svd.matrixV().transpose()).determinant();
            Matrix3d R = svd.matrixV() * id * svd.matrixU().transpose();

            Rs[i] = R;
        }
        std::cout << "Computed rotations" << std::endl;
    }

    void estimatePositions()
    {
        V1.setZero();
        computeTheRightSide();
        std::cout << "Computing the estimation positions..." << std::endl;

        V1.col(0) = sparse_solver.solve(rightSide.col(0));
        V1.col(1) = sparse_solver.solve(rightSide.col(1));
        V1.col(2) = sparse_solver.solve(rightSide.col(2));

        std::cout << "Computed positions" << std::endl;
    }

    const MatrixXd &V0;             // vertex coordinates of the original input mesh
    const MatrixXi &F;
    MatrixXd &V1;                   // vertex coordinates of the new subdivided mesh
    int nVertices;                  // number of vertices, faces

    SparseMatrix<double> LB;        // descrete Laplace-Beltrami operator
    SparseMatrix<double> leftSide;  // left side of the equation to estimate the points positions
    MatrixXd rightSide;             // right side of the equation to estimate the points positions
    SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> sparse_solver;

    std::vector<Matrix3d> Rs;

    map<int, Vector3d> constraints;
    bool initialized = false;
};
