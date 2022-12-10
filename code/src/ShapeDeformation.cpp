#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>

#include <map>

// #include "LaplacianEditing.cpp"

using namespace Eigen;
using namespace std;

class MeshDeformation
{

public:
    /** 
     * Initialize the data structures
     **/
    MeshDeformation(MatrixXd &V_original, MatrixXi &F_original, MatrixXd &V_initial_guess)
    {
        std::cout << "Starting initialization..." << std::endl;

        V0 = &V_original;
        F = &F_original;
        V1 = &V_initial_guess;

        // V1 = MatrixXd(nVertices, 3);
        // *V1 = *V0;
        
        nVertices = V_original.rows();
        nFaces = F_original.rows();

        std::cout << "Initialized shapes" << std::endl;

        // Compute the descrete Laplace-Beltrami operator
        igl::cotmatrix(V_original, F_original, LB);
        // std::cout << "laplacian:\n" << MatrixXd(LB) << std::endl << std::endl;
        // int zeroes = 0;
        // std::cout << "zeroes: " << (LB.size() - LB.nonZeros() + zeroes) / LB.size() * 100 << "%" << std::endl;

        std::cout << "Initialized Laplacian-Beltrami" << std::endl;

        // Initialize initial rotations
        Rs.resize(nVertices, Matrix3d::Identity());


        // MatrixXd change = MatrixXd(1, 3);

        /* Constraines for the bar */

        // change << 75., -90., 0.;
        // double cs = 1 / sqrt(2) * 100;
        // change << cs, -cs, 0.;
        // int i = 0;
        // while ((*V0)(i, 1) == 0) {
        //     constraints[i] = V0->row(i);
        //     ++i;
        // }
        
        // while ((*V0)(i, 1) > 150) {
        //     constraints[i] = V0->row(i) + change;
        //     ++i;
        // }

        // for (int i = 0; i < nVertices; ++i) {
        //     if ((*V1)(i, 1) > 140) {
        //         constraints[i] = V0->row(i) + change;
        //     }

        //     if ((*V1)(i, 1) < 8) {
        //         constraints[i] = V0->row(i);
        //     }
        // }


        /* Constraines for the square */

        // change << 0., -0.5, -0.75;
        // // change << 0., -0.75, 0.;
        // for (int i = 0; i < nVertices; ++i) {
        //     if ((*V1)(i, 1) > 0.9) {
        //         constraints[i] = V0->row(i) + change;
        //     }

        //     if ((*V1)(i, 1) < 0.1) {
        //         constraints[i] = V0->row(i);
        //     }
        // }

        rightSide.resize(nVertices, 3);
        // setTheLeftSide();
        // sparse_solver.analyzePattern(leftSide);
        // sparse_solver.factorize(leftSide);
        // std::cout << "Initialized solver" << std::endl;
    }

    void performOneIteration() {
        // if (iterations++ > 0) {
            estimateRotations();
        // }
        estimatePositions();
    }

    void setConstraints(const VectorXi &constraintsIds,
                        const MatrixXd &new_constraints) {
        constraints.clear();
        for (int i = 0; i < constraintsIds.size(); ++i) {
            constraints[constraintsIds[i]] = new_constraints.row(i);
        }

        setTheLeftSide();
        sparse_solver.analyzePattern(leftSide);
        sparse_solver.factorize(leftSide);
        std::cout << "Initialized solver" << std::endl;
    }

    /**
     * Return the initial vertices
     **/
    MatrixXd getVertexCoordinates()
    {
        return *V0;
    }

    /**
     * Return the new vertices
     **/
    MatrixXd getEstimatedVertexCoordinates()
    {
        return *V1;
    }

    /** 
     * Return the faces
     **/
    MatrixXi getFaces()
    {
        return *F;
    }

private:
    void setTheLeftSide() {
        leftSide = LB;
        for (const auto& [point, new_value] : constraints) {
            leftSide.row(point) *= 0;
            leftSide.coeffRef(point, point) = 1;
        }

        // std::cout << "The left side of points estimation equation:\n" << MatrixXd(leftSide) << std::endl << std::endl;
    }

    void computeTheRightSide() {
        std::cout << "Computing the right side..." << std::endl;
        rightSide.setZero();

        for (int i = 0; i < LB.outerSize(); ++i) {

            MatrixXd b_i = MatrixXd::Zero(1, 3);

            for (SparseMatrix<double>::InnerIterator it(LB, i); it; ++it)
            {
                int j = it.row();

                // skip zeros
                if (it.value() != 0 && i != j) {
                    MatrixXd temp = it.value() / 2. * (V0->row(j) - V0->row(i)) * (Rs[i] + Rs[j]).transpose();
                    // std::cout << "see for j = " << j << " is: " << temp << std::endl;
                    b_i += temp;
                }
            }

            rightSide.row(i) = b_i;
            // std::cout << "sum for i " << i << " is: " << b_i << std::endl << std::endl;
        }

        for (const auto& [point, new_value] : constraints) {
            rightSide.row(point) = new_value;
        }

        // std::cout << "The right side is:\n" << rightSide << std::endl;
        std::cout << "Computed the right side" << std::endl;
    }

    void estimateRotations() {
        std::cout << "Computing the estimation rotations..." << std::endl;
        for (int i = 0; i < LB.outerSize(); ++i) {
            Matrix3d C = Matrix3d::Zero();
            Matrix3d id = Matrix3d::Identity();

                for (SparseMatrix<double>::InnerIterator it(LB, i); it; ++it)
                {
                    int j = it.row();

                    // skip zeros
                    if (it.value() != 0 && j != i) {
                        C += (V0->row(j) - V0->row(i)).transpose() * (V1->row(j) - V1->row(i)) * it.value();
                    }
                }

            JacobiSVD<MatrixXd> svd(C, ComputeFullU | ComputeFullV);
            id(2, 2) = (svd.matrixU() * svd.matrixV().transpose()).determinant();
            Matrix3d R = svd.matrixV() * id * svd.matrixU().transpose();

            Rs[i] = R;

            // std::cout << "Rotation for " << i << " :\n" << R << std::endl << std::endl;
        }
        std::cout << "Computed rotations" << std::endl;
    }

    void estimatePositions() {
        V1->setZero();
        computeTheRightSide();
        std::cout << "Computing the estimation positions..." << std::endl;

        V1->col(0) = sparse_solver.solve(rightSide.col(0));
        // std::cout << V1 << std::endl << std::endl;
        V1->col(1) = sparse_solver.solve(rightSide.col(1));
        // std::cout << V1 << std::endl << std::endl;
        V1->col(2) = sparse_solver.solve(rightSide.col(2));
        // std::cout << V1 << std::endl << std::endl;

        std::cout << "Computed positions" << std::endl;
    }

    SparseMatrix<double> LB; 	// descrete Laplace-Beltrami operator
    SparseMatrix<double> leftSide;
    MatrixXd rightSide;
    SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > sparse_solver;

    std::vector<Matrix3d> Rs;

    int nVertices, nFaces;   	// number of vertices, faces in the new subdivided mesh
    MatrixXd *V0; 				// vertex coordinates of the original input mesh
    MatrixXi *F;
    MatrixXd *V1;		   		// vertex coordinates of the new subdivided mesh

    int iterations = 0;
    map<int, Vector3d> constraints;
};
