#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>

#include <map>

using namespace Eigen;
using namespace std;

class MeshDeformation
{

public:
    /** 
     * Initialize the data structures
     **/
    MeshDeformation(MatrixXd &V_original, MatrixXi &F_original)
    {
        // double c = cos(2.11);
        // double s = sin(2.11);
        // Rotation << c, -s, 0,
        //             s, c, 0,
        //             0, 0, 1;
        // std::cout << "Rotation:\n" << Rotation << std::endl << std::endl;

        std::cout << "Starting initialization..." << std::endl;

        V0 = &V_original;
        F = &F_original;

        // V1 = MatrixXd(nVertices, 3);
        V1 = *V0;
        
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


        MatrixXd change = MatrixXd(1, 3);

        /* Constraines for the bar */

        // change << 60., -75., 0.;
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
        //     if (V1(i, 1) > 145) {
        //         constraints[i] = V0->row(i) + change;
        //     }

        //     if (V1(i, 1) < 6) {
        //         constraints[i] = V0->row(i);
        //     }
        // }

        /* Constraines for the square */

        // change << 0., -0.5, -0.75;
        change << 0., -0.75, 0.;
        for (int i = 0; i < nVertices; ++i) {
            if (V1(i, 1) > 0.9) {
                constraints[i] = V0->row(i) + change;
            }

            if (V1(i, 1) < 0.1) {
                constraints[i] = V0->row(i);
            }
        }

        rightSide.resize(nVertices, 3);
        setTheLeftSide();
        sparse_solver.analyzePattern(leftSide);
        sparse_solver.factorize(leftSide);

        std::cout << "Initialized solver" << std::endl;
    }

    void performOneIteration() {
        // if (iterations++ > 0) {
            estimateRotations();
        // }

        estimatePositions();
    }

    /**
     * Return the vertices
     **/
    MatrixXd getVertexCoordinates()
    {
        return *V0;
    }

    MatrixXd getEstimatedVertexCoordinates()
    {
        return V1;
    }

    /** 
     * Return the faces
     **/
    MatrixXi getFaces()
    {
        return *F;
    }

    /** 
     * Print the combinatorial information of the subdivided mesh <b>
     * verbosity=0: print only the number of vertices and faces <b>
     * verbosity=1: print all incidence relations
     **/
    void print(int verbosity)
    {
        cout << "\tn=" << nVertices << ", f=" << nFaces << endl;

        if (verbosity > 0) // print all vertex coordinates and face/vertex incidence relations
        {
            for (int i = 0; i < nVertices; i++)
            {
                cout << "v" << i << ": " << V1.row(i) << endl;
            }

            std::cout << "new faces: " << nFaces << endl;
            for (int i = 0; i < nFaces; i++)
            {
                cout << "f" << i << ": " << F->row(i) << endl;
            }
        }
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
        std::cout << "Computed" << std::endl;
    }

    void estimateRotations() {
        std::cout << "Computing the estimation rotations..." << std::endl;
        for (int v_index = 0; v_index < nVertices; ++v_index) {
            Matrix3d C = Matrix3d::Zero();
            Matrix3d id = Matrix3d::Identity();

            for (int i = 0; i < LB.outerSize(); ++i) {

                for (SparseMatrix<double>::InnerIterator it(LB, i); it; ++it)
                {
                    int j = it.row();

                    // count each edge only once and skip zeros
                    if (it.value() != 0 && j > i) {
                        C += (V0->row(j) - V0->row(i)).transpose() * (V1.row(j) - V1.row(i)) * it.value();
                    }
                }
            }

            JacobiSVD<MatrixXd> svd(C, ComputeFullU | ComputeFullV);
            id(2, 2) = (svd.matrixU() * svd.matrixV().transpose()).determinant();
            Matrix3d R = svd.matrixV() * id * svd.matrixU().transpose();

            Rs[v_index] = R;

            // std::cout << "Rotation for " << v_index << " :\n" << R << std::endl << std::endl;
        }
        std::cout << "Computed" << std::endl;
    }

    void estimatePositions() {
        std::cout << "Computing the estimation positions..." << std::endl;
        V1.setZero();
        computeTheRightSide();

        V1.col(0) = sparse_solver.solve(rightSide.col(0));
        // std::cout << V1 << std::endl << endl;
        V1.col(1) = sparse_solver.solve(rightSide.col(1));
        // std::cout << V1 << std::endl << endl;
        V1.col(2) = sparse_solver.solve(rightSide.col(2));
        // std::cout << V1 << std::endl << endl;

        std::cout << "Computed" << std::endl;
    }

    SparseMatrix<double> LB; 	// descrete Laplace-Beltrami operator
    SparseMatrix<double> leftSide;
    MatrixXd rightSide;
    SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > sparse_solver;

    std::vector<Matrix3d> Rs;
    // Matrix3d Rotation;

    int nVertices, nFaces;   	// number of vertices, faces in the new subdivided mesh
    MatrixXd *V0; 				// vertex coordinates of the original input mesh
    MatrixXi *F;
    MatrixXd V1;		   		// vertex coordinates of the new subdivided mesh

    int iterations = 0;
    map<int, Vector3d> constraints;
};
