
#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>

#include <map>
#include <string>

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
        // std::cout << "As rigit as possible deformation class was created" << std::endl;
    }

    /**
     * Initialize the data structures
     **/
    void initialize(bool compute_initial_guess, const std::string &filename)
    {
        // std::cout << "Starting arap deformation initialization..." << std::endl;
        nVertices = V0.rows();
        should_compute_initial_guess = compute_initial_guess;

        // Compute the descrete Laplace-Beltrami operator
        igl::cotmatrix(V0, F, LB);
        // std::cout << "laplacian:\n" << MatrixXd(LB) << std::endl << std::endl;
        // std::cout << "Initialized Laplace-Beltrami" << std::endl;

        Rs.resize(nVertices, Matrix3d::Identity());

        initializeConstraints(filename);
        initializeEquation();
        initializeFirstGuess(should_compute_initial_guess);

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
        // std::cout << "Initialized constraints from outside" << std::endl;
        initializeEquation();
    }

    const std::map<int, Vector3d> getConstraints()
    {
        return constraints;
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
    void initializeConstraints(const std::string &filename)
    {
        constraints.clear();
        change = MatrixXd(1, 3);

        if (filename.find("bar1.off") != std::string::npos)
        {
            /* Constraints for the bar1 */

            change << 75., -90., 0.;

            for (int i = 0; i < nVertices; ++i)
            {
                if (V1(i, 1) > 145)
                {
                    constraints[i] = V0.row(i) + change;
                }

                if (V1(i, 1) < 6)
                {
                    constraints[i] = V0.row(i);
                }
            }
        }
        else if (filename.find("bar2.off") != std::string::npos)
        {
            /* Constraints for the bar2 */

            change << 75., -90., 0.;

            for (int i = 0; i < nVertices; ++i)
            {
                if (V1(i, 1) > 140)
                {
                    constraints[i] = V0.row(i) + change;
                }

                if (V1(i, 1) < 8)
                {
                    constraints[i] = V0.row(i);
                }
            }
        }
        else if (filename.find("square_21.off") != std::string::npos)
        {
            /* Constraints for the square */

            change << 0., -0.5, -0.75;
            // change << 0., -0.75, 0.;
            for (int i = 0; i < nVertices; ++i)
            {
                if (V1(i, 1) > 0.9)
                {
                    constraints[i] = V0.row(i) + change;
                }

                if (V1(i, 1) < 0.1)
                {
                    constraints[i] = V0.row(i);
                }
            }
        }
        else if (filename.find("square_21_spikes.off") != std::string::npos)
        {
            /* Constraints for the square with spikes */

            change << 0.25, -0.75, 0.;
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
        else if (filename.find("cactus_small.off") != std::string::npos)
        {
            /* Constraints for the small cactus */

            change << .5, 0., -.4;

            double maxy = 0;
            int highestVertexId;

            for (int i = 0; i < nVertices; ++i)
            {
                if (V0(i, 2) < .1)
                {
                    constraints[i] = V0.row(i);
                }
                if (V0(i, 2) > maxy)
                {
                    highestVertexId = i;
                    maxy = V0(i, 2);
                }
            }

            constraints[highestVertexId] = V0.row(highestVertexId) + change;
        }
        else if (filename.find("armadillo_1k.off") != std::string::npos)
        {
            /* Constraints for the armadillo */

            int maxy_i, miny_i;
            double maxy = V0.col(1).maxCoeff(&maxy_i);
            double miny = V0.col(1).minCoeff(&miny_i);
            double h = maxy - miny;
            change << 0., h / 3, 0.;

            for (int i = 0; i < nVertices; ++i)
            {
                if (V0(i, 1) <= 0.2 * h + miny)
                {
                    constraints[i] = V0.row(i);
                }
            }

            int maxx_i;
            double maxx = V0.col(0).maxCoeff(&maxx_i);
            constraints[maxx_i] = V0.row(maxx_i) + change;

        }
        else
        {
            /* Constraints for the rest of the shapes */

            int maxy_i, miny_i;
            double maxy = V0.col(1).maxCoeff(&maxy_i);
            double miny = V0.col(1).minCoeff(&miny_i);
            double h = maxy - miny;
            change << h / 2, -h / 2, -h / 2;

            for (int i = 0; i < nVertices; ++i)
            {
                if (V0(i, 1) <= 0.2 * h + miny)
                {
                    constraints[i] = V0.row(i);
                }
            }

            constraints[maxy_i] = V0.row(maxy_i) + change;
        }

        // std::cout << "Initialized constraints" << std::endl;
    }

    void initializeEquation()
    {
        setTheLeftSide();
        sparse_solver.isSymmetric(true);
        sparse_solver.analyzePattern(leftSide);
        sparse_solver.factorize(leftSide);
        rightSide.resize(nVertices, 3);
        // std::cout << "Initialized solver" << std::endl;
    }

    void initializeFirstGuess(bool compute_initial_guess)
    {
        if (compute_initial_guess)
        {
            V1 = laplacianEditingFromMap(V0, F, constraints, LB);
        }
    }

    void setTheLeftSide()
    {
        leftSide = LB;
        for (const auto &[point, new_value] : constraints)
        {
            leftSide.row(point) *= 0;
            leftSide.coeffRef(point, point) = 1;
        }

        // std::cout << "Initialized the left side" << std::endl;
    }

    void computeTheRightSide()
    {
        // std::cout << "Computing the right side..." << std::endl;
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

        // std::cout << "Computed the right side" << std::endl;
    }

    void estimateRotations()
    {
        // std::cout << "Computing the estimated rotations..." << std::endl;
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
        // std::cout << "Computed rotations" << std::endl;
    }

    void estimatePositions()
    {
        V1.setZero();
        computeTheRightSide();
        // std::cout << "Computing the estimated positions..." << std::endl;

        V1.col(0) = sparse_solver.solve(rightSide.col(0));
        V1.col(1) = sparse_solver.solve(rightSide.col(1));
        V1.col(2) = sparse_solver.solve(rightSide.col(2));

        // std::cout << "Computed positions" << std::endl;
    }

    const MatrixXd &V0; // vertex coordinates of the original input mesh
    const MatrixXi &F;
    MatrixXd &V1;  // vertex coordinates of the new subdivided mesh
    int nVertices; // number of vertices, faces

    SparseMatrix<double> LB;       // descrete Laplace-Beltrami operator
    SparseMatrix<double> leftSide; // left side of the equation to estimate the points positions
    MatrixXd rightSide;            // right side of the equation to estimate the points positions
    SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> sparse_solver;

    std::vector<Matrix3d> Rs;

    std::map<int, Vector3d> constraints;
    bool initialized = false;
    bool should_compute_initial_guess;
    MatrixXd change;
};
