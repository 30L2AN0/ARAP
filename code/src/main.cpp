#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>

#include "ShapeDeformation.cpp"
// #include "LaplacianEditing.cpp"

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V0;
MatrixXd V1;
MatrixXi F;

MeshDeformation deformation(V0, F, V1);

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
    if (key == '1')
    {
        HalfedgeBuilder *builder = new HalfedgeBuilder();
        HalfedgeDS he = (builder->createMeshWithFaces(V0.rows(), F)); // create the half-edge representation

        std::cout << "Setting constraints for laplacianEditing" << std::endl;
        //----- Constraints for octagon -----
        // Eigen::VectorXi constraintsIds(2);
        // Eigen::MatrixXd constraints(2, 3);
        // constraintsIds << 2, 4;
        // constraints << 1., 0.5, 0.,
        // 			   0., -1., 0.;
        //----- ~Constraints for octagon -----

        //----- Constraints for tower -----
        // const int nbConstraints = 25 * 2 + 16 * 2;
        // Eigen::VectorXi constraintsIds(nbConstraints);
        // Eigen::MatrixXd constraints(nbConstraints, 3);
        // for (size_t id = 0; id < 25; id++)
        // {
        // 	constraintsIds(id) = id;
        // 	constraints.row(id) = V0.row(id);
        // 	constraintsIds(25 + id) = 25 + id;
        // 	constraints.row(25 + id) = V0.row(25 + id) + Eigen::RowVector3d(75, -90, 0);
        // }
        // for (size_t id = 0; id < 16; id++)
        // {
        // 	constraintsIds(50 + id) = 50 + id;
        // 	constraints.row(50 + id) = V0.row(50 + id);
        // 	constraintsIds(50 + 16 + id) = 338 - 16 + id;
        // 	constraints.row(50 + 16 + id) = V0.row(338 - 16 + id) + Eigen::RowVector3d(75, -90, 0);
        // }
        //----- ~Constraints for tower -----

        //----- Constraints for cactus -----
   		const int nbVertices = V0.rows();
		std::vector<int> verticesIds;
		verticesIds.reserve(10);
		int nbLowVertices = 0;
		int highestVertexId;
		for (size_t vertexId = 0; vertexId < nbVertices; vertexId++)
		{
			if (V0(vertexId, 2) < .1)
			{
				nbLowVertices++;
				verticesIds.push_back(vertexId);
			}
			if (V0(vertexId, 2) > .89)
			{
				highestVertexId = vertexId;
			}
		}
		const int nbConstraints = nbLowVertices + 1;
		Eigen::VectorXi constraintsIds(nbConstraints);
		Eigen::MatrixXd constraints(nbConstraints, 3);
		for (size_t id = 0; id < nbLowVertices; id++)
		{
			constraintsIds(id) = verticesIds[id];
			constraints.row(id) = V0.row(verticesIds[id]);
		}
		constraintsIds(nbLowVertices) = highestVertexId;
		constraints.row(nbLowVertices) = V0.row(highestVertexId) + Eigen::RowVector3d(.5, 0, -.4);
        //----- ~Constraints for cactus -----
        
        std::cout << "Call laplacianEditing" << std::endl;
        V1 = laplacianEditing(V0, he, constraintsIds, constraints);

        if (!deformation.isInitialized()) {
            deformation.initialize(false);
        }
        deformation.setConstraints(constraintsIds, constraints);

        viewer.data().clear();
        viewer.data().set_mesh(V1, F);
        return true;
    }
    if (key == '2')
    {
        if (!deformation.isInitialized()) {
            deformation.initialize(true);
        } else {
            deformation.performOneIteration();
        }

        viewer.data().clear();
        viewer.data().set_mesh(V1, F);
        // viewer.append_mesh();
        // viewer.data().set_mesh(V0, F);
        // viewer.data(0).set_colors(Eigen::RowVector3d(0.3, 0.8, 0.3));
        // viewer.data(1).set_colors(Eigen::RowVector3d(0.8, 0.3, 0.3));

        return true;
    }

    if (key == 'S' || key == 's') // write the mesh to file (OFF format)
    {
        igl::writeOFF("../data/output.off", V1, F);
        return true;
    }

    return false;
}

/**
 * Create a triangle mesh corresponding to an octagon inscribed in the unit circle
 */
void createOctagon(MatrixXd &Vertices, MatrixXi &Faces)
{
    Vertices = MatrixXd(6, 3);
    Faces = MatrixXi(8, 3);

    Vertices << 0.0, 0.0, 1.0,
        1.000000, 0.000000, 0.000000,
        0.000000, 1.000000, 0.000000,
        -1.000000, 0.000000, 0.000000,
        0.000000, -1.000000, 0.000000,
        0.000000, 0.000000, -1.000000;

    Faces << 0, 1, 2,
        0, 2, 3,
        0, 3, 4,
        0, 4, 1,
        5, 2, 1,
        5, 3, 2,
        5, 4, 3,
        5, 1, 4;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Creating an octagon" << std::endl;
        createOctagon(V0, F);
    }
    else
    {
        std::cout << "reading input file: " << argv[1] << std::endl;
        igl::readOFF(argv[1], V0, F);
    }

    V1 = V0;

    igl::opengl::glfw::Viewer viewer; // create the 3d viewer
    std::cout << "Press '1' for one iteration of shape deformation" << std::endl
              << "Press 'S' save the current mesh to file" << std::endl;

    viewer.callback_key_down = &key_down;
    viewer.data().set_mesh(V0, F);

    viewer.core(0).align_camera_center(V0, F);
    viewer.launch();
}
