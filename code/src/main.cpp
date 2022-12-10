#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
// #include <igl/doublearea.h>
// #include <igl/massmatrix.h>
// #include <igl/invert_diag.h>
// #include <igl/jet.h>

// #include <igl/gaussian_curvature.h>
// #include <igl/per_vertex_normals.h>
// #include <igl/per_face_normals.h>

// #include "HalfedgeBuilder.cpp"
#include "ShapeDeformation.cpp"

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V0;
MatrixXd V1;
MatrixXi F;

MeshDeformation *deformation;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
	if (key == '1')
	{
		deformation->performOneIteration();
		// deformation->print(1);

		// update the current mesh
		V1 = deformation->getEstimatedVertexCoordinates(); // update vertex coordinates
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

// ------------ main program ----------------
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

	// we never change V0
	deformation = new MeshDeformation(V0, F);

	igl::opengl::glfw::Viewer viewer; // create the 3d viewer
	std::cout << "Press '1' for one iteration of shape deformation" << std::endl
			  << "Press 'S' save the current mesh to file" << std::endl;

	viewer.callback_key_down = &key_down;
	viewer.data().set_mesh(V0, F);

	viewer.core(0).align_camera_center(V0, F);
	viewer.launch();
}
