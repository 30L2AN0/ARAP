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

#include "HalfedgeBuilder.cpp"
// #include "SphereGeneration.cpp"
// #include "LoopSubdivision.cpp"
#include "Transformation.cpp"

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
		HalfedgeBuilder *builder = new HalfedgeBuilder();
		HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F)); // create the half-edge representation
		Eigen::VectorXi constraintsIds(2);
		constraintsIds << 0, 11;
		Eigen::MatrixXd constraints(constraintsIds.rows(), 3);
		constraints << V.row(0) + Eigen::RowVector3d(0, 0, 1), V.row(11) + Eigen::RowVector3d(0, 1, -1);
		std::cout << "faces: " << std::endl << F << std::endl << std::endl;
		Eigen::MatrixXd V1 = laplacianEditing(V, he, constraintsIds, constraints);
		std::cout << "constraintsIds: " << std::endl << constraintsIds << std::endl << std::endl;
		std::cout << "constraints: " << std::endl << constraints << std::endl << std::endl;
		std::cout << "nb vertices: " << std::endl << V.rows() << std::endl << std::endl;
		std::cout << "V before: " << std::endl << V << std::endl << std::endl;
		std::cout << "V after: " << std::endl << V1 << std::endl << std::endl;
		viewer.data().clear();
		viewer.data().set_mesh(V1, F);
		return true;
	}
	if (key == '2')
	{
		// HalfedgeBuilder *builder = new HalfedgeBuilder();
		// HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F)); // create the half-edge representation
		// LoopSubdivision *loop = new LoopSubdivision(V, F, he);		 //
		// loop->subdivide();											 // perform one round subdivision
		// loop->print(0);

		// // update the current mesh
		// V = loop->getVertexCoordinates(); // update vertex coordinates
		// F = loop->getFaces();
		// viewer.data().clear();
		// viewer.data().set_mesh(V, F);
		// return true;
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
