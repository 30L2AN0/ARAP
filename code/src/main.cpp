#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>

#include <igl/gaussian_curvature.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

#include "HalfedgeBuilder.cpp"
#include "SphereGeneration.cpp"
#include "LoopSubdivision.cpp"

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V;
MatrixXi F;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
	if (key == '1')
	{
		HalfedgeBuilder *builder = new HalfedgeBuilder();
		HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F)); // create the half-edge representation
		SphereGeneration *generator = new SphereGeneration(V, F, he);   //
		generator->subdivide();										 // perform one round subdivision
		generator->print(0);

		// update the current mesh
		V = generator->getVertexCoordinates(); // update vertex coordinates
		F = generator->getFaces();
		viewer.data().clear();
		viewer.data().set_mesh(V, F);
		return true;
	}
	if (key == '2')
	{
		HalfedgeBuilder *builder = new HalfedgeBuilder();
		HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F)); // create the half-edge representation
		LoopSubdivision *loop = new LoopSubdivision(V, F, he);   //
		loop->subdivide();										 // perform one round subdivision
		loop->print(0);

		// update the current mesh
		V = loop->getVertexCoordinates(); // update vertex coordinates
		F = loop->getFaces();
		viewer.data().clear();
		viewer.data().set_mesh(V, F);
		return true;
	}

	if (key == 'S' || key == 's') // write the mesh to file (OFF format)
	{
		igl::writeOFF("../data/output.off", V, F);
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
		createOctagon(V, F);
	}
	else
	{
		std::cout << "reading input file: " << argv[1] << std::endl;
		igl::readOFF(argv[1], V, F);
	}

	igl::opengl::glfw::Viewer viewer; // create the 3d viewer
	std::cout << "Press '1' for one round sphere generation" << std::endl
			  << "Press '2' for one round Loop subdivision" << std::endl
			  << "Press 'S' save the current mesh to file" << std::endl;

	viewer.callback_key_down = &key_down;
	viewer.data().set_mesh(V, F);

	viewer.core(0).align_camera_center(V, F);
	viewer.launch(); // run the viewer
}
