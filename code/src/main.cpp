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

MatrixXd V;
// MatrixXd V1;
MatrixXi F;

// MeshDeformation *deformation;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
	if (key == '1')
	{
		HalfedgeBuilder *builder = new HalfedgeBuilder();
		HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F)); // create the half-edge representation

		// Constraints for tower
		// const int nbConstraints = 25 * 2 + 16 * 2;
		// Eigen::VectorXi constraintsIds(nbConstraints);
		// Eigen::MatrixXd constraints(nbConstraints, 3);
		// for (size_t id = 0; id < 25; id++)
		// {
		// 	constraintsIds(id) = id;
		// 	constraints.row(id) = V.row(id);
		// 	constraintsIds(25 + id) = 25 + id;
		// 	constraints.row(25 + id) = V.row(25 + id) + Eigen::RowVector3d(75, -90, 0);
		// }
		// for (size_t id = 0; id < 16; id++)
		// {
		// 	constraintsIds(50 + id) = 50 + id;
		// 	constraints.row(50 + id) = V.row(50 + id);
		// 	constraintsIds(50 + 16 + id) = 338 - 16 + id;
		// 	constraints.row(50 + 16 + id) = V.row(338 - 16 + id) + Eigen::RowVector3d(75, -90, 0);
		// }
		// ~Constraints for tower

		// Constraints for cactus
		const int nbVertices = V.rows();
		std::vector<int> verticesIds;
		verticesIds.reserve(10);
		int nbLowVertices = 0;
		int highestVertexId;
		for (size_t vertexId = 0; vertexId < nbVertices; vertexId++)
		{
			if (V(vertexId, 2) < .1)
			{
				nbLowVertices++;
				verticesIds.push_back(vertexId);
			}
			if (V(vertexId, 2) > .89)
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
			constraints.row(id) = V.row(verticesIds[id]);
		}
		constraintsIds(nbLowVertices) = highestVertexId;
		constraints.row(nbLowVertices) = V.row(highestVertexId) + Eigen::RowVector3d(.5, 0, -.4);
		// ~Constraints for cactus

		std::cout << "faces: " << std::endl
				  << F << std::endl
				  << std::endl;
		Eigen::MatrixXd V1 = laplacianEditing(V, he, constraintsIds, constraints);
		std::cout << "constraintsIds: " << std::endl
				  << constraintsIds << std::endl
				  << std::endl;
		std::cout << "constraints: " << std::endl
				  << constraints << std::endl
				  << std::endl;
		std::cout << "nb vertices: " << std::endl
				  << V.rows() << std::endl
				  << std::endl;
		std::cout << "V before: " << std::endl
				  << V << std::endl
				  << std::endl;
		std::cout << "V after: " << std::endl
				  << V1 << std::endl
				  << std::endl;
		// std::cout << "highest vertex: " << std::endl
		// 		  << V.row(highestVertexId) << std::endl
		// 		  << std::endl;
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

	// we never change V0
	// deformation = new MeshDeformation(V0, F);

	igl::opengl::glfw::Viewer viewer; // create the 3d viewer
	std::cout << "Press '1' for one iteration of shape deformation" << std::endl
			  << "Press 'S' save the current mesh to file" << std::endl;

	viewer.callback_key_down = &key_down;
	viewer.data().set_mesh(V, F);

	viewer.core(0).align_camera_center(V, F);
	viewer.launch();
}
