/* ========================================================================= *
 *                                                                           *
 *                       Luca Castelli Aleardi                       		 *
 *           Copyright (c) 2019, Ecole Polytechnique                		 *
 *           Department of Computer Science                  				 *
 *                          All rights reserved.                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of the course material developed for		             *
 *   INF574 Digital Representation and Analysis of Shapes (2019/20)			 *
 * ========================================================================= */
#include <igl/opengl/glfw/Viewer.h>

#include <map>

#ifndef HALFEDGE_DS_HEADER
  #define HALFEDGE_DS_HEADER
  #include "HalfedgeDS.cpp"
#endif

using namespace Eigen;
using namespace std;

/**
 * @author Luca Castelli Aleardi (2019)
 */
class SphereGeneration
{

public:
	/** 
	 * Initialize the data structures
	 **/
	SphereGeneration(MatrixXd &V_original, MatrixXi &F_original, HalfedgeDS &mesh)
	{
		he = &mesh;
		V = &V_original;
		// F = &F_original; // NOT NEEDED if using the half-edge data structure
		int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
		int n = V_original.rows();		   // number of vertices in the original mesh
		int f = F_original.rows();

		// TO BE COMPLETED (initialize arrays V1 and F1)
		V1 = new MatrixXd(n + e, 3);
		F1 = new MatrixXi(f * 4, 3);
	}

	~SphereGeneration() {
		delete V1;
		delete F1;
	}

	/**
	 * Perform the subdivision of the mesh (just perform one round subdivision). <b>
	 * As result, the resulting subdivided mesh is stored in arrays 'V1' and 'F1' 
	 **/
	void subdivide()
	{
		std::cout << "Performing one round subdivision" << endl;
		int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
		int n = he->sizeOfVertices();	   // number of vertices in the original mesh
		int f = he->sizeOfFaces();		   // number of faces in the original mesh

		V1->bottomRows(n) = *V;
		V1->rowwise().normalize();

		int v_i = 0;
		for (int i = 0; i < he->sizeOfHalfedges(); ++i) {
			if (he->getOpposite(i) > i) {
				V1->row(v_i) = computeEdgePoint(i);
				edge_indexes[i] = v_i;
				edge_indexes[he->getOpposite(i)] = v_i++;
			}
		}

		int f_i = 0;
		for (int i = 0; i < f; ++i) {
			int e0 = he->getEdgeInFace(i);
			int e1 = he->getNext(e0);
			int e2 = he->getNext(e1);

			int m0 = edge_indexes[e0];
			int m1 = edge_indexes[e1];
			int m2 = edge_indexes[e2];

			int v0 = he->getTarget(e0) + e;
			int v1 = he->getTarget(e1) + e;
			int v2 = he->getTarget(e2) + e;
			
			F1->row(f_i++) << v0, m1, m0;
			F1->row(f_i++) << m1, m2, m0;
			F1->row(f_i++) << m1, v1, m2;
			F1->row(f_i++) << m0, m2, v2;
		}
	}

	/**
	 * Return the number of half-edges
	 **/
	MatrixXd getVertexCoordinates()
	{
		return *V1;
	}

	/** 
	 * Return the number of faces
	 **/
	MatrixXi getFaces()
	{
		return *F1;
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
				cout << "v" << i << ": " << V1->row(i) << endl;
			}

			std::cout << "new faces: " << nFaces << endl;
			for (int i = 0; i < nFaces; i++)
			{
				cout << "f" << i << ": " << F1->row(i) << endl;
			}
		}
	}

private:
	/**
	 * Compute the midpoint of the given half-edge 'h=(u,v)'
	 */
	MatrixXd computeEdgePoint(int h)
	{
		int v = he->getTarget(h);
		int u = he->getTarget(he->getOpposite(h));
		MatrixXd point = (V->row(v) + V->row(u)) / 2;
		return point.normalized();
	}

	/** Half-edge representation of the original input mesh */
	HalfedgeDS *he;
	MatrixXd *V; // vertex coordinates of the original input mesh

	/** faces/vertex incidence relations in the original mesh */
	// MatrixXi *F; // REMARK: not needed if using the half-edge data structure

	int nVertices, nFaces; // number of vertices, faces in the new subdivided mesh
	MatrixXd *V1;		   // vertex coordinates of the new subdivided mesh
	MatrixXi *F1;		   // faces of the new subdivided mesh

	map<int, int> edge_indexes;
};
