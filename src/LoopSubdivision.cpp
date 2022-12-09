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

#ifndef HALFEDGE_DS_HEADER
  #define HALFEDGE_DS_HEADER
  #include "HalfedgeDS.cpp"
#endif

using namespace Eigen;
using namespace std;

/**
 * @author Luca Castelli Aleardi (2019)
 */
class LoopSubdivision
{

public:
    /** 
	 * Initialize the data structures
	 **/
    LoopSubdivision(MatrixXd &V_original, MatrixXi &F_original, HalfedgeDS &mesh)
    {
        he = &mesh;
        V = &V_original;
		F = &F_original; // NOT NEEDED if using the half-edge data structure
        int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
        int n = V_original.rows();         // number of vertices in the original mesh
		int f = F_original.rows();

		// TO BE COMPLETED (initialize arrays V1 and F1)
		V1 = new MatrixXd(n + e, 3);
		F1 = new MatrixXi(f * 4, 3);

        he->print();
	}

	~LoopSubdivision() {
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
        int n = he->sizeOfVertices();      // number of vertices in the original mesh
        int f = he->sizeOfFaces();         // number of vertices in the original mesh

        std::cout << "Computing new points" << std::endl;
		int v_i = 0;
		for (int i = 0; i < he->sizeOfHalfedges(); ++i) {
			if (he->getOpposite(i) > i) {
				V1->row(v_i) = computeEdgePoint(i);
				edge_indexes[i] = v_i;
				edge_indexes[he->getOpposite(i)] = v_i++;
			}
		}

        // нужно апдейтить с учетом новых соседей
        std::cout << "Updating old points" << std::endl;
        for (int i = 0; i < V->rows(); ++i) {
            V1->row(i + e) = updatedOriginalPoint(i);
        }

        std::cout << "Filling faces" << std::endl;
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
        int op = he->getOpposite(h);

        int v0 = he->getTarget(h);
		int v1 = he->getTarget(op);
        int v2 = he->getTarget(he->getNext(h));
        int v3 = he->getTarget(he->getNext(op));

        if (he->getFace(op) != -1) {
            // std::cout << "new point on edge " << h << ": " << 3. / 8. * (V->row(v0) + V->row(v1)) + 1. / 8. * (V->row(v2) + V->row(v3)) << std::endl;
            return 3. / 8. * (V->row(v0) + V->row(v1)) + 1. / 8. * (V->row(v2) + V->row(v3));
        } else {
            return 1. / 2. * (V->row(v0) + V->row(v1));
        }
    }

    /**
	 * Given a vertex 'v' of the original mesh, compute and return its new coordinates
	 */
    MatrixXd updatedOriginalPoint(int v)
    {
        int edge = he->getEdge(v);
        int edge_op = he->getOpposite(edge);
        int edge_next = he->getNext(edge);

        // std::cout << "0 neighbor: " << edge << " points to " << V->row(he->getTarget(edge_op)) << std::endl;

        MatrixXd boundary_sum = 3. / 4. * V->row(v);
        bool boundary = false;
        if (he->getFace(edge) == -1 || he->getFace(edge_op) == -1) {
            boundary = true;
            boundary_sum += 1. / 8. * V->row(he->getTarget(edge_op));
        }

        double d = 1;
        MatrixXd nsum = V->row(he->getTarget(edge_op));
        edge_op = he->getOpposite(edge_next);
        
        while (edge_op != edge) {
            // std::cout << d << " neighbor: " << edge_next << " points to " << V->row(he->getTarget(edge_next)) << std::endl;
            nsum += V->row(he->getTarget(edge_next));

            if (he->getFace(edge_next) == -1 || he->getFace(edge_op) == -1) {
                boundary = true;
                boundary_sum += 1. / 8. * V->row(he->getTarget(edge_next));
            }

            edge_next = he->getNext(edge_op);
            edge_op = he->getOpposite(edge_next);
            d++;
        }

        double alpha = 3. / 16.;
        if (d > 3) {
            alpha = 3. / (8. * d);
        }

        // std::cout << "updating vertex " << v << "; degree " << d << "; alpha " << alpha << "; boundary: " << boundary << std::endl;

        nsum *= alpha;
        nsum += (1. - alpha * d) * V->row(v);

        // std::cout << "nsum: " << nsum << std::endl;

        if (boundary) {
            return boundary_sum;
        }

        return nsum;
    }

    /** Half-edge representation of the original input mesh */
    HalfedgeDS *he;
    MatrixXd *V; // vertex coordinates of the original input mesh

	/** faces/vertex incidence relations in the original mesh */
	MatrixXi *F; // REMARK: not needed if using the half-edge data structure

    int nVertices, nFaces; // number of vertices, faces in the new subdivided mesh
    MatrixXd *V1;          // vertex coordinates of the new subdivided mesh
    MatrixXi *F1;          // faces of the new subdivided mesh

    map<int, int> edge_indexes;
};
