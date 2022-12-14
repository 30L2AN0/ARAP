#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <string>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOFF.h>

#include <set>

#include "ShapeDeformation.cpp"
// #include "LaplacianEditing.cpp"

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V0;
MatrixXd V1;
MatrixXi F;

MeshDeformation deformation(V0, F, V1);
string filename;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
    if (key == '1')
    {
        if (!deformation.isInitialized())
        {
            deformation.initialize(true, filename);
        } else {
            deformation.performOneIteration();
        }

        viewer.data().clear();
        viewer.data().point_size = 10;

        for (const auto &[point, _] : deformation.getConstraints()) {
            viewer.data().add_points(V1.row(point), RowVector3d(1., 0., 0.));
        }

        viewer.data().set_mesh(V1, F);

        return true;
    }

    if (key == 'R' || key == 'r')
    {
        deformation.reset();
        
        viewer.data().clear();
        viewer.data().set_mesh(V1, F);
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
        filename = argv[1];
        std::cout << "reading input file: " << argv[1] << std::endl;
        if (filename.find(".off") != string::npos) {
            igl::readOFF(argv[1], V0, F);
        } else if (filename.find(".obj") != string::npos) {
            igl::readOBJ(argv[1], V0, F);
        }
    }

    V1 = V0;

    igl::opengl::glfw::Viewer viewer; // create the 3d viewer
    std::cout << "Press '1' for one iteration of shape deformation" << std::endl
              << "Press 'S' for save the current mesh to file" << std::endl
              << "Press 'R' to reset the shape to initial state" << std::endl;

    viewer.core().background_color << 1., 1., 1., 1.;

    viewer.callback_key_down = &key_down;
    viewer.data().set_mesh(V0, F);

    viewer.core(0).align_camera_center(V0, F);
    viewer.launch();
}
