#include <iostream>
#include "CImg.h"

using namespace std;
using namespace cimg_library;

int main(int argc, const char* argv[]) {
    cout << "Hello World!" << endl;

    cout << endl;

    // creation d'une image vide
    CImg<> img, visu(512, 512, 24);
    // creation d'un tableau de float vide pour recuperer la taille des voxels
    float voxelSize[3] = {0};

    // lecture de l'image passee en parametre et initialisation de la taille des voxel
    img.load_analyze("knix.hdr", voxelSize);

    // affichage des données
    // EXERCICE 1
    /*
    cout << "\ttaille voxel x : " << voxelSize[0] << endl;
    cout << "\ttaille voxel y : " << voxelSize[1] << endl;
    cout << "\ttaille voxel z : " << voxelSize[2] << endl;

    cout << "\twidth : " << img.width() << endl;
    cout << "\theight : " << img.height() << endl;
    cout << "\tdepth : " << img.depth() << endl;

    cout << "\tvaleur min : " << img.min() << endl;
    cout << "\tvaleur max : " << img.max() << endl;

    cout << "\tvoxel (256, 256, 12) : " << img(256, 256, 12, 0) << endl;

    cout << endl;
    //*/

    // EXERCICE 2-3
    /*

    CImgDisplay disp(512, 512, "");
    bool redraw = true;
    int displayedSlice = 12;
    CImg<> img_tmp = img;

    while(!disp.is_closed() && !disp.is_keyESC()) {
        if(disp.wheel()) {
            if(displayedSlice != 23 && disp.wheel() > 0) {
                displayedSlice++;
                redraw = true;
            }
            if(displayedSlice != 0 && disp.wheel() < 0) {
                displayedSlice--;
                redraw = true;
            }
        }

        if(disp.is_keyL()) {
            img = img_tmp;
            redraw = true;
        }

        if(disp.is_keyM()) {
            img.blur(1, 1, 1, false, false);
            redraw = true;
        }

        //reinitialise l'event wheel (set to 0)
        disp.set_wheel();

        if(redraw) {
            disp.display(img.get_slice(displayedSlice));
            redraw = false;
        }
    }
    //*/

    // EXERCICE 4
    //*
    CImgDisplay dispX(img.width(), img.height(), "")
    CImgDisplay dispY(img.width(), img.depth(), "")
    CImgDisplay dispZ(img.height(), img.depth(), "")

    for(int x = 0; x < img.depth(); x++) {
    for(int x = 0; x < img.depth(); x++) {
    for(int x = 0; x < img.depth(); x++) {
        for(y) {
            for(z) {
                // cherche min z
                // cherche max z
                // calcule moyen z
            }
        }
    }
    }
    }
    //*/

    cout << "Goodbye World!" << endl;

    return 0;
}
