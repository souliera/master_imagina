// g++ -o a.out p2.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11

#include <iostream>
#include "CImg.h"

using namespace std;
using namespace cimg_library;


int main(int argc, const char* argv[]) {
    if(argc != 3) {
        cout << "Usage: " << argv[0] << "<image.hdr> <seuil>" << endl;
        exit(EXIT_FAILURE);
    }

    // initialisation
    CImg<> img, img_bin;
    float voxelSize[3] = {0};

    const unsigned char white[] = {255, 255, 255};

    img.load_analyze(argv[1], voxelSize);
    unsigned int width = img.width();
    unsigned int height = img.height();
    unsigned int depth = img.depth();
    img.resize(width, height, depth, 3);

    img_bin.load_analyze(argv[1], voxelSize);
    CImg<> img_tmp(width, height, depth);
    CImgDisplay disp(width, height, "");

    // seuil binaire
    unsigned int initThreshold = atoi(argv[2]);
    unsigned int threshold = atoi(argv[2]);
    img_bin.threshold(threshold);
    img_bin.save_analyze("bin.hdr", voxelSize);

    bool redraw = true;
    int displayedSlice = depth / 2;
    int x, y, z;
    int sum = 0;

    FILE* file = fopen("profil.dat", "w");

    while(!disp.is_closed() && !disp.is_keyESC()) {
        // img.display(disp);
        // disp.display(img);
        disp.show().flush();
        img.display(disp);
        // if(disp.wheel()) {
        //     if(displayedSlice != depth-1 && disp.wheel() > 0) {
        //         displayedSlice++;
        //         redraw = true;
        //     }
        //     if(displayedSlice != 0 && disp.wheel() < 0) {
        //         displayedSlice--;
        //         redraw = true;
        //     }
        //
        //     //reinitialise l'event wheel (set to 0)
        //     disp.set_wheel();
        // }
        //
        // if(disp.button()&1) {
        //     x = disp.mouse_x();
        //     y = disp.mouse_y();
        //     z = displayedSlice;
        //     // for(int l = 0; l < 100; l++) {
        //     //     sum = 0;
        //     //     img_bin.load_analyze(argv[1], voxelSize);
        //     //     img_bin.threshold(threshold);
        //         img_bin.draw_fill(x, y, z, white, 1.0f, img_tmp, 0.0f);
        //     //     for(int i = 0; i < width; i++) {
        //     //         for(int j = 0; j < height; j++) {
        //     //             for(int k = 0; k < depth; k++) {
        //     //                 if(img_tmp(i, j, k) > 0.0f) {
        //     //                     sum++;
        //     //                 }
        //     //             }
        //     //         }
        //     //     }
        //     //     fprintf(file, "%i, %i \n", initThreshold - l, sum);
        //     //     threshold--;
        //     //
        //     // }
        //     img_tmp.save_analyze("mask.hdr", voxelSize);
        //     break;
        // }
        //
        // if(redraw) {
        //     disp.display(img_bin.get_slice(displayedSlice));
        //     redraw = false;
        // }
    }

    fclose(file);

    exit(EXIT_SUCCESS);
}
