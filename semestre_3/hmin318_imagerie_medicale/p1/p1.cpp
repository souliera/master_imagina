// g++ -o a.out p1.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11

#include <iostream>
#include "CImg.h"

using namespace std;
using namespace cimg_library;

// 256 124 256
// 333 320 256

int main(int argc, const char* argv[]) {
    if(argc != 5) {
        cout << "Usage: " << argv[0] << " threshold erodeIteration k-connectivity dilateIteration" << endl;
        return 1;
    }

    // initialisation
    CImg<> img;
    float voxelSize[3] = {0};

    img.load_analyze("accuracy_adult01.hdr", voxelSize);
    unsigned int width = img.width();
    unsigned int height = img.height();
    unsigned int depth = img.depth();

    // seuil binaire
    unsigned int threshold = atoi(argv[1]);
    img.threshold(threshold);

    // p erosion
    unsigned int ep = atoi(argv[2]);
    for(unsigned int i = 0; i < ep; i++) {
        img.erode(2, 2, 2);
    }

    // label
    if(atoi(argv[3]) == 4) {
        img.label(false, 0);
    } else if(atoi(argv[3]) == 8) {
        img.label(true, 0);
    } else {
        cout << "k-connectivity = 4 or 8" << endl;
        return 1;
    }

    // dilatation
    unsigned int dp = atoi(argv[4]);
    for(unsigned int i = 0; i < dp; i++) {
        img.dilate(2, 2, 2);
    }

    // save
    img.save_analyze("accuracy_adult01_result.hdr", voxelSize);

    return 0;
}

