#include <iostream>
#include <cstdlib>
#include <cmath>

#include "image_ppm.h"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: <ImageIn.pgm>" << endl;
        exit(EXIT_FAILURE);
    }

    cout << "Hello World" << endl;

    // ************
    // * Variable *
    // ************

    int height = 0;
    int width = 0;
    int size = 0;

    OCTET *imgIn;

    int histo[256] = {0};
    float sum = 0.0f;

    // ******************
    // * Initialisation *
    // ******************

    lire_nb_lignes_colonnes_image_pgm(argv[1], &height, &width);
    size = height * width;

    allocation_tableau(imgIn, OCTET, size);
    lire_image_pgm(argv[1], imgIn, height * width);

    // *************
    // * Programme *
    // *************

    // histogram
    for(int i = 0; i < size; i++) {
        histo[imgIn[i]]++;
    }

    // sum
    for(int i = 0; i < 256; i++) {
        if((float)histo[i] / (float)size != 0) {
            sum += ((float)histo[i] / (float)size) * (log2((float)histo[i] / (float)size));
        }
    }

    float entropie = -sum;

    cout << "Entropie ordre zéro de " << argv[1] << " : " << entropie << " bpp" << endl;

    // ********
    // * Save *
    // ********

    free(imgIn);

    cout << "Goodbye World" << endl;

    return 0;
}
