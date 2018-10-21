#include <iostream>
#include <cmath>
#include <typeinfo>

#include "image_ppm.h"

using namespace std;

bool is_prime(unsigned int n) {
    for(unsigned int i = 2; i < sqrt(n); i++) {
        if(n % i == 0) {
            return false;
        }
    }

    return true;
}

bool is_prime_between_them(unsigned int p, unsigned int q) {
    if(q ? is_prime_between_them(q, p%q) : p == 1) {
        return true;
    } else {
        return false;
    }
}

void list_exposant(unsigned int phi) {
    for(unsigned int i = 2; i < phi; i++) {
        if(is_prime_between_them(i, phi)) {
            cout << i << " ";
        }
    }
}

int modpow(int pixel, unsigned int e, unsigned int n) {
    int c = 1;

    while(e > 0) {
        c = (pixel * c) % n;
        e--;
    }

    return c;
}

unsigned int invers_modulo(unsigned int e, unsigned int n) {
    for(int i = 1; i < n; i++) {
        if((e * i) % n == 1) {
            return i;
        }
    }
}

int main(int argc, char* argv[]) {
    if(argc != 2) {
        cout << "Usage: " << argv[0] << " <image>" << endl;
        exit(EXIT_FAILURE);
    }

    cout << "Hello World!" << endl;

    unsigned int p = 11;
    unsigned int q = 23;
    unsigned int n = p*q;
    unsigned int phi = (p - 1) * (q - 1);
    unsigned int e = 17;

    OCTET *imgIn, *imgOut;
    int height, width, size;

    // cout << p << " : " << is_prime(p) << endl;
    // cout << q << " : " << is_prime(q) << endl;
    // cout << "10 : " << is_prime(10) << endl;

    // cout << p << " et " << q << " : " << is_prime_between_them(p, q) << endl;
    // cout << "2 et 5 " << is_prime_between_them(2, 5) << endl;
    // cout << "10 et 20 " << is_prime_between_them(10, 20) << endl;

    // cout << "Liste exposant :" << endl;
    // list_exposant((p-1)*(q-1));
    // cout << endl;

    lire_nb_lignes_colonnes_image_pgm(argv[1], &height, &width);
    size = height * width;
    allocation_tableau(imgIn, OCTET, size);
    allocation_tableau(imgOut, OCTET, size);
    lire_image_pgm(argv[1], imgIn, size);

    // threshold
    for(int i = 0; i < size; i++) {
        if(imgIn[i] < 128) {
            imgIn[i] = 0;
        } else {
            imgIn[i] = 172;
        }
    }

    // chiffrement
    // public key (e, p*q)
    for(int i = 0; i < size; i++) {
        imgOut[i] = modpow(imgIn[i], e, n);
    }

    ecrire_image_pgm((char*)"rsc/chiff_seuil.pgm", imgOut, height, width);

    // dechiffrement
    // private key (d)
    unsigned int d = invers_modulo(e, phi);

    for(int i = 0; i < size; i++) {
        imgIn[i] = modpow(imgOut[i], d, n);
    }

    ecrire_image_pgm((char*)"rsc/dechiff_seuil.pgm", imgIn, height, width);

    free(imgIn);
    free(imgOut);

    cout << "Goodbye World!" << endl;

    exit(EXIT_SUCCESS);
}
