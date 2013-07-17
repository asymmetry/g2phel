#include "TROOT.h"

Bool_t isexist(Char_t* fname) {
    FILE *temp;
    Bool_t isopen;

    if ((temp = fopen(fname, "r")) == NULL) {
        isopen = false;
    }
    else {
        isopen = true;
        fclose(temp);
    }

    return isopen;
}
