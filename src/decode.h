#define MAXROC 32

int evlen, evtype, evnum;
int rocpos[MAXROC], roclen[MAXROC];
int irn[MAXROC];

struct rocinfo {
    int roc;
    unsigned header;
    int index;
};
rocinfo info_HEL;
rocinfo info_RIN;
rocinfo info_HAP;
rocinfo info_TIM;
