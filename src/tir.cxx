#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <libconfig.h>

#include "TROOT.h"

#include "hel.h"

#define M1 0x55555555
#define M2 0x33333333
#define M4 0x0f0f0f0f
#define M8 0x00ff00ff
#define M16 0x0000ffff

FILE *fp1, *fp2;

//Global variables
Int_t gHelicity_rep[LEN];
Int_t gHelicity_act[LEN];
Int_t gQRT[LEN];
Int_t gMPS[LEN];
Int_t gPairSync[LEN];
Int_t gTimeStamp[LEN];
Int_t gSeed_rep[LEN];
Int_t gSeedRing_rep[LEN];
Int_t gEventRing[LEN];
Int_t gError[LEN];
Int_t gN, gNRing;

Int_t readin(Int_t nrun, Bool_t usering);
Int_t predicttir(Int_t nrun, Bool_t usering);
Int_t printout(Int_t nrun);
Int_t RanBit30(Int_t &runseed);
Int_t popcount(Int_t x);
void usage(int argc, char** argv);

Char_t CFGFILE[300] = "./config.cfg";
Char_t INDIR[300] = ".";
Char_t OUTDIR[300] = ".";

int main(int argc, char** argv) {
    int c;
    Bool_t usering = kFALSE;

    while (1) {
        static struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"ring", no_argument, 0, 'r'},
            {"cfgfile", required_argument, 0, 'c'},
            {"indir", required_argument, 0, 'i'},
            {"outdir", required_argument, 0, 'o'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "c:hi:o:r", long_options, &option_index);

        if (c == -1) break;

        switch (c) {
        case 'c':
            strcpy(CFGFILE, optarg);
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'i':
            strcpy(INDIR, optarg);
            break;
        case 'o':
            strcpy(OUTDIR, optarg);
            break;
        case 'r':
            usering = kTRUE;
            break;
        case '?':
            // getopt_long already printed an error message
            break;
        default:
            usage(argc, argv);
        }
    }

    Int_t nrun;

    if (optind < argc) {
        nrun = atoi(argv[optind++]);
    }
    else {
        usage(argc, argv);
        exit(-1);
    }

    config_t cfg;
    config_setting_t *setting;

    config_init(&cfg);

    if (!config_read_file(&cfg, CFGFILE)) {
        fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg), config_error_line(&cfg), config_error_text(&cfg));
        config_destroy(&cfg);
        exit(EXIT_FAILURE);
    }

    Bool_t configerror = kFALSE;

    if (!(config_lookup_int(&cfg, "ndelay", &NDELAY) && config_lookup_int(&cfg, "maxbit", &MAXBIT) && config_lookup_float(&cfg, "windowlength", &WT)))
        configerror = kTRUE;

    if (usering) {
        setting = config_lookup(&cfg, "ringinfo.data");
        if (setting != NULL)
            NRING = config_setting_length(setting);
        else
            configerror = kTRUE;
    }

    if (configerror) {
        fprintf(stderr, "Invalid cfg file\n");
        exit(-1);
    }

    readin(nrun, usering);
    predicttir(nrun, usering);
    printout(nrun);

    return 0;
}

Int_t readin(Int_t nrun, Bool_t usering) {
    printf("Reading TIR helicity information ...\n");

    if ((fp1 = fopen(Form("%s/helTIR_%d.decode.dat", INDIR, nrun), "r")) == NULL) {
        fprintf(stderr, "Can not open %s/helTIR_%d.decode.dat", INDIR, nrun);
        exit(-1);
    }

    Int_t temp[10];

    gN = 0;

    fscanf(fp1, "%d", &temp[0]);
    while (!feof(fp1)) {
        fscanf(fp1, "%d%d%d%d%d%d%d", &gHelicity_rep[gN], &gQRT[gN], &gPairSync[gN], &gMPS[gN], &gTimeStamp[gN], &temp[0], &temp[1]);
        fscanf(fp1, "%d", &temp[0]);
        gN++;
    }

    fclose(fp1);

    if (usering) {
        if ((fp1 = fopen(Form("%s/helRIN_%d.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helRIN_%d.dat", INDIR, nrun);
            exit(-1);
        }

        fscanf(fp1, "%d", &gNRing);

        for (Int_t i = 0; i < gNRing; i++) {
            fscanf(fp1, "%d%d%d%d%x%d", &gEventRing[i], &temp[1], &temp[2], &temp[3], &gSeedRing_rep[i], &temp[4]);
            for (Int_t k = 0; k < NRING; k++)
                fscanf(fp1, "%d", &temp[0]);
        }

        fclose(fp1);
    }

    return 0;
}

Int_t predicttir(Int_t nrun, Bool_t usering) {
    printf("Predicting TIR helicity information ...\n");

    Int_t fPolarityTIR_rep = 0, fPolarityTIR_act = 0;
    Int_t fPhaseTIR_rep = 0;
    Int_t fSeedTIR_rep = 0, fSeedTIR_act = 0;
    Int_t fSeedTIR_fake = 0, fMask = 0;
    Int_t fNSeedTIR = 0;
    Int_t pLastEvent = 0, fLastQRT = 0;
    Int_t fGapEvent = 0, fGapQRT = 0;

    for (Int_t i = 0; i < gN; i++) {
        if (gMPS[i] == 0) {
            // calculate time stamp
            fGapEvent = gTimeStamp[i] - gTimeStamp[pLastEvent];
            fGapQRT = gTimeStamp[i] - fLastQRT;

            // predict helicity
            if (fNSeedTIR == MAXBIT) {
                if (gQRT[i] == 1) {
                    if (gQRT[pLastEvent] == 0) {
                        Int_t MissedQRT = Int_t(fGapEvent / (4 * WT));
                        if (MissedQRT < MAXMISSED) {
                            for (Int_t j = 0; j <= MissedQRT; j++) {
                                fPolarityTIR_rep = RanBit30(fSeedTIR_rep);
                                fPolarityTIR_act = RanBit30(fSeedTIR_act);
                            }
#ifdef DEBUG
                            printf("%d\tD1\t%d\t%d\t%d\t%d\t%d\t%d\n", i + 1, gQRT[i], gHelicity_rep[i], fGapEvent, fGapQRT, MissedQRT, Int_t(fGapEvent / WT));
#endif
                            if (fPolarityTIR_rep == gHelicity_rep[i])
                                fPhaseTIR_rep = 0;
                            else
                                fPhaseTIR_rep = 4;
                        }
                        else {
                            fPhaseTIR_rep = 4;
                        }
                    }
                    else {
                        Int_t j;
                        for (j = 0; j < MAXMISSED; j++) {
                            if (fGapEvent < ((2.0 + j * 4) * WT)) break;
                            fPolarityTIR_rep = RanBit30(fSeedTIR_rep);
                            fPolarityTIR_act = RanBit30(fSeedTIR_act);
                        }
#ifdef DEBUG
                        printf("%d\tD2\t%d\t%d\t%d\t%d\t%d\t%d\n", i + 1, gQRT[i], gHelicity_rep[i], fGapEvent, fGapQRT, j, Int_t(fGapEvent / WT));
#endif
                        if (j < MAXMISSED) {
                            if (fPolarityTIR_rep == gHelicity_rep[i])
                                fPhaseTIR_rep = 0;
                            else
                                fPhaseTIR_rep = 4;
                        }
                        else {
                            fPhaseTIR_rep = 4;
                        }
                    }
                }
                else {
                    Int_t MissedQRT = Int_t(fGapQRT / (4 * WT));
                    Int_t MissedWindow = Int_t(fGapEvent / WT);
#ifdef DEBUG
                    printf("%d\tD3\t%d\t%d\t%d\t%d\t%d\t%d\n", i + 1, gQRT[i], gHelicity_rep[i], fGapEvent, fGapQRT, MissedQRT, MissedWindow);
#endif
                    if ((MissedQRT == 1) && (MissedWindow == 0)
                            && (gQRT[pLastEvent] == 0)) MissedQRT = 0;
                    if (MissedQRT < MAXMISSED) {
                        for (Int_t j = 0; j < MissedQRT; j++) {
                            fPolarityTIR_rep = RanBit30(fSeedTIR_rep);
                            fPolarityTIR_act = RanBit30(fSeedTIR_act);
                        }
                        if (gHelicity_rep[i] == fPolarityTIR_rep)
                            fPhaseTIR_rep = 3;
                        else if (gPairSync[i] == 0)
                            fPhaseTIR_rep = 1;
                        else
                            fPhaseTIR_rep = 2;
                    }
                    else {
                        fPhaseTIR_rep = 4;
                    }
                    fLastQRT = Int_t(fLastQRT + MissedQRT * 4.0 * WT);
                }
                if (fPhaseTIR_rep >= 4) {
                    fNSeedTIR = 0;
                    fSeedTIR_rep = 0;
                }
            }

            //Catch random seed
            if (fNSeedTIR < MAXBIT) {
                if (gQRT[i] == 1) {
                    Int_t k = 0;
                    if (gQRT[pLastEvent] == 0) {
                        if ((fGapEvent < fGapQRT) && (fGapQRT < 6.0 * WT)) {
                            fSeedTIR_rep = ((fSeedTIR_rep << 1 & 0x3FFFFFFF) | gHelicity_rep[i]);
                            fNSeedTIR++;
                        }
                        else {
                            fSeedTIR_rep = 0;
                            fNSeedTIR = 0;
                            if (usering) {
                                for (k = 0; k < MAXMISSED; k++) {
                                    if (fGapQRT < (6.0 + k * 4) * WT) break;
                                    fSeedTIR_fake = (fSeedTIR_fake << 1 & 0x3FFFFFFF) | 0x0;
                                    fMask = (fMask << 1 & 0x3FFFFFFF) | 0x0;
                                }
                            }
                        }
                        if (usering) {
                            fSeedTIR_fake = (fSeedTIR_fake << 1 & 0x3FFFFFFF) | gHelicity_rep[i];
                            fMask = (fMask << 1 & 0x3FFFFFFF) | 0x1;
                        }
                        if (k >= MAXMISSED) fSeedTIR_fake = 0;
                    }
                    else if (gQRT[pLastEvent] == 1) {
                        if (fGapEvent < 2.0 * WT) {
                        }
                        else {
                            if (fGapEvent < 6.0 * WT) {
                                fSeedTIR_rep = ((fSeedTIR_rep << 1 & 0x3FFFFFFF) | gHelicity_rep[i]);
                                fNSeedTIR++;
                            }
                            else {
                                fNSeedTIR = 0;
                                fSeedTIR_rep = 0;
                                if (usering) {
                                    for (k = 0; k < MAXMISSED; k++) {
                                        if (fGapQRT < (2.0 + k * 4) * WT) break;
                                        fSeedTIR_fake = (fSeedTIR_fake << 1 & 0x3FFFFFFF) | 0x0;
                                        fMask = (fMask << 1 & 0x3FFFFFFF) | 0x0;
                                    }
                                }
                            }
                            if (usering) {
                                fSeedTIR_fake = (fSeedTIR_fake << 1 & 0x3FFFFFFF) | gHelicity_rep[i];
                                fMask = (fMask << 1 & 0x3FFFFFFF) | 0x1;
                            }
                            if (k >= MAXMISSED) fSeedTIR_fake = 0;
                        }
                    }
                    if ((usering) && (popcount(fMask) >= 10)) {
                        Int_t j;
                        for (j = 0; j < gNRing; j++) {
                            if (gEventRing[j] >= i + 1) break;
                        }
                        for (Int_t k = j; k < j + 40; k++) {
                            if ((gSeedRing_rep[k] & fMask) == fSeedTIR_fake) {
                                fSeedTIR_rep = gSeedRing_rep[k];
                                fNSeedTIR = MAXBIT;
                                break;
                            }
                        }
                    }
                    if (fNSeedTIR == MAXBIT) {
                        fPolarityTIR_rep = fSeedTIR_rep & 0x01;
                        fSeedTIR_act = fSeedTIR_rep;
                        for (Int_t j = 0; j < NDELAY; j++)
                            fPolarityTIR_act = RanBit30(fSeedTIR_act);
                        gError[i] = 0;
                    }
#ifdef DEBUG
                    printf("%d\tS1\t%d\t%d\t%d\t%d\t%d\t%08x\n", i + 1, gQRT[i], gHelicity_rep[i], fGapEvent, fGapQRT, fNSeedTIR, fMask);
#endif
                }
            }

            // assign actual helicity to event
            if (fNSeedTIR == MAXBIT) {
                if (fPolarityTIR_act == 1) {
                    if (fPhaseTIR_rep == 0 || fPhaseTIR_rep == 3)
                        gHelicity_act[i] = 1;
                    else
                        gHelicity_act[i] = -1;
                }
                else {
                    if (fPhaseTIR_rep == 0 || fPhaseTIR_rep == 3)
                        gHelicity_act[i] = -1;
                    else
                        gHelicity_act[i] = 1;
                }
                gSeed_rep[i] = fSeedTIR_rep;
            }
            else {
                gError[i] |= 0x01;
                gSeed_rep[i] = 0;
                gHelicity_act[i] = 0;
            }
            if (gQRT[i] == 1) fLastQRT = gTimeStamp[i];
            pLastEvent = i;
        }
        else {
            // non-physics trigger event
            gError[i] |= 0x08;
            gSeed_rep[i] = 0;
            gHelicity_act[i] = 0;
        }
    }

    return 0;
}

Int_t printout(Int_t nrun) {
    if ((fp1 = fopen(Form("%s/helTIR_%d.decode.dat", INDIR, nrun), "r")) == NULL) {
        fprintf(stderr, "Can not open %s/helTIR_%d.decode.dat", INDIR, nrun);
        exit(-1);
    }
    if ((fp2 = fopen(Form("%s/helTIR_%d.noring.dat", OUTDIR, nrun), "w")) == NULL) {
        fprintf(stderr, "Can not open %s/helTIR_%d.noring.dat", OUTDIR, nrun);
        exit(-1);
    };

    Int_t temp[10];

    fprintf(fp2, "%d\n", gN);
    for (Int_t i = 0; i < gN; i++) {
        fscanf(fp1, "%d%d%d%d%d%d%d%d", &temp[0], &temp[1], &temp[2], &temp[3], &temp[4], &temp[5], &temp[6], &temp[7]);
        fprintf(fp2, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%08x\t%d\t%d\t%d\n", temp[0], gHelicity_act[i], gHelicity_rep[i], gQRT[i], gPairSync[i], gMPS[i], gTimeStamp[i], gSeed_rep[i], gError[i], temp[6], temp[7]);
    }

    fclose(fp1);
    fclose(fp2);

    return 0;
}

Int_t RanBit30(Int_t &ranseed) {
    // Take 7,28,29,30 bit of ranseed out
    UInt_t bit7 = ((ranseed & 0x00000040) != 0);
    UInt_t bit28 = ((ranseed & 0x08000000) != 0);
    UInt_t bit29 = ((ranseed & 0x10000000) != 0);
    UInt_t bit30 = ((ranseed & 0x20000000) != 0);

    UInt_t newbit = (bit30 ^ bit29 ^ bit28 ^ bit7) & 0x1;

    if (ranseed <= 0) {
        newbit = 0;
    }

    ranseed = ((ranseed << 1) | newbit) & 0x3FFFFFFF;

    return newbit;
}

Int_t popcount(Int_t x) {
    x = (x & M1) + ((x >> 1) & M1);
    x = (x & M2) + ((x >> 2) & M2);
    x = (x & M4) + ((x >> 4) & M4);
    x = (x & M8) + ((x >> 8) & M8);
    x = (x & M16) + ((x >> 16) & M16);
    return x;
}

void usage(int argc, char** argv) {
    printf("usage: %s [options] RUN_NUMBER\n", argv[0]);
    printf("  -c, --cfgfile=config.cfg   Set configuration file name\n");
    printf("  -h, --help                 This small usage guide\n");
    printf("  -i, --indir=.              Set input directory\n");
    printf("  -o, --outdir=.             Set output directory\n");
    printf("  -r, --ring                 Use ring info to help predict\n");
}
