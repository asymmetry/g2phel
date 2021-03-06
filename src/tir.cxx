#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <getopt.h>

#include <libconfig.h>

#include "TROOT.h"

#define M1 0x55555555
#define M2 0x33333333
#define M4 0x0f0f0f0f
#define M8 0x00ff00ff
#define M16 0x0000ffff

#define MAXMISSED 300

FILE *fp1, *fp2;

// global variables
Int_t NRING = 0;
Int_t NHAPPEX = 0;
Int_t NDELAY;
Int_t MAXBIT;
Double_t WT, MPST;
Bool_t USEHAPPEX;

Int_t *gHelicity_act;
Int_t *gHelicity_rep;
Int_t *gPairSync;
Int_t *gQRT;
Int_t *gMPS;
Int_t *gTimeStamp;
Int_t *gSeed_rep;
Int_t *gError;
Int_t gN;

Int_t *gSeedRing_rep;
Int_t *gEventRing;
Int_t gNRing;

Int_t readin(Int_t nrun, Bool_t usering);
Int_t predicttir(Bool_t usering);
Int_t printout(Int_t nrun, Bool_t usering);
Int_t RanBit30(Int_t &runseed);
Int_t popcount(Int_t x);
void usage(Int_t argc, Char_t **argv);

Char_t CFGFILE[300] = "./config.cfg";
Char_t INDIR[300] = ".";
Char_t OUTDIR[300] = ".";

Int_t main(Int_t argc, Char_t **argv)
{
    Int_t c;
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

        Int_t option_index = 0;
        c = getopt_long(argc, argv, "c:hi:o:r", long_options, &option_index);

        if (c == -1)
            break;

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

    if (optind < argc)
        nrun = atoi(argv[optind++]);
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

    if (!(config_lookup_int(&cfg, "ndelay", &NDELAY) && config_lookup_int(&cfg, "maxbit", &MAXBIT) && config_lookup_float(&cfg, "windowlength", &WT) && config_lookup_float(&cfg, "mpslength", &MPST)))
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

    clock_t start = clock();
    readin(nrun, usering);
    predicttir(usering);
    printout(nrun, usering);
    clock_t end = clock();

    printf("Actual helicity calculated in %5.3f s\n", (Double_t)(end - start) / (Double_t) CLOCKS_PER_SEC);

    return 0;
}

Int_t readin(Int_t nrun, Bool_t usering)
{
    printf("Reading TIR helicity information ...\n");

    if ((fp1 = fopen(Form("%s/helTIR_%d.decode.dat", INDIR, nrun), "r")) == NULL) {
        fprintf(stderr, "Can not open %s/helTIR_%d.decode.dat\n", INDIR, nrun);
        exit(-1);
    }

    Int_t tempi[10];
    Char_t tempc[300];

    gN = 0;

    fscanf(fp1, "%d", &tempi[0]);

    while (!feof(fp1)) {
        fgets(tempc, 300, fp1);
        gN++;
        fscanf(fp1, "%d", &tempi[0]);
    }

    gHelicity_rep = new Int_t[gN];
    gHelicity_act = new Int_t[gN];
    gQRT = new Int_t[gN];
    gPairSync = new Int_t[gN];
    gMPS = new Int_t[gN];
    gTimeStamp = new Int_t[gN];
    gSeed_rep = new Int_t[gN];
    gError = new Int_t[gN];

    memset(gHelicity_rep, 0, sizeof(Int_t) * gN);
    memset(gHelicity_act, 0, sizeof(Int_t) * gN);
    memset(gQRT, 0, sizeof(Int_t) * gN);
    memset(gPairSync, 0, sizeof(Int_t) * gN);
    memset(gMPS, 0, sizeof(Int_t) * gN);
    memset(gTimeStamp, 0, sizeof(Int_t) * gN);
    memset(gSeed_rep, 0, sizeof(Int_t) * gN);
    memset(gError, 0, sizeof(Int_t) * gN);

    rewind(fp1);

    for (Int_t i = 0; i < gN; i++)
        fscanf(fp1, "%d%d%d%d%d%d%d%d", &tempi[0], &gHelicity_rep[i], &gQRT[i], &gPairSync[i], &gMPS[i], &gTimeStamp[i], &tempi[6], &tempi[7]);

    fclose(fp1);

    if (usering) {
        if ((fp1 = fopen(Form("%s/helRIN_%d.nalign.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helRIN_%d.dat\n", INDIR, nrun);
            exit(-1);
        }

        fscanf(fp1, "%d", &gNRing);

        gEventRing = new Int_t[gNRing];
        gSeedRing_rep = new Int_t[gNRing];

        memset(gEventRing, 0, sizeof(Int_t) * gNRing);
        memset(gSeedRing_rep, 0, sizeof(Int_t) * gNRing);

        for (Int_t i = 0; i < gNRing; i++) {
            fscanf(fp1, "%d%d%d%d%x%d", &gEventRing[i], &tempi[1], &tempi[2], &tempi[3], &gSeedRing_rep[i], &tempi[5]);

            for (Int_t k = 0; k < NRING; k++)
                fscanf(fp1, "%d", &tempi[0]);
        }

        fclose(fp1);
    }

    return 0;
}

Int_t predicttir(Bool_t usering)
{
    printf("Calculating actual helicity information ...\n");

    Int_t fSeedTIR_rep = 0, fPolarityTIR_rep = 0, fPhaseTIR_rep = 0;
    Int_t fSeedTIR_act = 0, fPolarityTIR_act = 0;
    Int_t fSeedTIR_fake = 0, fMask = 0;
    Int_t fNSeedTIR = 0;
    Int_t pLastEvent = 0, pLastMPS = 0;
    Int_t fPhaseLastTIR = -1;
    Int_t fTimeLastQRT = 0;
    Double_t fTimeGap = 0, fTimeGapQRT = 0;
    Double_t fTimeRef = -1, fTimeRefQRT = -1; // Use MPS==1 event to set a reference

    for (Int_t i = 0; i < gN; i++) {
        if (gMPS[i] == 0) {
            // calculate time stamp
            if (fTimeRef > 0) {
                fTimeGap = gTimeStamp[i] - fTimeRef;
                fTimeGapQRT = gTimeStamp[i] - fTimeRefQRT;
            } else {
                fTimeGap = gTimeStamp[i] - gTimeStamp[pLastEvent];
                fTimeGapQRT = gTimeStamp[i] - fTimeLastQRT;
            }

#ifdef DEBUG
            printf("%8d  %1d  %1d  %6d  %6d  %10d  ", i + 1, gQRT[i], gHelicity_rep[i], Int_t(fTimeGap), Int_t(fTimeGapQRT), Int_t(fTimeRefQRT));
#endif

            // predict helicity
            if (fNSeedTIR == MAXBIT) {
                if (gQRT[i] == 1) { // possibly a new pattern
                    if (gQRT[pLastEvent] == 0) { // the first event of a new pattern
                        Int_t MissedQRT = Int_t(fTimeGap / (4 * WT));

                        if (MissedQRT < MAXMISSED) {
                            for (Int_t j = 0; j <= MissedQRT; j++) {
                                fPolarityTIR_rep = RanBit30(fSeedTIR_rep);
                                fPolarityTIR_act = RanBit30(fSeedTIR_act);

                                if (fTimeRef > 0) {
                                    fTimeRef += 4 * WT;
                                    fTimeRefQRT += 4 * WT;
                                }
                            }

#ifdef DEBUG
                            printf("D1  %2d  %08x\n", MissedQRT, fSeedTIR_rep);
#endif

                            if (gHelicity_rep[i] == fPolarityTIR_rep)
                                fPhaseTIR_rep = 0;
                            else
                                fPhaseTIR_rep = 4;
                        } else
                            fPhaseTIR_rep = 4;
                    } else {
                        // 2 possibilities
                        // 1) missed 3 or 3+4n windows between 2 QRT windows
                        // 2) 2 events in the same QRT window
                        Int_t MissedQRT;

                        for (MissedQRT = 0; MissedQRT < MAXMISSED; MissedQRT++) {
                            if (fTimeGapQRT < ((2.0 + MissedQRT * 4) * WT))
                                break;

                            fPolarityTIR_rep = RanBit30(fSeedTIR_rep);
                            fPolarityTIR_act = RanBit30(fSeedTIR_act);

                            if (fTimeRef > 0) {
                                fTimeRef += 4 * WT;
                                fTimeRefQRT += 4 * WT;
                            }
                        }

#ifdef DEBUG
                        printf("D2  %2d  %08x\n", MissedQRT, fSeedTIR_rep);
#endif

                        if (MissedQRT < MAXMISSED) {
                            if (gHelicity_rep[i] == fPolarityTIR_rep)
                                fPhaseTIR_rep = 0;
                            else
                                fPhaseTIR_rep = 4;
                        } else
                            fPhaseTIR_rep = 4;
                    }
                } else {
                    Int_t MissedQRT = Int_t(fTimeGapQRT / (4 * WT));

                    if (fTimeRefQRT > 0) {
                        Double_t fRefQRT_tmp = fTimeRefQRT + (MissedQRT + 1) * 4 * WT - 0.4 * WT;

                        if (fRefQRT_tmp < gTimeStamp[i])
                            MissedQRT += 1;
                    }

                    if (MissedQRT < MAXMISSED) {
                        for (Int_t j = 0; j < MissedQRT; j++) {
                            fPolarityTIR_rep = RanBit30(fSeedTIR_rep);
                            fPolarityTIR_act = RanBit30(fSeedTIR_act);
                        }

                        if (gHelicity_rep[i] == fPolarityTIR_rep) {
                            if (gPairSync[i] == 0)
                                fPhaseTIR_rep = 3;
                            else
                                fPhaseTIR_rep = 4;
                        } else {
                            if (gPairSync[i] == 0)
                                fPhaseTIR_rep = 1;
                            else
                                fPhaseTIR_rep = 2;
                        }

                        if ((fSeedTIR_rep == gSeed_rep[pLastEvent]) && (fPhaseTIR_rep < fPhaseLastTIR)) {
                            MissedQRT += 1;
                            fPolarityTIR_rep = RanBit30(fSeedTIR_rep);
                            fPolarityTIR_act = RanBit30(fSeedTIR_act);

                            if (gHelicity_rep[i] == fPolarityTIR_rep) {
                                if (gPairSync[i] == 0)
                                    fPhaseTIR_rep = 3;
                                else
                                    fPhaseTIR_rep = 4;
                            } else {
                                if (gPairSync[i] == 0)
                                    fPhaseTIR_rep = 1;
                                else
                                    fPhaseTIR_rep = 2;
                            }
                        }
                    } else
                        fPhaseTIR_rep = 4;

#ifdef DEBUG
                    printf("D3  %2d  %08x\n", MissedQRT, fSeedTIR_rep);
#endif
                    fTimeLastQRT = fTimeLastQRT + MissedQRT * 4.0 * WT;

                    if (fTimeRef > 0) {
                        fTimeRef += MissedQRT * 4 * WT;
                        fTimeRefQRT += MissedQRT * 4 * WT;
                    }
                }

                if (fPhaseTIR_rep >= 4) {
                    fNSeedTIR = 0;
                    fSeedTIR_rep = 0;
                    fTimeRef = -1;
                    fTimeRefQRT = -1;
                    fTimeGap = gTimeStamp[i] - gTimeStamp[pLastEvent];
                    fTimeGapQRT = gTimeStamp[i] - fTimeLastQRT;
#ifdef DEBUG
                    printf("%8d  %1d  %1d  %6d  %6d  %10d  ", i + 1, gQRT[i], gHelicity_rep[i], Int_t(fTimeGap), Int_t(fTimeGapQRT), Int_t(fTimeRefQRT));
#endif
                }
            }

            //Catch random seed
            if (fNSeedTIR < MAXBIT) {
                if (gQRT[i] == 1) {
                    Int_t k = 0;

                    if (gQRT[pLastEvent] == 0) {
                        if ((fTimeGap < fTimeGapQRT) && (fTimeGapQRT < 6.0 * WT)) {
                            fSeedTIR_rep = ((fSeedTIR_rep << 1 & 0x3FFFFFFF) | gHelicity_rep[i]);
                            fNSeedTIR++;
                        } else {
                            fSeedTIR_rep = 0;
                            fNSeedTIR = 0;

                            if (usering) {
                                for (k = 0; k < MAXMISSED; k++) {
                                    if (fTimeGapQRT < (6.0 + k * 4) * WT)
                                        break;

                                    fSeedTIR_fake = (fSeedTIR_fake << 1 & 0x3FFFFFFF) | 0x0;
                                    fMask = (fMask << 1 & 0x3FFFFFFF) | 0x0;
                                }
                            }
                        }

                        if (usering) {
                            fSeedTIR_fake = (fSeedTIR_fake << 1 & 0x3FFFFFFF) | gHelicity_rep[i];
                            fMask = (fMask << 1 & 0x3FFFFFFF) | 0x1;
                        }

                        if (k >= MAXMISSED)
                            fSeedTIR_fake = 0;
                    } else if (gQRT[pLastEvent] == 1) {
                        if (fTimeGap > 2.0 * WT) {
                            if (fTimeGap < 6.0 * WT) {
                                fSeedTIR_rep = ((fSeedTIR_rep << 1 & 0x3FFFFFFF) | gHelicity_rep[i]);
                                fNSeedTIR++;
                            } else {
                                fNSeedTIR = 0;
                                fSeedTIR_rep = 0;

                                if (usering) {
                                    for (k = 0; k < MAXMISSED; k++) {
                                        if (fTimeGapQRT < (2.0 + k * 4) * WT)
                                            break;

                                        fSeedTIR_fake = (fSeedTIR_fake << 1 & 0x3FFFFFFF) | 0x0;
                                        fMask = (fMask << 1 & 0x3FFFFFFF) | 0x0;
                                    }
                                }
                            }

                            if (usering) {
                                fSeedTIR_fake = (fSeedTIR_fake << 1 & 0x3FFFFFFF) | gHelicity_rep[i];
                                fMask = (fMask << 1 & 0x3FFFFFFF) | 0x1;
                            }

                            if (k >= MAXMISSED)
                                fSeedTIR_fake = 0;
                        }
                    }

                    if ((usering) && (popcount(fMask) >= 10)) {
                        Int_t j;

                        for (j = 0; j < gNRing; j++) {
                            if (gEventRing[j] >= i + 1)
                                break;
                        }

                        for (Int_t k = j; k < j + 40; k++) {
                            if ((gSeedRing_rep[k] & fMask) == fSeedTIR_fake) {
                                fSeedTIR_rep = gSeedRing_rep[k];
                                fNSeedTIR = MAXBIT;
                                fMask = 0;
                                break;
                            }
                        }
                    }

                    if (fNSeedTIR == MAXBIT) {
                        fPolarityTIR_rep = fSeedTIR_rep & 0x01;
                        fSeedTIR_act = fSeedTIR_rep;

                        for (Int_t j = 0; j < NDELAY; j++)
                            fPolarityTIR_act = RanBit30(fSeedTIR_act);

                        fPhaseTIR_rep = 0;
                        gError[i] = 0;
                    }

#ifdef DEBUG
                    printf("S1  %2d  %08x  %08x\n", fNSeedTIR, fSeedTIR_rep, fMask);
#endif
                } else {
#ifdef DEBUG
                    printf("S2  %2d  %08x  %08x\n", fNSeedTIR, fSeedTIR_rep, fMask);
#endif
                }
            }

            // assign actual helicity to event
            if (fNSeedTIR == MAXBIT) {
                if (fPhaseTIR_rep <= 3 || fPhaseTIR_rep >= 0) {
                    if (fPolarityTIR_act == 1) {
                        if (fPhaseTIR_rep == 0 || fPhaseTIR_rep == 3)
                            gHelicity_act[i] = 1;
                        else
                            gHelicity_act[i] = -1;
                    } else {
                        if (fPhaseTIR_rep == 0 || fPhaseTIR_rep == 3)
                            gHelicity_act[i] = -1;
                        else
                            gHelicity_act[i] = 1;
                    }

                    if (pLastMPS > pLastEvent) {
                        if (((gTimeStamp[i] - gTimeStamp[pLastMPS]) < 0.75 * WT) && (gTimeStamp[pLastMPS] < gTimeStamp[i]) && (gHelicity_rep[pLastMPS] == gHelicity_rep[i]) && (gQRT[pLastMPS] == gQRT[i]) && (gPairSync[pLastMPS] == gPairSync[i])) {
                            fTimeRef = gTimeStamp[pLastMPS] + (3 - fPhaseTIR_rep) * WT;
                            fTimeRefQRT = fTimeRef - 2 * WT;
                        }
                    }
                } else {
                    gError[i] |= 0x02;
                    gHelicity_act[i] = 0;
                }

                gSeed_rep[i] = fSeedTIR_rep;
                fPhaseLastTIR = fPhaseTIR_rep;
            } else {
                gError[i] |= 0x01;
                gSeed_rep[i] = 0;
                gHelicity_act[i] = 0;
                fPhaseLastTIR = -1;
            }

            if (gQRT[i] == 1)
                fTimeLastQRT = gTimeStamp[i];

            pLastEvent = i;
        } else {
            // non-physics trigger event
            gError[i] |= 0x08;
            gSeed_rep[i] = 0;
            gHelicity_act[i] = 0;
            pLastMPS = i;
        }
    }

    return 0;
}

Int_t printout(Int_t nrun, Bool_t usering)
{
    if ((fp1 = fopen(Form("%s/helTIR_%d.decode.dat", INDIR, nrun), "r")) == NULL) {
        fprintf(stderr, "Can not open %s/helTIR_%d.decode.dat\n", INDIR, nrun);
        exit(-1);
    }

    if ((fp2 = fopen(Form("%s/helTIR_%d.noring.dat", OUTDIR, nrun), "w")) == NULL) {
        fprintf(stderr, "Can not open %s/helTIR_%d.noring.dat\n", OUTDIR, nrun);
        exit(-1);
    };

    Int_t tempi[10];

    fprintf(fp2, "%d\n", gN);

    for (Int_t i = 0; i < gN; i++) {
        fscanf(fp1, "%d%d%d%d%d%d%d%d", &tempi[0], &tempi[1], &tempi[2], &tempi[3], &tempi[4], &tempi[5], &tempi[6], &tempi[7]);
        fprintf(fp2, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%08x\t%d\t%d\t%d\n", tempi[0], gHelicity_act[i], gHelicity_rep[i], gQRT[i], gPairSync[i], gMPS[i], gTimeStamp[i], gSeed_rep[i], gError[i], tempi[6], tempi[7]);
    }

    delete[] gHelicity_act;
    delete[] gHelicity_rep;
    delete[] gPairSync;
    delete[] gQRT;
    delete[] gMPS;
    delete[] gTimeStamp;
    delete[] gSeed_rep;
    delete[] gError;

    if (usering) {
        delete[] gSeedRing_rep;
        delete[] gEventRing;
    }

    fclose(fp1);
    fclose(fp2);

    return 0;
}

Int_t RanBit30(Int_t &ranseed)
{
    // Take 7,28,29,30 bit of ranseed out

    UInt_t bit7 = ((ranseed & 0x00000040) != 0);
    UInt_t bit28 = ((ranseed & 0x08000000) != 0);
    UInt_t bit29 = ((ranseed & 0x10000000) != 0);
    UInt_t bit30 = ((ranseed & 0x20000000) != 0);
    UInt_t newbit = (bit30 ^ bit29 ^ bit28 ^ bit7) & 0x1;

    if (ranseed <= 0)
        newbit = 0;

    ranseed = ((ranseed << 1) | newbit) & 0x3FFFFFFF;

    return newbit;
}

Int_t popcount(Int_t x)
{
    x = (x & M1) + ((x >> 1) & M1);
    x = (x & M2) + ((x >> 2) & M2);
    x = (x & M4) + ((x >> 4) & M4);
    x = (x & M8) + ((x >> 8) & M8);
    x = (x & M16) + ((x >> 16) & M16);

    return x;
}

void usage(Int_t argc, Char_t **argv)
{
    printf("usage: %s [options] RUN_NUMBER\n", argv[0]);
    printf("  -c, --cfgfile=config.cfg   Set configuration file name\n");
    printf("  -h, --help                 This small usage guide\n");
    printf("  -i, --indir=.              Set input directory\n");
    printf("  -o, --outdir=.             Set output directory\n");
    printf("  -r, --ring                 Use ring info to help predict\n");
}
