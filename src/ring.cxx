#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <getopt.h>

#include <libconfig.h>

#include "TROOT.h"

FILE *fp1, *fp2;

// global variables
Int_t NRING = 0;
Int_t NHAPPEX = 0;
Int_t NDELAY;
Int_t MAXBIT;
Double_t WT;
Bool_t USEHAPPEX;

Int_t *gHelicity_rep;
Int_t *gHelicity_act;
Int_t *gQRT;
Int_t *gSeed_rep;
Int_t *gError;
Int_t gN;

Int_t readin(Int_t nrun, Int_t nring, Int_t select);
Int_t predictring(Int_t select);
Int_t delayring(Int_t ndelay, Int_t select);
Int_t printout(Int_t nrun, Int_t nring, Int_t select);
Int_t RanBit30(Int_t &runseed);
Int_t BitRan30(Int_t &runseed); // reversal prediction
void usage(Int_t argc, Char_t **argv);

Char_t CFGFILE[300] = "./config.cfg";
Char_t INDIR[300] = ".";
Char_t OUTDIR[300] = ".";

Int_t main(Int_t argc, Char_t **argv)
{
    Int_t c;

    while (1) {
        static struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"cfgfile", required_argument, 0, 'c'},
            {"indir", required_argument, 0, 'i'},
            {"outdir", required_argument, 0, 'o'},
            {0, 0, 0, 0}
        };

        Int_t option_index = 0;
        c = getopt_long(argc, argv, "c:hi:o:", long_options, &option_index);

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

    if (!(config_lookup_int(&cfg, "ndelay", &NDELAY) && config_lookup_int(&cfg, "maxbit", &MAXBIT) && config_lookup_float(&cfg, "windowlength", &WT)))
        configerror = kTRUE;

    setting = config_lookup(&cfg, "ringinfo.data");

    Int_t delayRIN;

    if (setting != NULL) {
        NRING = config_setting_length(setting);
        config_lookup_int(&cfg, "ringinfo.delay", &delayRIN);
    } else
        configerror = kTRUE;


    setting = config_lookup(&cfg, "happexinfo.data");

    Int_t delayHAP;

    if (setting != NULL) {
        USEHAPPEX = kTRUE;
        NHAPPEX = config_setting_length(setting);
        config_lookup_int(&cfg, "happexinfo.delay", &delayHAP);
    } else
        USEHAPPEX = kFALSE;

    if (configerror) {
        fprintf(stderr, "Invalid cfg file\n");
        exit(-1);
    }

    clock_t start, end;

    start = clock();
    readin(nrun, NRING, 1);
    predictring(1);
    delayring(delayRIN, 1);
    printout(nrun, NRING, 1);
    end = clock();

    printf("Actual helicity calculated in %5.3f s\n", (Double_t)(end - start) / (Double_t) CLOCKS_PER_SEC);

    if (USEHAPPEX) {
        start = clock();
        readin(nrun, NHAPPEX, 2);
        predictring(2);
        delayring(delayHAP, 2);
        printout(nrun, NHAPPEX, 2);
        end = clock();

        printf("Actual helicity calculated in %5.3f s\n", (Double_t)(end - start) / (Double_t) CLOCKS_PER_SEC);
    }

    config_destroy(&cfg);

    return 0;
}

Int_t readin(Int_t nrun, Int_t nring, Int_t select)
{
    if (select == 1) {
        printf("Reading scaler ring buffer helicity information ...\n");

        if ((fp1 = fopen(Form("%s/helRIN_%d.decode.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helRIN_%d.decode.dat\n", INDIR, nrun);
            exit(-1);
        }
    } else if (select == 2) {
        printf("Reading happex ring buffer helicity information ...\n");

        if ((fp1 = fopen(Form("%s/helHAP_%d.decode.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helHAP_%d.decode.dat\n", INDIR, nrun);
            exit(-1);
        }
    }

    Int_t tempi[20];
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
    gSeed_rep = new Int_t[gN];
    gError = new Int_t[gN];

    memset(gHelicity_rep, 0, sizeof(Int_t) * gN);
    memset(gHelicity_act, 0, sizeof(Int_t) * gN);
    memset(gQRT, 0, sizeof(Int_t) * gN);
    memset(gSeed_rep, 0, sizeof(Int_t) * gN);
    memset(gError, 0, sizeof(Int_t) * gN);

    rewind(fp1);

    for (Int_t i = 0; i < gN; i++) {
        fscanf(fp1, "%d%d%d", &tempi[0], &gHelicity_rep[i], &gQRT[i]);
        fgets(tempc, 300, fp1);
    }

    fclose(fp1);

    return 0;
}

Int_t predictring(Int_t select)
{
    printf("Calculating actual helicity information ...\n");

    Int_t fPhaseRing_rep = 0;
    Int_t fPolarityRing_rep = 0, fPolarityRing_act = 0;
    Int_t fSeedRing_rep = 0, fSeedRing_act = 0;
    Int_t fNSeedRing = 0;

    for (Int_t i = 0; i < gN; i++) {
        gError[i] = 0;
        gSeed_rep[i] = 0;
        gHelicity_act[i] = 0;

        if (gQRT[i] == 1)
            fPhaseRing_rep = 0;
        else
            fPhaseRing_rep += 1;

        if (fPhaseRing_rep >= 4) {
            fNSeedRing = 0;
            gError[i] |= (0x02 << (select * 4));
        }

        if ((fNSeedRing == MAXBIT) && gQRT[i] == 1) {
            fPolarityRing_rep = RanBit30(fSeedRing_rep);
            fPolarityRing_act = RanBit30(fSeedRing_act);

            if (fPolarityRing_rep != gHelicity_rep[i]) {
                fNSeedRing = 0;
                gError[i] |= (0x02 << (select * 4));
            }
        }

        if ((fNSeedRing < MAXBIT) && gQRT[i] == 1) {
            fSeedRing_rep = ((fSeedRing_rep << 1) & 0x3FFFFFFF) | gHelicity_rep[i];
            fNSeedRing += 1;

            if (fNSeedRing == MAXBIT) {
                fSeedRing_act = fSeedRing_rep;

                for (Int_t j = 0; j < NDELAY; j++)
                    fPolarityRing_act = RanBit30(fSeedRing_act);

                gError[i] = 0;
            }
        }

        if (fNSeedRing == MAXBIT) {
            if (fPolarityRing_act == 1) {
                if (fPhaseRing_rep == 0 || fPhaseRing_rep == 3)
                    gHelicity_act[i] = 1;
                else
                    gHelicity_act[i] = -1;
            } else {
                if (fPhaseRing_rep == 0 || fPhaseRing_rep == 3)
                    gHelicity_act[i] = -1;
                else
                    gHelicity_act[i] = 1;
            }

            gSeed_rep[i] = fSeedRing_rep;
        } else {
            gError[i] |= (0x01 << (select * 4));
            gSeed_rep[i] = 0;
            gHelicity_act[i] = 0;
        }
    }

    for (Int_t i = gN - 2; i >= 0; i--) {
        if (gQRT[i + 1] == 1)
            fPhaseRing_rep = 3;
        else
            fPhaseRing_rep -= 1;

        if (fPhaseRing_rep >= 0) {
            if (gError[i] == 0) {
                if (gQRT[i] == 1)
                    fSeedRing_rep = gSeed_rep[i];
            } else {
                if (gQRT[i] == 1) {
                    if (fSeedRing_rep != 0) {
                        fPolarityRing_rep = BitRan30(fSeedRing_rep);

                        if (gHelicity_rep[i] == fPolarityRing_rep) {
                            fSeedRing_act = fSeedRing_rep;

                            for (Int_t j = 0; j < NDELAY; j++)
                                fPolarityRing_act = RanBit30(fSeedRing_act);

                            if (fPolarityRing_act == 1) {
                                gHelicity_act[i] = gHelicity_act[i + 3] = 1;
                                gHelicity_act[i + 1] = gHelicity_act[i + 2] = -1;

                                for (Int_t j = 0; j < 4; j++) {
                                    gError[i + j] = 0;
                                    gSeed_rep[i + j] = fSeedRing_rep;
                                }
                            } else {
                                gHelicity_act[i] = gHelicity_act[i + 3] = -1;
                                gHelicity_act[i + 1] = gHelicity_act[i + 2] = 1;

                                for (Int_t j = 0; j < 4; j++) {
                                    gError[i + j] = 0;
                                    gSeed_rep[i + j] = fSeedRing_rep;
                                }
                            }
                        } else
                            fSeedRing_rep = 0;
                    }
                }
            }
        } else
            fSeedRing_rep = 0;
    }

    return 0;
}

Int_t delayring(Int_t ndelay, Int_t select)
{
    if (ndelay > 0) {
        for (Int_t i = 0; i < gN - ndelay; i++) {
            gHelicity_rep[i] = gHelicity_rep[i + ndelay];
            gHelicity_act[i] = gHelicity_act[i + ndelay];
            gQRT[i] = gQRT[i + ndelay];
            gSeed_rep[i] = gSeed_rep[i + ndelay];
            gError[i] = gError[i + ndelay];
        }

        for (Int_t i = gN - ndelay; i < gN; i++) {
            gHelicity_rep[i] = 0;
            gHelicity_act[i] = 0;
            gQRT[i] = 0;
            gSeed_rep[i] = 0;
            gError[i] |= (0x02 << (select * 4));
        }
    } else if (ndelay < 0) {
        for (Int_t i = gN - 1; i >= -ndelay; i--) {
            gHelicity_rep[i] = gHelicity_rep[i + ndelay];
            gHelicity_act[i] = gHelicity_act[i + ndelay];
            gQRT[i] = gQRT[i + ndelay];
            gSeed_rep[i] = gSeed_rep[i + ndelay];
            gError[i] = gError[i + ndelay];
        }

        for (Int_t i = -ndelay - 1; i >= 0; i--) {
            gHelicity_rep[i] = 0;
            gHelicity_act[i] = 0;
            gQRT[i] = 0;
            gSeed_rep[i] = 0;
            gError[i] |= (0x02 << (select * 4));
        }
    }

    return 0;
}

Int_t printout(Int_t nrun, Int_t nring, Int_t select)
{
    if (select == 1) {
        if ((fp1 = fopen(Form("%s/helRIN_%d.decode.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helRIN_%d.decode.dat\n", INDIR, nrun);
            exit(-1);
        }

        if ((fp2 = fopen(Form("%s/helRIN_%d.nalign.dat", OUTDIR, nrun), "w")) == NULL) {
            fprintf(stderr, "Can not open %s/helRIN_%d.dat\n", OUTDIR, nrun);
            exit(-1);
        }
    } else if (select == 2) {
        if ((fp1 = fopen(Form("%s/helHAP_%d.decode.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helHAP_%d.decode.dat\n", INDIR, nrun);
            exit(-1);
        }

        if ((fp2 = fopen(Form("%s/helHAP_%d.nalign.dat", OUTDIR, nrun), "w")) == NULL) {
            fprintf(stderr, "Can not open %s/helHAP_%d.dat\n", OUTDIR, nrun);
            exit(-1);
        }
    }

    Int_t tempi[10];

    fprintf(fp2, "%d\n", gN);

    for (Int_t i = 0; i < gN; i++) {
        fscanf(fp1, "%d%d%d", &tempi[0], &tempi[1], &tempi[2]);
        fprintf(fp2, "%d\t%d\t%d\t%d\t%08x\t%d", tempi[0], gHelicity_act[i], gHelicity_rep[i], gQRT[i], gSeed_rep[i], gError[i]);

        for (Int_t k = 0; k < nring; k++) {
            fscanf(fp1, "%d", &tempi[0]);
            fprintf(fp2, "\t%d", tempi[0]);
        }

        fprintf(fp2, "\n");
    }

    fclose(fp1);
    fclose(fp2);

    delete[] gHelicity_rep;
    delete[] gHelicity_act;
    delete[] gQRT;
    delete[] gSeed_rep;
    delete[] gError;

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

Int_t BitRan30(Int_t &ranseed)
{
    // Take 1,8,29,30 bit of ranseed out
    // Backward predict

    UInt_t bit1 = ((ranseed & 0x00000001) != 0);
    UInt_t bit8 = ((ranseed & 0x00000080) != 0);
    UInt_t bit29 = ((ranseed & 0x10000000) != 0);
    UInt_t bit30 = ((ranseed & 0x20000000) != 0);
    UInt_t newbit = (bit30 ^ bit29 ^ bit8 ^ bit1) & 0x1;

    if (ranseed <= 0)
        newbit = 0;

    ranseed = ((ranseed >> 1) | (newbit << 29)) & 0x3FFFFFFF;
    newbit = ((ranseed & 0x1) != 0);

    return newbit;
}

void usage(Int_t argc, Char_t **argv)
{
    printf("usage: %s [options] RUN_NUMBER\n", argv[0]);
    printf("  -c, --cfgfile=config.cfg   Set configuration file name\n");
    printf("  -h, --help                 This small usage guide\n");
    printf("  -i, --indir=.              Set input directory\n");
    printf("  -o, --outdir=.             Set output directory\n");
}
