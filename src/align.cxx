#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <libconfig.h>

#include "TROOT.h"
#include "TSystem.h"

#include "hel.h"

#define NDATA 8
#define NBUFFER 400
#define CUTOFF1 120
#define CUTOFF2 120

FILE *fp1, *fp2;

//Global variables
Int_t gHel_act[LENTIR];
Int_t gSeed_rep[LENTIR];
Int_t gData[NDATA][LENTIR];
Int_t gError[LENTIR];
Int_t gHelRing_act[LENRIN];
Int_t gSeedRing_rep[LENRIN];
Int_t gDataRing[NDATA][LENRIN];
Int_t gErrorRing[LENRIN];
Int_t gN, gNRing;

Int_t readin(Int_t nrun, Int_t nring, Int_t select);
Int_t check(Int_t index, Int_t threshold, Int_t select);
Int_t align(Int_t nring, Int_t select);
Int_t printout(Int_t nrun, Int_t nring, Int_t select);
void usage(int argc, char** argv);

Char_t CFGFILE[300] = "./config.cfg";
Char_t INDIR[300] = ".";
Char_t OUTDIR[300] = ".";

int main(int argc, char** argv)
{
    int c;

    while (1) {
        static struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"cfgfile", required_argument, 0, 'c'},
            {"indir", required_argument, 0, 'i'},
            {"outdir", required_argument, 0, 'o'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "c:hi:o:", long_options, &option_index);

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
    } else {
        usage(argc, argv);
        exit(-1);
    }

    Int_t fIndexRing, fIndexHappex;
    Int_t fThrRing, fThrHappex;

    config_t cfg;
    config_setting_t *setting;

    config_init(&cfg);

    if (!config_read_file(&cfg, CFGFILE)) {
        fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg), config_error_line(&cfg), config_error_text(&cfg));
        config_destroy(&cfg);
        exit(EXIT_FAILURE);
    }

    Bool_t configerror = kFALSE;

    setting = config_lookup(&cfg, "ringinfo.data");
    if (setting != NULL) {
        NRING = config_setting_length(setting);
        config_lookup_int(&cfg, "ringinfo.bcm.index", &fIndexRing);
        config_lookup_int(&cfg, "ringinfo.bcm.threshold", &fThrRing);
    } else
        configerror = kTRUE;

    setting = config_lookup(&cfg, "happexinfo.data");
    if (setting != NULL) {
        USEHAPPEX = kTRUE;
        NHAPPEX = config_setting_length(setting);
        config_lookup_int(&cfg, "happexinfo.bcm.index", &fIndexHappex);
        config_lookup_int(&cfg, "happexinfo.bcm.threshold", &fThrHappex);
    } else {
        USEHAPPEX = kFALSE;
    }

    if (configerror) {
        fprintf(stderr, "Invalid cfg file\n");
        exit(-1);
    }

    readin(nrun, NRING, 1);
    check(fIndexRing, fThrRing, 1);
    align(NRING, 1);
    printout(nrun, NRING, 1);

    if (USEHAPPEX) {
        gSystem->Rename(Form("%s/hel_%d.dat", OUTDIR, nrun), Form("%s/helTIR_%d.nohapp.dat", INDIR, nrun));
        readin(nrun, NHAPPEX, 2);
        check(fIndexHappex, fThrHappex, 2);
        align(NHAPPEX, 2);
        printout(nrun, NHAPPEX, 2);
    }

    return 0;
}

Int_t readin(Int_t nrun, Int_t nring, Int_t select)
{
    if (select == 1) {
        if ((fp1 = fopen(Form("%s/helTIR_%d.noring.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helTIR_%d.noring.dat", INDIR, nrun);
            exit(-1);
        }
        if ((fp2 = fopen(Form("%s/helRIN_%d.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helRIN_%d.dat", INDIR, nrun);
            exit(-1);
        }
    } else if (select == 2) {
        if ((fp1 = fopen(Form("%s/helTIR_%d.nohapp.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helTIR_%d.nohapp.dat", INDIR, nrun);
            exit(-1);
        }
        if ((fp2 = fopen(Form("%s/helHAP_%d.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helHAP_%d.dat", INDIR, nrun);
            exit(-1);
        }
    }

    Int_t tempi[10] = {0};
    Char_t tempc[300];

    printf("Reading TIR helicity information ...\n");

    for (Int_t i = 0; i < LENTIR; i++) {
        gHel_act[i] = 0;
        gSeed_rep[i] = 0;
        gError[i] = 0;
        for (Int_t k = 0; k < NDATA; k++) gData[k][i] = 0;
    }

    fscanf(fp1, "%d", &gN);
    for (Int_t i = 0; i < gN; i++) {
        fscanf(fp1, "%d%d%d%d%d%d%d%x%d", &tempi[0], &gHel_act[i], &tempi[2], &tempi[3], &tempi[4], &tempi[5], &tempi[6], &gSeed_rep[i], &gError[i]);
        fgets(tempc, 300, fp1);
    }

    if (select == 1) {
        printf("Reading scaler ring buffer helicity information ...\n");
    } else if (select == 2) {
        printf("Reading happex ring buffer helicity information ...\n");
    }

    for (Int_t i = 0; i < LENRIN; i++) {
        gHelRing_act[i] = 0;
        gSeedRing_rep[i] = 0;
        gErrorRing[i] = 0;
        for (Int_t k = 0; k < NDATA; k++) gDataRing[k][i] = 0;
    }

    Int_t fPhaseRing = 0;
    Int_t fLastSeedRing = 0;
    fscanf(fp2, "%d", &gNRing);
    for (Int_t i = 0; i < gNRing; i++) {
        fscanf(fp2, "%d%d%d%d%x%d", &tempi[0], &gHelRing_act[i], &tempi[2], &tempi[3], &gSeedRing_rep[i], &gErrorRing[i]);
        for (Int_t k = 0; k < nring; k++)
            fscanf(fp2, "%d", &gDataRing[k][i]);

        if (gSeedRing_rep[i] == fLastSeedRing) fPhaseRing++;
        else fPhaseRing = 0;
        if (fPhaseRing > 3) gErrorRing[i] |= (0x2 << (4 * select));
        fLastSeedRing = gSeedRing_rep[i];
    }
    // last pattern may have some problem, correct it
    if (fPhaseRing < 3) gErrorRing[gNRing - 1] |= (0x2 << (4 * select));

    fclose(fp1);
    fclose(fp2);

    return 0;
}

Int_t check(Int_t index, Int_t threshold, Int_t select)
{
    printf("Checking TIR helicity information ...\n");

    if (select == 1) {
        printf("Checking scaler ring buffer helicity information ...\n");
    } else if (select == 2) {
        printf("Checking happex ring buffer helicity information ...\n");
    }

    // check ring buffer helicity information
    Int_t l = 0;
    while (l < gNRing) {
        if (gDataRing[index][l] < threshold) { // check whether beam is off
            Int_t m = l;
            for (Int_t i = 0; (i < CUTOFF1)&&(m >= 0); i++, m--) gErrorRing[m] |= ((0x4 << (4 * select))+(0x1000 * select));
            Int_t msave = m;
            while ((gSeedRing_rep[m] == gSeedRing_rep[msave])&&(m >= 0)) {
                gErrorRing[m] |= ((0x4 << (4 * select))+(0x1000 * select));
                m--;
            }
            Int_t p = l;
            while ((gDataRing[index][p] < threshold)&&(p < gNRing)) {
                gErrorRing[p] |= ((0x4 << (4 * select))+(0x1000 * select));
                p++;
            }
            for (Int_t i = 0; (i < CUTOFF2)&&(p < gNRing); i++, p++) gErrorRing[p] |= ((0x4 << (4 * select))+(0x1000 * select));
            Int_t psave = p;
            while ((gSeedRing_rep[p] == gSeedRing_rep[psave])&&(p < gNRing)) {
                gErrorRing[p] |= ((0x4 << (4 * select))+(0x1000 * select));
                p++;
            }
            l = p;
        } else if (gErrorRing[l] != 0) { // if a quartet contains bad event, mark the whole quartet to be bad
            Int_t m = l;
            while ((gSeedRing_rep[m] == gSeedRing_rep[l])&&(m >= 0)) {
                gErrorRing[m] |= (gErrorRing[l] | (0x1000 * select));
                m--;
            }
            Int_t p = l + 1;
            while ((gSeedRing_rep[p] == gSeedRing_rep[l])&&(p < gNRing)) {
                gErrorRing[p] |= (gErrorRing[l] | (0x1000 * select));
                p++;
            }
            l = p;
        } else {
            l++;
        }
    }

    // check tir helicity information
    // if a quartet contains bad event, mark the whole quartet to be bad
    l = 0;
    Int_t pRing = 0;
    while (l < gN) {
        if (gSeed_rep[l] != 0) {
            while ((gSeedRing_rep[pRing] != gSeed_rep[l])&&(pRing < gNRing)) pRing++;
            if (pRing < gNRing) {
                if (gErrorRing[pRing] != 0) gError[l] |= (gErrorRing[pRing] | (0x1000 * select));
            } else {
                pRing = 0;
                gError[l] |= (0x1000 * select);
            }
        }
        if ((gError[l]&0xFF7) != 0) {
            Int_t m = l;
            while ((gSeed_rep[m] == gSeed_rep[l])&&(m >= 0)) {
                gError[m] |= (gError[l] | (0x1000 * select));
                m--;
            }
            Int_t p = l + 1;
            while ((gSeed_rep[p] == gSeed_rep[l])&&(p < gN)) {
                gError[p] |= (gError[l] | (0x1000 * select));
                p++;
            }
            l = p;
        } else {
            if ((gError[l]&0x008) == 0x0008) {
                gError[l] |= (0x1000 * select);
            }
            l++;
        }
    }

    return 0;
}

Int_t align(Int_t nring, Int_t select)
{
    if (select == 1) {
        printf("Align TIR helicity information with scaler helicity information ...\n");
    } else if (select == 2) {
        printf("Align TIR helicity information with HAPPEX helicity information ...\n");
    }

    Int_t pStart = 0, pStartRing = 0;
    Int_t l = 0;
    while ((pStartRing == 0)&&(pStart < gN)) {
        while (((gSeed_rep[pStart] == gSeed_rep[l]) || ((gError[pStart]&0xFFF) != 0))&&(pStart < gNRing)) pStart++;
        Int_t m = 0;
        while ((gSeedRing_rep[m] != gSeed_rep[pStart])&&(m < gNRing)) m++;
        if (m < gNRing) pStartRing = m;
        else l = pStart;
    }

    for (Int_t i = 0; i < pStart; i++) gError[i] |= (0x1000 * select);
    for (Int_t i = 0; i < pStartRing; i++) gErrorRing[i] |= (0x1000 * select);

    Int_t pGoodP = 0, pGoodN = 0;
    Int_t fDataP[NDATA], fDataN[NDATA];
    Int_t fSaved = 0;
    Int_t fSaveP[NDATA] = {0}, fSaveN[NDATA] = {0};

    Int_t pRing = pStartRing, pLastRing = pStartRing;
    while (pRing < gNRing) {
        if (gErrorRing[pRing] == 0) {
            // normal events
            pLastRing = pRing;
            while ((gSeedRing_rep[pRing] == gSeedRing_rep[pLastRing])&&(pRing < gNRing)) pRing++;
            if ((pRing - pLastRing) == 4) {
                // good pattern
                for (Int_t i = 0; i < NDATA; i++) {
                    fDataP[i] = 0;
                    fDataN[i] = 0;
                }
                for (Int_t i = pLastRing; i < pRing; i++) {
                    if (gHelRing_act[i] == 1) {
                        for (Int_t j = 0; j < nring; j++) fDataP[j] += gDataRing[j][i];
                    } else if (gHelRing_act[i] == -1) {
                        for (Int_t j = 0; j < nring; j++) fDataN[j] += gDataRing[j][i];
                    }
                }

                // check if tir has same pattern
                Int_t IsGood = 0;
                Int_t m = pStart;
                while ((((gError[m]&0xFFF) != 0) || (gSeed_rep[m] != gSeedRing_rep[pLastRing]))&&(m < gN)) m++;
                if (m < gN) {
                    // found the same pattern
                    Int_t n = m;

                    while ((((gError[m]&0xFFF) != 0) || (gSeed_rep[n] == gSeed_rep[m]))&&(n < gN)) {
                        if (gHel_act[n] == 1) IsGood |= 0x01;
                        if (gHel_act[n] == -1) IsGood |= 0x10;
                        n++;
                    }
                    if (IsGood == 0x11) {
                        // good tir pattern
                        Int_t Found = 0;
                        for (Int_t i = m; i < n; i++) {
                            if (((Found & 0x01) == 0)&&(gHel_act[i] == 1)) {
                                pGoodP = i;
                                Found |= 0x01;
                            }
                            if (((Found & 0x10) == 0)&&(gHel_act[i] == -1)) {
                                pGoodN = i;
                                Found |= 0x10;
                            }
                            if (Found == 0x11) break;
                        }
                        for (Int_t i = 0; i < NDATA; i++) {
                            gData[i][pGoodP] += (fDataP[i] + fSaveP[i]);
                            fSaveP[i] = 0;
                            gData[i][pGoodN] += (fDataN[i] + fSaveN[i]);
                            fSaveN[i] = 0;
                            fSaved = 0;
                        }
                    }
#ifdef DEBUG
                    printf("%8d  %8d  %08x  %8d  %8d  %8d\n", m, n, gSeed_rep[m], pRing, pGoodP, pGoodN);
#endif
                }
                if (IsGood != 0x11) {
                    for (Int_t i = 0; i < NDATA; i++) {
                        fSaveP[i] += fDataP[i];
                        fSaveN[i] += fDataN[i];
                    }
                    fSaved++;
                }
            } else {
                // something very bad happened, should never reach here after check()
            }
        } else {
            pRing++;
        }
        if (fSaved > 0) {
            for (Int_t i = 0; i < NDATA; i++) {
                gData[i][pGoodP] += fSaveP[i];
                fSaveP[i] = 0;
                gData[i][pGoodN] += fSaveN[i];
                fSaveN[i] = 0;
                fSaved = 0;
            }
        }
    }

    return 0;
}

Int_t printout(Int_t nrun, Int_t nring, Int_t select)
{
    if (select == 1) {
        if ((fp1 = fopen(Form("%s/helTIR_%d.noring.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helTIR_%d.noring.dat", INDIR, nrun);
            exit(-1);
        }
    } else if (select == 2) {
        if ((fp1 = fopen(Form("%s/helTIR_%d.nohapp.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helTIR_%d.nohapp.dat", INDIR, nrun);
            exit(-1);
        }
    }
    if ((fp2 = fopen(Form("%s/hel_%d.dat", OUTDIR, nrun), "w")) == NULL) {
        fprintf(stderr, "Can not open %s/hel_%d.dat", OUTDIR, nrun);
        exit(-1);
    }

    Int_t tempi[10];
    Char_t tempc[300];

    fscanf(fp1, "%d", &tempi[0]);
    fprintf(fp2, "%d\n", gN);
    for (Int_t i = 0; i < gN; i++) {
        fscanf(fp1, "%d%d%d%d%d%d%d%x%d", &tempi[0], &tempi[1], &tempi[2], &tempi[3], &tempi[4], &tempi[5], &tempi[6], &tempi[7], &tempi[8]);
        fgets(tempc, 300, fp1);
        tempc[strlen(tempc) - 1] = '\0';
        fprintf(fp2, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%08x\t%d\t%s", tempi[0], tempi[1], tempi[2], tempi[3], tempi[4], tempi[5], tempi[6], tempi[7], gError[i], tempc);
        for (Int_t j = 0; j < nring; j++) {
            fprintf(fp2, "\t%d", gData[j][i]);
        }
        fprintf(fp2, "\n");
    }

    fclose(fp1);
    fclose(fp2);

    if (select == 1) {
        if ((fp1 = fopen(Form("%s/helRIN_%d.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helRIN_%d.dat", INDIR, nrun);
            exit(-1);
        }
        if ((fp2 = fopen(Form("%s/helRIN_%d.appcor.dat", OUTDIR, nrun), "w")) == NULL) {
            fprintf(stderr, "Can not open %s/helRIN_%d.appcor.dat", OUTDIR, nrun);
            exit(-1);
        }
    } else if (select == 2) {
        if ((fp1 = fopen(Form("%s/helHAP_%d.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helHAP_%d.dat", INDIR, nrun);
            exit(-1);
        }
        if ((fp2 = fopen(Form("%s/helHAP_%d.appcor.dat", OUTDIR, nrun), "w")) == NULL) {
            fprintf(stderr, "Can not open %s/helHAP_%d.appcor.dat", OUTDIR, nrun);
            exit(-1);
        }
    }

    fscanf(fp1, "%d", &tempi[0]);
    fprintf(fp2, "%d\n", gNRing);
    for (Int_t i = 0; i < gNRing; i++) {
        fscanf(fp1, "%d%d%d%d%x%d", &tempi[0], &tempi[1], &tempi[2], &tempi[3], &tempi[4], &tempi[5]);
        fgets(tempc, 300, fp1);
        tempc[strlen(tempc) - 1] = '\0';
        fprintf(fp2, "%d\t%d\t%d\t%d\t%08x\t%d\t%s\n", tempi[0], tempi[1], tempi[2], tempi[3], tempi[4], gErrorRing[i], tempc);
    }

    fclose(fp1);
    fclose(fp2);

    if (select == 1)
        gSystem->Rename(Form("%s/helRIN_%d.appcor.dat", OUTDIR, nrun), Form("%s/helRIN_%d.dat", OUTDIR, nrun));
    else if (select == 2)
        gSystem->Rename(Form("%s/helHAP_%d.appcor.dat", OUTDIR, nrun), Form("%s/helHAP_%d.dat", OUTDIR, nrun));

    return 0;
}

void usage(int argc, char** argv)
{
    printf("usage: %s [options] RUN_NUMBER\n", argv[0]);
    printf("  -c, --cfgfile=config.cfg   Set configuration file name\n");
    printf("  -h, --help                 This small usage guide\n");
    printf("  -i, --indir=.              Set input directory\n");
    printf("  -o, --outdir=.             Set output directory\n");
}

