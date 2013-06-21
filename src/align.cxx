#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <libconfig.h>

#include "TROOT.h"
#include "TSystem.h"

#include "hel.h"

#define NDATA 8

FILE *fp1, *fp2;

//Global variables
Int_t gHelicity_act[LEN];
Int_t gSeed_rep[LEN];
Int_t gPhase[LEN];
Int_t gDATA[NDATA][LEN];
Int_t gError[LEN];
Int_t gN;

Int_t readin(Int_t nrun, Int_t nring, Int_t select);
Int_t align(Int_t nrun, Int_t nring, Int_t select);
void usage(int argc, char** argv);

Char_t CFGFILE[300] = "./config.cfg";
Char_t INDIR[300] = ".";
Char_t OUTDIR[300] = ".";

int main(int argc, char** argv) {
    int c;

    while (1) {
        static struct option long_options[] = {
                { "help", no_argument, 0, 'h' },
                { "cfgfile", required_argument, 0, 'c' },
                { "indir", required_argument, 0, 'i' },
                { "outdir", required_argument, 0, 'o' },
                { 0, 0, 0, 0 } };

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
    }
    else {
        usage(argc, argv);
        exit(-1);
    }

    config_t cfg;
    config_setting_t *setting;

    config_init(&cfg);

    if (!config_read_file(&cfg, CFGFILE)) {
        fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
                config_error_line(&cfg), config_error_text(&cfg));
        config_destroy(&cfg);
        exit(EXIT_FAILURE);
    }

    Bool_t configerror = kFALSE;

    setting = config_lookup(&cfg, "ringinfo.data");
    if (setting != NULL) {
        NRING = config_setting_length(setting);
    }
    else
        configerror = kTRUE;

    setting = config_lookup(&cfg, "happexinfo.data");
    if (setting != NULL) {
        USEHAPPEX = kTRUE;
        NHAPPEX = config_setting_length(setting);
    }
    else {
        USEHAPPEX = kFALSE;
    }

    if (configerror) {
        fprintf(stderr, "Invalid cfg file\n");
        exit(-1);
    }

    readin(nrun, NRING, 1);
    align(nrun, NRING, 1);

    if (USEHAPPEX) {
        gSystem->Exec(
                Form("mv -vf %s/hel_%d.dat %s/helTIR_%d.nohapp.dat", OUTDIR,
                        nrun, INDIR, nrun));
        readin(nrun, NHAPPEX, 2);
        align(nrun, NHAPPEX, 2);
    }

    return 0;
}

Int_t readin(Int_t nrun, Int_t nring, Int_t select) {
    if (select == 1) {
        if ((fp1 = fopen(Form("%s/helRIN_%d.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helRIN_%d.dat", INDIR, nrun);
            exit(-1);
        }
    }
    else if (select == 2) {
        if ((fp1 = fopen(Form("%s/helHAP_%d.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helHAP_%d.dat", INDIR, nrun);
            exit(-1);
        }
    }

    Int_t fHelRing_rep, fQRTRing, fPhase = 0, temp1;

    for (Int_t i = 0; i < LEN; i++) {
        gHelicity_act[i] = 0;
        gSeed_rep[i] = 0;
        gError[i] = 0;
        for (Int_t k = 0; k < NDATA; k++)
            gDATA[k][i] = 0;
        gPhase[i] = 0;
    }

    fscanf(fp1, "%d", &gN);
    for (Int_t i = 0; i < gN; i++) {
        fscanf(fp1, "%d%d%d%d%x%d", &temp1, &gHelicity_act[i], &fHelRing_rep,
                &fQRTRing, &gSeed_rep[i], &gError[i]);
        if (fQRTRing == 1)
            fPhase = 0;
        else
            fPhase++;
        if (gError[i] > 0) fPhase = 0;
        gPhase[i] = fPhase;
        for (Int_t k = 0; k < nring; k++)
            fscanf(fp1, "%d", &gDATA[k][i]);
    }

    fclose(fp1);

    return 0;
}

Int_t align(Int_t nrun, Int_t nring, Int_t select) {
    printf("Constructing new TIR helicity information ...\n");

    if (select == 1) {
        if ((fp1 = fopen(Form("%s/helTIR_%d.noring.dat", INDIR, nrun), "r"))
                == NULL) {
            fprintf(stderr, "Can not open %s/helTIR_%d.noring.dat", INDIR,
                    nrun);
            exit(-1);
        }
    }
    else if (select == 2) {
        if ((fp1 = fopen(Form("%s/helTIR_%d.nohapp.dat", INDIR, nrun), "r"))
                == NULL) {
            fprintf(stderr, "Can not open %s/helTIR_%d.nohapp.dat", INDIR,
                    nrun);
            exit(-1);
        }
    }
    if ((fp2 = fopen(Form("%s/hel_%d.dat", OUTDIR, nrun), "w")) == NULL) {
        fprintf(stderr, "Can not open %s/hel_%d.dat", OUTDIR, nrun);
        exit(-1);
    }

    Int_t fHelicity_rep = 0, fHelicity_act = 0, fQRT = 0, fMPS = 0, fPairSync =
            0, fTimeStamp = 0, fError = 0;
    Int_t fEvNum;
    Char_t temp[300];
    Int_t fSeedTIR_rep = 0, fPhaseTIR = 0, fPolarityTIR_rep = 0;
    Int_t pLastp = 0, pLastm = 0;
    Int_t fDATA[8];
    Bool_t fNewFlag = kTRUE;
    Int_t N = 0;

    fscanf(fp1, "%d", &N);
    fprintf(fp2, "%d\n", N);
    for (Int_t k = 0; k < N; k++) {
        fscanf(fp1, "%d%d%d%d%d%d%d%x%d", &fEvNum, &fHelicity_act,
                &fHelicity_rep, &fQRT, &fPairSync, &fMPS, &fTimeStamp,
                &fSeedTIR_rep, &fError);
        fgets(temp, 300, fp1);
        temp[strlen(temp) - 1] = '\0';
        for (Int_t j = 0; j < NDATA; j++)
            fDATA[j] = 0;
        if (fError == 0) {
            fPolarityTIR_rep = fSeedTIR_rep & 0x01;
            if (fQRT == 1) {
                fPhaseTIR = 0;
            }
            else if (fHelicity_rep == fPolarityTIR_rep) {
                fPhaseTIR = 3;
            }
            else if (fPairSync == 1) {
                fPhaseTIR = 2;
            }
            else {
                fPhaseTIR = 1;
            }
            if (fNewFlag) {
                for (Int_t i = (pLastp < pLastm) ? pLastp : pLastm + 1; i < gN;
                        i++) {
                    if ((gPhase[i] == fPhaseTIR)
                            && (gSeed_rep[i] == fSeedTIR_rep)) {
                        pLastp = i - 1;
                        pLastm = i - 1;
                        fNewFlag = kFALSE;
                        break;
                    }
                }
                if (!((gPhase[pLastp + 1] == fPhaseTIR)
                        && (gSeed_rep[pLastp + 1] == fSeedTIR_rep))) {
                    fError = fError | 0x20;
                    fNewFlag = kTRUE;
                    fprintf(fp2, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%08x\t%d\t%s\t",
                            fEvNum, fHelicity_act, fHelicity_rep, fQRT,
                            fPairSync, fMPS, fTimeStamp, fSeedTIR_rep, fError,
                            temp);
                    if (nring > 0) {
                        for (Int_t j = 0; j < nring - 1; j++)
                            fprintf(fp2, "%d\t", fDATA[j]);
                        fprintf(fp2, "%d\n", fDATA[nring - 1]);
                    }
                    continue;
                }
            }
            if (fHelicity_act == 1) {
                if (!((gPhase[pLastp] == fPhaseTIR)
                        && (gSeed_rep[pLastp] == fSeedTIR_rep))) {
                    for (Int_t i = pLastp + 1; (i < pLastp + 1000) && (i < gN);
                            i++) {
                        if (gHelicity_act[i] == 1) {
                            for (Int_t j = 0; j < nring; j++)
                                fDATA[j] += gDATA[j][i];
                        }
                        if ((gPhase[i] == fPhaseTIR)
                                && (gSeed_rep[i] == fSeedTIR_rep)) {
                            pLastp = i;
                            break;
                        }
                    }
                    if ((gPhase[pLastp] == fPhaseTIR)
                            && (gSeed_rep[pLastp] == fSeedTIR_rep)) {
                    }
                    else {
                        fError = fError | 0x10;
                        fNewFlag = kTRUE;
                    }
                }
            }
            else {
                if (!((gPhase[pLastm] == fPhaseTIR)
                        && (gSeed_rep[pLastm] == fSeedTIR_rep))) {
                    for (Int_t i = pLastm + 1; (i < pLastm + 1000) && (i < gN);
                            i++) {
                        if (gHelicity_act[i] == -1) {
                            for (Int_t j = 0; j < nring; j++)
                                fDATA[j] += gDATA[j][i];
                        }
                        if ((gPhase[i] == fPhaseTIR)
                                && (gSeed_rep[i] == fSeedTIR_rep)) {
                            pLastm = i;
                            break;
                        }
                    }
                    if ((gPhase[pLastm] == fPhaseTIR)
                            && (gSeed_rep[pLastm] == fSeedTIR_rep)) {
                    }
                    else {
                        fError = fError | 0x10;
                        fNewFlag = kTRUE;
                    }
                }
            }
        }
        else if (fError == 8) {
        }
        else {
            fNewFlag = kTRUE;
        }

        fprintf(fp2, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%08x\t%d\t%s\t", fEvNum,
                fHelicity_act, fHelicity_rep, fQRT, fPairSync, fMPS, fTimeStamp,
                fSeedTIR_rep, fError, temp);
        if (nring > 0) {
            for (Int_t j = 0; j < nring - 1; j++)
                fprintf(fp2, "%d\t", fDATA[j]);
            fprintf(fp2, "%d\n", fDATA[nring - 1]);
        }
    }

    fclose(fp1);
    fclose(fp2);

    return 0;
}

void usage(int argc, char** argv) {
    printf("usage: %s [options] RUN_NUMBER\n", argv[0]);
    printf("  -c, --cfgfile=config.cfg   Set configuration file name\n");
    printf("  -h, --help                 This small usage guide\n");
    printf("  -i, --indir=.              Set input directory\n");
    printf("  -o, --outdir=.             Set output directory\n");
}

