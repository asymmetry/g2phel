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

FILE *fp1, *fp2;

//Global variables
Int_t gHelicity_act[LENRIN];
Int_t gSeed_rep[LENRIN];
Int_t gDATA[NDATA][LENRIN];
Int_t gError[LENRIN];
Int_t gUsed[LENRIN];
Int_t gN;

Int_t readin(Int_t nrun, Int_t nring, Int_t select);
Int_t align(Int_t nrun, Int_t nring, Int_t select);
Int_t printout(Int_t nrun, Int_t nring, Int_t select);
void usage(int argc, char** argv);

Char_t CFGFILE[300] = "./config.cfg";
Char_t INDIR[300] = ".";
Char_t OUTDIR[300] = ".";

int main(int argc, char** argv) {
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
    printout(nrun, NRING, 1);

    if (USEHAPPEX) {
        gSystem->Rename(Form("%s/hel_%d.dat", OUTDIR, nrun), Form("%s/helTIR_%d.nohapp.dat", INDIR, nrun));
        readin(nrun, NHAPPEX, 2);
        align(nrun, NHAPPEX, 2);
        printout(nrun, NHAPPEX, 2);
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

    Int_t fHelRing_rep, fQRTRing, temp1;

    for (Int_t i = 0; i < LENRIN; i++) {
        gHelicity_act[i] = 0;
        gSeed_rep[i] = 0;
        gError[i] = 0;
        for (Int_t k = 0; k < NDATA; k++)
            gDATA[k][i] = 0;
        gUsed[i] = 0;
    }

    fscanf(fp1, "%d", &gN);
    for (Int_t i = 0; i < gN; i++) {
        fscanf(fp1, "%d%d%d%d%x%d", &temp1, &gHelicity_act[i], &fHelRing_rep, &fQRTRing, &gSeed_rep[i], &gError[i]);
        for (Int_t k = 0; k < nring; k++)
            fscanf(fp1, "%d", &gDATA[k][i]);
    }

    fclose(fp1);

    return 0;
}

Int_t align(Int_t nrun, Int_t nring, Int_t select) {
    printf("Constructing new TIR helicity information ...\n");

    if (select == 1) {
        if ((fp1 = fopen(Form("%s/helTIR_%d.noring.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helTIR_%d.noring.dat", INDIR, nrun);
            exit(-1);
        }
    }
    else if (select == 2) {
        if ((fp1 = fopen(Form("%s/helTIR_%d.nohapp.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helTIR_%d.nohapp.dat", INDIR, nrun);
            exit(-1);
        }
    }
    if ((fp2 = fopen(Form("%s/hel_%d.dat", OUTDIR, nrun), "w")) == NULL) {
        fprintf(stderr, "Can not open %s/hel_%d.dat", OUTDIR, nrun);
        exit(-1);
    }

    Int_t fHelicity_rep[NBUFFER], fHelicity_act[NBUFFER], fQRT[NBUFFER], fMPS[NBUFFER], fPairSync[NBUFFER];
    Int_t fTimeStamp[NBUFFER], fSeedTIR_rep[NBUFFER], fError[NBUFFER];
    Int_t fEvNum[NBUFFER];
    Char_t fCharTemp[NBUFFER][300];
    Char_t tempc[300];
    Int_t fDataP[NDATA], fDataN[NDATA];
    Int_t tempi[20];
    Int_t N = 0;
    Int_t NBuff = 0;
    Int_t fLastSeed = 0;
    Int_t Index = 0, IStart = 0;

    fscanf(fp1, "%d", &N);
    for (Int_t k = 0; k < N; k++) {
        fscanf(fp1, "%d%d%d%d%d%d%d%x%d", &tempi[0], &tempi[1], &tempi[2], &tempi[3], &tempi[4], &tempi[5], &tempi[6], &tempi[7], &tempi[8]);
        fgets(tempc, 300, fp1);
        tempc[strlen(tempc) - 1] = '\0';

        if (tempi[8] == 0) {
            if ((tempi[7] != fLastSeed)&&(NBuff > 0)) { // found a pattern in TIR
                Int_t Filled = 0;
                Index = IStart;
                while ((gSeed_rep[Index] != fLastSeed)&&(Index < gN)) Index++; // search this pattern in ring buffer
                if (Index < gN) { // found this pattern in ring buffer
                    Int_t NPattern = 0;
                    IStart = Index;
                    for (Int_t i = 0; i < NDATA; i++) {
                        fDataP[i] = 0;
                        fDataN[i] = 0;
                    }
                    while ((gSeed_rep[Index] == fLastSeed)&&(Index < gN)) {
                        if (gHelicity_act[Index] == 1) {
                            for (Int_t j = 0; j < nring; j++) fDataP[j] += gDATA[j][Index];
                        }
                        else if (gHelicity_act[Index] == -1) {
                            for (Int_t j = 0; j < nring; j++) fDataN[j] += gDATA[j][Index];
                        }
                        else {
                            NPattern = 10;
                            break;
                        }
                        Index++;
                        NPattern++;
                    }
                    if (NPattern == 4) {
                        Int_t GoodPattern = 0;
                        for (Int_t i = 0; i < NBuff; i++) {
                            if (fHelicity_act[i] == 1) GoodPattern |= 0x01;
                            if (fHelicity_act[i] == -1) GoodPattern |= 0x10;
                        }
                        if (GoodPattern == 0x11) {
                            for (Int_t i = 0; i < NBuff; i++) {
                                fprintf(fp2, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%08x\t%d\t%s\t", fEvNum[i], fHelicity_act[i], fHelicity_rep[i], fQRT[i], fPairSync[i], fMPS[i], fTimeStamp[i], fSeedTIR_rep[i], fError[i], fCharTemp[i]);
                                if (((Filled & 0x01) == 0)&&(fHelicity_act[i] == 1)) {
                                    if (nring > 0) {
                                        for (Int_t j = 0; j < nring - 1; j++) fprintf(fp2, "%d\t", fDataP[j]);
                                        fprintf(fp2, "%d\n", fDataP[nring - 1]);
                                    }
                                    Filled |= 0x01;
                                }
                                else if (((Filled & 0x10) == 0)&&(fHelicity_act[i] == -1)) {
                                    if (nring > 0) {
                                        for (Int_t j = 0; j < nring - 1; j++) fprintf(fp2, "%d\t", fDataN[j]);
                                        fprintf(fp2, "%d\n", fDataN[nring - 1]);
                                    }
                                    Filled |= 0x10;
                                }
                                else {
                                    if (nring > 0) {
                                        for (Int_t j = 0; j < nring - 1; j++) fprintf(fp2, "0\t");
                                        fprintf(fp2, "0\n");
                                    }
                                }
                            }
                            for (Int_t i = IStart; i < Index; i++) gUsed[i] = 1;
                        }
                    }
                }
                if (Filled == 0) {
                    for (Int_t i = 0; i < NBuff; i++) {
                        fprintf(fp2, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%08x\t%d\t%s\t", fEvNum[i], fHelicity_act[i], fHelicity_rep[i], fQRT[i], fPairSync[i], fMPS[i], fTimeStamp[i], fSeedTIR_rep[i], 0x1000 * select, fCharTemp[i]);
                        if (nring > 0) {
                            for (Int_t j = 0; j < nring - 1; j++) fprintf(fp2, "0\t");
                            fprintf(fp2, "0\n");
                        }
                    }
                }
                NBuff = 0;
            }

            // possible new pattern, store info into buff
            fEvNum[NBuff] = tempi[0];
            fHelicity_act[NBuff] = tempi[1];
            fHelicity_rep[NBuff] = tempi[2];
            fQRT[NBuff] = tempi[3];
            fPairSync[NBuff] = tempi[4];
            fMPS[NBuff] = tempi[5];
            fTimeStamp[NBuff] = tempi[6];
            fSeedTIR_rep[NBuff] = tempi[7];
            fError[NBuff] = tempi[8];
            strcpy(fCharTemp[NBuff], tempc);
            fLastSeed = tempi[7];
            NBuff++;
        }
        else {
            fprintf(fp2, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%08x\t%d\t%s\t", tempi[0], tempi[1], tempi[2], tempi[3], tempi[4], tempi[5], tempi[6], tempi[7], tempi[8], tempc);
            if (nring > 0) {
                for (Int_t j = 0; j < nring - 1; j++) fprintf(fp2, "0\t");
                fprintf(fp2, "0\n");
            }
        }
    }

    if (NBuff != 0) {
        for (Int_t i = 0; i < NBuff; i++) {
            fprintf(fp2, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%08x\t%d\t%s\t", fEvNum[i], fHelicity_act[i], fHelicity_rep[i], fQRT[i], fPairSync[i], fMPS[i], fTimeStamp[i], fSeedTIR_rep[i], 0x1000 * select, fCharTemp[i]);
            if (nring > 0) {
                for (Int_t j = 0; j < nring - 1; j++) fprintf(fp2, "0\t");
                fprintf(fp2, "0\n");
            }
        }
        NBuff = 0;
    }

    fclose(fp1);
    fclose(fp2);

    gSystem->Exec(Form("sort -n \"%s/hel_%d.dat\" >> temp.dat", OUTDIR, nrun));
    gSystem->Exec(Form("echo \"%d\" > \"%s/hel_%d.dat\"", N, OUTDIR, nrun));
    gSystem->Exec(Form("cat temp.dat >> \"%s/hel_%d.dat\"", OUTDIR, nrun));
    gSystem->Exec("rm temp.dat");

    return 0;
}

Int_t printout(Int_t nrun, Int_t nring, Int_t select) {
    if (select == 1) {
        if ((fp1 = fopen(Form("%s/helRIN_%d.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helRIN_%d.dat", INDIR, nrun);
            exit(-1);
        }
        if ((fp2 = fopen(Form("%s/helRIN_%d.appcor.dat", OUTDIR, nrun), "w")) == NULL) {
            fprintf(stderr, "Can not open %s/helRIN_%d.appcor.dat", OUTDIR, nrun);
            exit(-1);
        }
    }
    else if (select == 2) {
        if ((fp1 = fopen(Form("%s/helHAP_%d.dat", INDIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/helHAP_%d.dat", INDIR, nrun);
            exit(-1);
        }
        if ((fp2 = fopen(Form("%s/helHAP_%d.appcor.dat", OUTDIR, nrun), "w")) == NULL) {
            fprintf(stderr, "Can not open %s/helHAP_%d.appcor.dat", OUTDIR, nrun);
            exit(-1);
        }
    }

    Int_t temp[10];

    fscanf(fp1, "%d", &temp[0]);
    fprintf(fp2, "%d\n", gN);
    for (Int_t i = 0; i < gN; i++) {
        fscanf(fp1, "%d%d%d%d%x%d", &temp[0], &temp[1], &temp[2], &temp[3], &temp[4], &temp[5]);
        if (gUsed[i] == 1)
            fprintf(fp2, "%d\t%d\t%d\t%d\t%08x\t%d", temp[0], temp[1], temp[2], temp[3], temp[4], gError[i]);
        else
            fprintf(fp2, "%d\t%d\t%d\t%d\t%08x\t%d", temp[0], temp[1], temp[2], temp[3], temp[4], 0x80);
        for (Int_t k = 0; k < nring; k++) {
            fscanf(fp1, "%d", &temp[6]);
            fprintf(fp2, "\t%d", temp[6]);
        }
        fprintf(fp2, "\n");
    }

    fclose(fp1);
    fclose(fp2);

    if (select == 1)
        gSystem->Rename(Form("%s/helRIN_%d.appcor.dat", OUTDIR, nrun), Form("%s/helRIN_%d.dat", OUTDIR, nrun));
    else if (select == 2)
        gSystem->Rename(Form("%s/helHAP_%d.appcor.dat", OUTDIR, nrun), Form("%s/helHAP_%d.dat", OUTDIR, nrun));

    return 0;
}

void usage(int argc, char** argv) {
    printf("usage: %s [options] RUN_NUMBER\n", argv[0]);
    printf("  -c, --cfgfile=config.cfg   Set configuration file name\n");
    printf("  -h, --help                 This small usage guide\n");
    printf("  -i, --indir=.              Set input directory\n");
    printf("  -o, --outdir=.             Set output directory\n");
}

