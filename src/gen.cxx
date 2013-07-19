#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <libconfig.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "THaEvent.h"

#include "hel.h"

#define NDATA 16

struct datainfo {
    int index;
    char name[300];
};

struct treeinfo {
    char name[300];
    char prefix[300];
    datainfo data[NDATA];
};

treeinfo tree_HEL, tree_RIN, tree_HAP;

FILE *fp1, *fp2, *fp3;

Int_t inserttir(Int_t nrun, Int_t ntir);
Int_t insertring(Int_t nrun, Int_t nring, Int_t select);
void usage(int argc, char** argv);

#include "isexist.h"

Char_t CFGFILE[300] = "./config.cfg";
Char_t INFODIR[300] = ".";
Char_t ROOTDIR[300] = ".";

int main(int argc, char** argv) {
    int c;

    while (1) {
        static struct option long_options[] = {
            { "help", no_argument, 0, 'h'},
            { "cfgfile", required_argument, 0, 'c'},
            { "infodir", required_argument, 0, 'i'},
            { "rootdir", required_argument, 0, 'r'},
            { 0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "c:hi:r:", long_options, &option_index);

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
            strcpy(INFODIR, optarg);
            break;
        case 'r':
            strcpy(ROOTDIR, optarg);
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
    config_setting_t *dataelem;

    config_init(&cfg);

    if (!config_read_file(&cfg, CFGFILE)) {
        fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg), config_error_line(&cfg), config_error_text(&cfg));
        config_destroy(&cfg);
        exit(EXIT_FAILURE);
    }

    Bool_t configerror = kFALSE;
    const char *temp;

    Int_t ntir = 0;
    setting = config_lookup(&cfg, "tirinfo");
    if (setting != NULL) {
        config_lookup_string(&cfg, "tirinfo.name", &temp);
        strcpy(tree_HEL.name, temp);
        config_lookup_string(&cfg, "tirinfo.prefix", &temp);
        strcpy(tree_HEL.prefix, temp);
        setting = config_lookup(&cfg, "tirinfo.data");
        ntir = config_setting_length(setting);
        for (int i = 0; i < ntir; i++) {
            dataelem = config_setting_get_elem(setting, i);
            config_setting_lookup_string(dataelem, "name", &temp);
            strcpy(tree_HEL.data[i].name, temp);
            config_setting_lookup_int(dataelem, "index", &tree_HEL.data[i].index);
        }
    }
    else
        configerror = kTRUE;

    setting = config_lookup(&cfg, "ringinfo");
    if (setting != NULL) {
        config_lookup_string(&cfg, "ringinfo.name", &temp);
        strcpy(tree_RIN.name, temp);
        config_lookup_string(&cfg, "ringinfo.prefix", &temp);
        strcpy(tree_RIN.prefix, temp);
        setting = config_lookup(&cfg, "ringinfo.data");
        NRING = config_setting_length(setting);
        for (int i = 0; i < NRING; i++) {
            strcpy(tree_RIN.data[i].name, config_setting_get_string_elem(setting, i));
            tree_RIN.data[i].index = -1;
        }
    }
    else
        configerror = kTRUE;

    setting = config_lookup(&cfg, "happexinfo");
    if (setting != NULL) {
        USEHAPPEX = kTRUE;
        config_lookup_string(&cfg, "happexinfo.name", &temp);
        strcpy(tree_HAP.name, temp);
        config_lookup_string(&cfg, "happexinfo.prefix", &temp);
        strcpy(tree_HAP.prefix, temp);
        setting = config_lookup(&cfg, "happexinfo.data");
        NHAPPEX = config_setting_length(setting);
        for (int i = 0; i < NHAPPEX; i++) {
            strcpy(tree_HAP.data[i].name, config_setting_get_string_elem(setting, i));
            tree_HAP.data[i].index = -1;
        }
    }
    else {
        USEHAPPEX = kFALSE;
    }

    if (configerror) {
        fprintf(stderr, "Invalid cfg file\n");
        exit(-1);
    }

    inserttir(nrun, ntir);
    insertring(nrun, NRING, 1);
    if (USEHAPPEX) insertring(nrun, NHAPPEX, 2);

    return 0;
}

Int_t inserttir(Int_t nrun, Int_t ntir) {
    printf("Opening existed rootfile ...\n");

    Int_t filecount = 0;
    Char_t filename[300];

    sprintf(filename, "%s/g2p_%d.root", ROOTDIR, nrun);

    while (isexist(filename)) {
        TFile *f = new TFile(filename, "UPDATE");

        if ((fp1 = fopen(Form("%s/hel_%d.dat", INFODIR, nrun), "r")) == NULL) {
            fprintf(stderr, "Can not open %s/hel_%d.dat", INFODIR, nrun);
            exit(-1);
        }

        TTree *t = (TTree *) f->Get("T");

        THaEvent *event = new THaEvent();

        t->SetBranchAddress("Event_Branch", &event);

        Int_t fHelicity_rep = 0, fHelicity_act = 0;
        Int_t fQRT = 0, fMPS = 0, fPairSync = 0;
        Int_t fTimeStamp = 0, fSeed = 0, fError = 0;
        Int_t fIRing = 0, fIHappex = 0;
        Int_t fDATA[NDATA];
        Int_t fEvNum = -100;

        TList newBranch;

        newBranch.Add(t->Branch(Form("%shel_act", tree_HEL.prefix), &fHelicity_act, "hel_act/I"));
        newBranch.Add(t->Branch(Form("%shel_rep", tree_HEL.prefix), &fHelicity_rep, "hel_rep/I"));
        newBranch.Add(t->Branch(Form("%sqrt", tree_HEL.prefix), &fQRT, "qrt/I"));
        newBranch.Add(t->Branch(Form("%smps", tree_HEL.prefix), &fMPS, "mps/I"));
        newBranch.Add(t->Branch(Form("%spairsync", tree_HEL.prefix), &fPairSync, "pairsync/I"));
        newBranch.Add(t->Branch(Form("%stimestamp", tree_HEL.prefix), &fTimeStamp, "timestamp/I"));
        newBranch.Add(t->Branch(Form("%sseed", tree_HEL.prefix), &fSeed, "seed/I"));
        newBranch.Add(t->Branch(Form("%serror", tree_HEL.prefix), &fError, "error/I"));
        newBranch.Add(t->Branch(Form("%snring", tree_HEL.prefix), &fIRing, "nring/I"));
        newBranch.Add(t->Branch(Form("%snhappex", tree_HEL.prefix), &fIHappex, "nhappex/I"));
        for (Int_t i = 0; i < ntir; i++)
            newBranch.Add(t->Branch(Form("%s%s", tree_HEL.prefix, tree_HEL.data[i].name), &fDATA[tree_HEL.data[i].index], Form("%s/I", tree_HEL.data[i].name)));

        Int_t nentries;
        Int_t gEvNum = 0, N;

        fscanf(fp1, "%d", &N);
        nentries = t->GetEntries();
        for (Int_t k = 0; k < nentries; k++) {
            t->GetEntry(k);
            gEvNum = Int_t(event->GetHeader()->GetEvtNum());
            if (gEvNum % 1000 == 0) printf("%d\n", gEvNum);
            while ((fEvNum < gEvNum) && (!feof(fp1))) {
                fscanf(fp1, "%d%d%d%d%d%d%d%x%d%d%d", &fEvNum, &fHelicity_act, &fHelicity_rep, &fQRT, &fPairSync, &fMPS, &fTimeStamp, &fSeed, &fError, &fIRing, &fIHappex);
                for (Int_t l = 0; l < NRING + NHAPPEX; l++) {
                    fscanf(fp1, "%d", &fDATA[l]);
                }
            }
            if (fEvNum != gEvNum) {
                fHelicity_rep = 0;
                fHelicity_act = 0;
                fQRT = 0;
                fPairSync = 0;
                fMPS = 0;
                fTimeStamp = 0;
                fSeed = 0;
                fError = 0x8000;
                fIRing = 0;
                fIHappex = 0;
                for (Int_t l = 0; l < NRING + NHAPPEX; l++) {
                    fDATA[l] = 0;
                }
            }
            TIter next(&newBranch);
            while (TBranch * workBranch = (TBranch*) next()) {
                workBranch->Fill();
            }
        }

        t->Write("", TObject::kOverwrite);

        f->Close();

        filecount++;
        sprintf(filename, "%s/g2p_%d_%d.root", ROOTDIR, nrun, filecount);

        fclose(fp1);
    }

    if (filecount == 0) return -1;

    return 0;
}

Int_t insertring(Int_t nrun, Int_t nring, Int_t select) {
    printf("Opening existed rootfile ...\n");

    Int_t filecount = 0;
    Char_t filename[300];
    treeinfo *ringtree = NULL;

    sprintf(filename, "%s/g2p_%d.root", ROOTDIR, nrun);
    while (isexist(filename)) {
        TFile *f = new TFile(filename, "UPDATE");
        if (select == 1) {
            if ((fp1 = fopen(Form("%s/helRIN_%d.dat", INFODIR, nrun), "r")) == NULL) {
                fprintf(stderr, "Can not open %s/helRIN_%d.dat", INFODIR, nrun);
                exit(-1);
            }
            ringtree = &tree_RIN;
        }
        else if (select == 2) {
            if ((fp1 = fopen(Form("%s/helHAP_%d.dat", INFODIR, nrun), "r")) == NULL) {
                fprintf(stderr, "Can not open %s/helHAP_%d.dat", INFODIR, nrun);
                exit(-1);
            }
            ringtree = &tree_HAP;
        }

        Int_t fHelicity_rep = 0, fHelicity_act = 0, fQRT = 0;
        Int_t fSeed = 0, fError = 0;
        Int_t fDATA[NDATA];
        Int_t fEvNum;

        TTree *t = new TTree(ringtree->name, ringtree->name);

        t->Branch(Form("%sevnum", ringtree->prefix), &fEvNum, "evnum/I");
        t->Branch(Form("%shel_act", ringtree->prefix), &fHelicity_act, "hel_act/I");
        t->Branch(Form("%shel_rep", ringtree->prefix), &fHelicity_rep, "hel_rep/I");
        t->Branch(Form("%sqrt", ringtree->prefix), &fQRT, "qrt/I");
        t->Branch(Form("%sseed", ringtree->prefix), &fSeed, "seed/I");
        t->Branch(Form("%serror", ringtree->prefix), &fError, "error/I");
        for (Int_t i = 0; i < nring; i++)
            t->Branch(Form("%s%s", ringtree->prefix, ringtree->data[i].name), &fDATA[i], Form("%s/I", ringtree->data[i].name));

        Int_t nentries;
        Int_t N, gEvNumMax = 0, gEvNumMin = 0;

        TTree *ori = (TTree *) f->Get("T");

        THaEvent *event = new THaEvent();
        ori->SetBranchAddress("Event_Branch", &event);

        nentries = ori->GetEntries();

        ori->GetEntry(0);
        gEvNumMin = Int_t(event->GetHeader()->GetEvtNum());
        ori->GetEntry(nentries - 10);
        gEvNumMax = Int_t(event->GetHeader()->GetEvtNum());
        for (Int_t i = 9; i >= 1; i--) {
            ori->GetEntry(nentries - i);
            Int_t temp = Int_t(event->GetHeader()->GetEvtNum());
            if (temp == gEvNumMax + 1) gEvNumMax = temp;
        }

        ori = NULL;

        fscanf(fp1, "%d", &N);
        for (Int_t k = 0; k < N; k++) {
            if (k % 1000 == 0) printf("%d\n", k);
            fscanf(fp1, "%d%d%d%d%x%d", &fEvNum, &fHelicity_act, &fHelicity_rep, &fQRT, &fSeed, &fError);
            for (Int_t l = 0; l < nring; l++)
                fscanf(fp1, "%d", &fDATA[l]);
            if ((fEvNum >= gEvNumMin) && (fEvNum <= gEvNumMax)) t->Fill();
            if (fEvNum > gEvNumMax) break;
        }

        f->Write("", TObject::kOverwrite);

        f->Close();

        filecount++;
        sprintf(filename, "%s/g2p_%d_%d.root", ROOTDIR, nrun, filecount);

        fclose(fp1);
    }

    if (filecount == 0) return -1;

    return 0;
}

void usage(int argc, char** argv) {
    printf("usage: %s [options] RUN_NUMBER\n", argv[0]);
    printf("  -c, --cfgfile=config.cfg   Set configuration file name\n");
    printf("  -h, --help                 This small usage guide\n");
    printf("  -i, --infodir=.            Set helicity info directory\n");
    printf("  -r, --rootdir=.            Set rootfile directory\n");
}
