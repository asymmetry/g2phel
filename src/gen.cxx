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

using namespace std;

#define NDATA 16

struct datainfo
{
    int index;
    char name[300];
};

struct treeinfo
{
    char name[300];
    char prefix[300];
    datainfo data[NDATA];
};

treeinfo tree_HEL, tree_RIN, tree_HAP;

FILE *fp1,*fp2,*fp3;

Int_t inserttir(Int_t nrun, Int_t n);
Int_t construct(Int_t nrun);
void usage(int argc, char** argv);

#include "isexist.h"

Char_t CFGFILE[300] = "./config.cfg";
Char_t INFODIR[300] = ".";
Char_t ROOTDIR[300] = ".";

int main(int argc,char** argv)
{
    int c;

    while (1) {
        static struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"cfgfile", required_argument, 0, 'c'},
            {"infodir", required_argument, 0, 'i'},
            {"rootdir", required_argument, 0, 'r'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "c:hi:o:", long_options, &option_index);

        if (c==-1) break;
        
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

    if (optind<argc) {
        nrun = atoi(argv[optind++]);
    }
    else {
        usage(argc, argv);
        exit(-1);
    }

    config_t cfg;
    config_setting_t *setting;
    config_setting_t *datasetting;

    config_init(&cfg);

    if (!config_read_file(&cfg, CFGFILE)) {
        fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg), config_error_line(&cfg), config_error_text(&cfg));
        config_destroy(&cfg);
        exit(EXIT_FAILURE);
    }

    Bool_t configerror = kFALSE;

    Int_t n_data;
    setting = config_lookup(&cfg, "tirinfo");
    if (setting!=NULL) {
        strcpy(tree_HEL.name, config_setting_get_string(&setting, "name"));
        strcpy(tree_HEL.prefix, config_setting_get_string(&setting, "prefix"));
        datasetting = config_setting_get_member(&setting, "data");
        n_data = config_setting_length(datasetting);
        for (int i=0; i<n; i++) {
            strcpy(tree_HEL.data[i].name, config_setting_get_string_elem(&datasetting, i));
            tree_HEL.data[i].index = config_setting_get_int_elem(&datasetting, i);
        }
    }
    else configerror = kTRUE;
    
    setting = config_lookup(&cfg, "ringinfo.data");
    if (setting!=NULL) {
        strcpy(tree_RIN.name, config_setting_get_string(&setting, "name"));
        strcpy(tree_RIN.prefix, config_setting_get_string(&setting, "prefix"));
        datasetting = config_setting_get_member(&setting, "data");
        NRING = config_setting_length(setting);
        for (int i=0; i<NRING; i++) {
            strcpy(tree_RIN.data[i].name, config_setting_get_string_elem(&datasetting, i));
            tree_RIN.data[i].index = config_setting_get_int_elem(&datasetting, i);
        }
    }
    else configerror = kTRUE;

    setting = config_lookup(&cfg, "happexinfo.data");
    if (setting!=NULL) {
        USEHAPPEX = kTRUE;
        strcpy(tree_HAP.name, config_setting_get_string(&setting, "name"));
        strcpy(tree_HAP.prefix, config_setting_get_string(&setting, "prefix"));
        datasetting = config_setting_get_member(&setting, "data");
        NHAPPEX = config_setting_length(setting);
        for (int i=0; i<NHAPPEX; i++) {
            strcpy(tree_HAP.data[i].name, config_setting_get_string_elem(&datasetting, i));
            tree_HAP.data[i].index = config_setting_get_int_elem(&datasetting, i);
        }
    }
    else {
        USEHAPPEX = kFALSE;
    }
    
    if (configerror) {
        fprintf(stderr, "Invalid cfg file\n");
        exit(-1);
    }

    inserttir(nrun, n_data);
    insertring(nrun, NRING, 1);
    if (USEHAPPEX) insertring(nrun, NHAPPEX, 2);
    
    Int_t filecount=1;
    Char_t filename[300];

    // if(select==1){
    //     sprintf(filename,"g2p_%d.root",nrun);
    //     insert(nrun,filename);
    //     sprintf(filename,"g2p_%d_%d.root",nrun,filecount);
    //     while(isexist(filename)){
    //         insert(nrun,filename);
    //         filecount++;
    //         sprintf(filename,"g2p_%d_%d.root",nrun,filecount);
    //     }
    // }
    // else if(select==2)construct(nrun);
    
    return 0;
}

Int_t inserttir (Int_t nrun, Int_t n_data)
{
    printf("Opening existed rootfile ...\n");

    if ((fp1 = fopen(Form("%s/hel_%d.dat", INFODIR, nrun), "r"))==NULL){
        fprintf(stderr, "Can not open %s/hel_%d.dat", INFODIR, nrun);
        exit(-1);
    }
    
    Int_t filecount = 0;

    TFile *f = new TFile(Form("%s/g2p_%d.root", ROOTDIR, nrun), "UPDATE");
    while (f.IsOpen()) {
        TTree *t = (TTree *)f->Get("T");

        THaEvent *event = new THaEvent();

        t->SetBranchAddress("Event_Branch", &event);

        Int_t gEvNum;

        TBranch *bHelicity_act, *bHelicity_rep, *bQRT, *bMPS, *bPairSync;
        TBranch *bTimeStamp, *bSeed, *bError, *bIRing, *bIHappex;
        TBranch *bDATA[NDATA];

        Int_t fHelicity_rep = 0, fHelicity_act = 0;
        Int_t fQRT = 0, fMPS = 0, fPairSync=0;
        Int_t fTimeStamp = 0, fSeed = 0, fError = 0;
        Int_t fIRing = 0, fIHappex = 0;
        Int_t fDATA[NDATA];
        Int_t fEvNum;

        TList newBranch;

        newBranch.Add(bHelicity_act);
        newBranch.Add(bHelicity_rep);
        newBranch.Add(bQRT);
        newBranch.Add(bMPS);
        newBranch.Add(bPairSync);
        newBranch.Add(bTimeStamp);
        newBranch.Add(bSeed);
        newBranch.Add(bError);
        newBranch.Add(bIRing);
        newBranch.Add(bIHappex);
        for (Int_t i=0; i<n_data; i++) newBranch.Add(bDATA[i]);

        map<TBranch*,Int_t> indexmap;
        for (Int_t i=0; i<n_data; i++) branchmap[bDATA[i]] = tree_HEL.data[i].index;

        map<TBranch*,Int_t*> branchmap;
        branchmap[bHelicity_act] = &fHelicity_act;
        branchmap[bHelicity_rep] = &fHelicity_rep;
        branchmap[bQRT] = &fQRT;
        branchmap[bMPS] = &fMPS;
        branchmap[bPairSync] = &fPairSync;
        branchmap[bTimeStamp] = &fTimeStamp;
        branchmap[bSeed] = &fSeed;
        branchmap[bError] = &fError;
        branchmap[bIRing] = &fIRing;
        branchmap[bIHappex] = &fIHappex;
        for (Int_t i=0; i<n_data; i++) branchmap[bDATA[i]] = &fDATA[indexmap[bDATA[i]]];
        
        TIter next(&newBranch);
        TBranch* workBranch;
        while ((workBranch = (TBranch*)next())) {
            
        }
        
        

        
        
        filecount++;
        sprintf(filename,"g2p_%d_%d.root",nrun,filecount);
    }
    
    Int_t fEvNum=0;
    Int_t fT1=0,fT2=0;
    Int_t fL1A=0;
    Int_t N=0;
    Int_t temp[10];
    Int_t gEvNum=0,gEvNumMin=100000000,gEvNumMax=0;
    Int_t nentries;
    THaEvent *event=new THaEvent();

    t->SetBranchAddress("Event_Branch",&event);
    
    t1->Branch(Form("hel.%s.hel_rep",arm),&fHelicity_rep,"TIR Reported Helicity/I");
    t1->Branch(Form("hel.%s.hel_act",arm),&fHelicity_act,"TIR Actual Helicity/I");
    t1->Branch(Form("hel.%s.qrt",arm),&fQRT,"TIR QRT/I");
    t1->Branch(Form("hel.%s.mps",arm),&fMPS,"TIR MPS/I");
    t1->Branch(Form("hel.%s.pairsync",arm),&fPairSync,"TIR PairSync/I");
    t1->Branch(Form("hel.%s.timestamp",arm),&fTimeStamp,"TIR TimeStamp/I");
    t1->Branch(Form("hel.%s.seed",arm),&fSeed,"TIR Seed/I");   
    t1->Branch(Form("hel.%s.error",arm),&fError,"TIR Decode Error/I");
    t1->Branch(Form("hel.%s.numring",arm),&fIRing,"Event num in the ring/I");   
    t1->Branch(Form("hel.%s.numhappex",arm),&fIHappex,"Event num in HAPPEX/I");  
    t1->Branch(Form("hel.%s.bcmup",arm),&fBCMup,"Helicity Gated Upstream BCM/I");
    t1->Branch(Form("hel.%s.bcmdown",arm),&fBCMdown,"Helicity Gated Downstream BCM/I");

    if((arm[0]=='L')||(arm[0]=='R')){
        t1->Branch(Form("hel.%s.bcmuphappex",arm),&fBCMuph,"Helicity Gated Upstream BCM/I");
        t1->Branch(Form("hel.%s.bcmdownhappex",arm),&fBCMdownh,"Helicity Gated Downstream BCM/I");
    }
    t1->Branch(Form("hel.%s.time",arm),&fTime,"Helicity Gated Time/I");

    fscanf(fp1,"%d",&N);
    nentries=t->GetEntries();
    for(Int_t k=0;k<nentries;k++){
        t->GetEntry(k);
        gEvNum=Int_t(event->GetHeader()->GetEvtNum());
        printf("%d\n",gEvNum);
        if(gEvNum<gEvNumMin)gEvNumMin=gEvNum;
        if(gEvNum>gEvNumMax)gEvNumMax=gEvNum;
        do{
            if((arm[0]=='L')||(arm[0]=='R')){
                fscanf(fp1,"%d%d%d%d%d%d%d%x%d%d%d%d%d%d%d%d",&temp[0],&fHelicity_act,&fHelicity_rep,&fQRT,&fPairSync,&fMPS,&fTimeStamp,&fSeed,&fError,&fIRing,&fIHappex,&fBCMup,&fBCMdown,&fTime,&fBCMuph,&fBCMdownh);
            }
            else{
                fscanf(fp1,"%d%d%d%d%d%d%d%x%d%d%d%d%d%d",&temp[0],&fHelicity_act,&fHelicity_rep,&fQRT,&fPairSync,&fMPS,&fTimeStamp,&fSeed,&fError,&fIRing,&fIHappex,&fBCMup,&fBCMdown,&fTime);
            }
        }while((temp[0]!=gEvNum)&&(!feof(fp1)));
        if(temp[0]==gEvNum){
            t1->Fill();
        }
        else{
            fHelicity_rep=0;
            fHelicity_act=0;
            fQRT=0;
            fPairSync=0;
            fMPS=0;
            fTimeStamp=0;
            fSeed=0;
            fError=0;
            fIRing=0;
            fIHappex=0;
            fBCMup=0;
            fBCMdown=0;
            fTime=0;
            fBCMuph=0;
            fBCMdownh=0;
            t1->Fill();
        }
    }

    t1->Write();
    
    t2->Branch(Form("hel.%s.evnum",arm),&fEvNum,"Event Number/I");
    t2->Branch(Form("hel.%s.hel_rep",arm),&fHelicity_rep,"Ring Reported Helicity/I");
    t2->Branch(Form("hel.%s.hel_act",arm),&fHelicity_act,"Ring Actual Helicity/I");
    t2->Branch(Form("hel.%s.qrt",arm),&fQRT,"Ring QRT/I");
    t2->Branch(Form("hel.%s.time",arm),&fTime,"Ring Time/I");
    t2->Branch(Form("hel.%s.seed",arm),&fSeed,"Ring Seed/I");   
    t2->Branch(Form("hel.%s.error",arm),&fError,"Ring Decode Error/I");
    t2->Branch(Form("hel.%s.bcmup",arm),&fBCMup,"Helicity Gated Upstream BCM/I");
    t2->Branch(Form("hel.%s.bcmdown",arm),&fBCMdown,"Helicity Gated Downstream BCM/I");
    t2->Branch(Form("hel.%s.L1A",arm),&fL1A,"Helicity Gated L1A/I");

    if(arm[0]=='L'){
        t2->Branch(Form("hel.%s.T3",arm),&fT1,"Helicity Gated T3/I");
        t2->Branch(Form("hel.%s.T4",arm),&fT2,"Helicity Gated T4/I");
    }
    else if(arm[0]=='R'){
        t2->Branch(Form("hel.%s.T1",arm),&fT1,"Helicity Gated T1/I");
        t2->Branch(Form("hel.%s.T2",arm),&fT2,"Helicity Gated T2/I");
    }
    else if(arm[0]=='T'){
        t2->Branch(Form("hel.%s.T1",arm),&fT1,"Helicity Gated T1/I");
    }

    fscanf(fp2,"%d",&N);
    for(Int_t k=0;k<N;k++){
        if((arm[0]=='L')||(arm[0]=='R')){
            fscanf(fp2,"%d%d%d%d%x%d%d%d%d%d%d%d%d",&fEvNum,&fHelicity_act,&fHelicity_rep,&fQRT,&fSeed,&fError,&fTime,&fBCMup,&fBCMdown,&fL1A,&fT1,&fT2,&temp[0]);
        }
        else if(arm[0]=='T'){
            fscanf(fp2,"%d%d%d%d%x%d%d%d%d%d%d%d",&fEvNum,&fHelicity_act,&fHelicity_rep,&fQRT,&fSeed,&fError,&fTime,&fBCMup,&fBCMdown,&fL1A,&fT1,&temp[0]);
        }
        if((fEvNum>=gEvNumMin)&&(fEvNum<=gEvNumMax))t2->Fill();
    }

    t2->Write();

    if((arm[0]=='L')||(arm[0]=='R')){
        t3->Branch(Form("hel.%s.evnum",arm),&fEvNum,"Event Number/I");
        t3->Branch(Form("hel.%s.hel_rep",arm),&fHelicity_rep,"HAPPEX Reported Helicity/I");
        t3->Branch(Form("hel.%s.hel_act",arm),&fHelicity_act,"HAPPEX Actual Helicity/I");
        t3->Branch(Form("hel.%s.qrt",arm),&fQRT,"HAPPEX QRT/I");
        t3->Branch(Form("hel.%s.seed",arm),&fSeed,"HAPPEX Seed/I");
        t3->Branch(Form("hel.%s.error",arm),&fError,"HAPPEX Decode Error/I");
        t3->Branch(Form("hel.%s.bcmup",arm),&fBCMup,"Helicity Gated Upstream BCM/I");
        t3->Branch(Form("hel.%s.bcmdown",arm),&fBCMdown,"Helicity Gated Downstream BCM/I");

        fscanf(fp3,"%d",&N);
        for(Int_t k=0;k<N;k++){
            fscanf(fp3,"%d%d%d%d%x%d%d%d",&fEvNum,&fHelicity_act,&fHelicity_rep,&fQRT,&fSeed,&fError,&fBCMup,&fBCMdown);
            if((fEvNum>=gEvNumMin)&&(fEvNum<=gEvNumMax))t3->Fill();
        }

        t3->Write();
    }
    
    f->Close();

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    
    return 0;
}

Int_t insert(Int_t nrun,Char_t *filename)
{
    printf("Opening existed rootfile ...\n");

    fp1=fopen(Form("hel_%d.dat",nrun),"r");
    fp2=fopen(Form("helRIN_%d.dat",nrun),"r");
    fp3=fopen(Form("helHAP_%d.dat",nrun),"r");
    
    TFile *f=new TFile(filename,"UPDATE");

    TTree *t=(TTree *)f->Get("T");

    TTree *t1=new TTree("hel","TIR Helicity Tree");
    TTree *t2=new TTree("hel_ring","Ring Helicity Tree");
    TTree *t3=new TTree("hel_happ","HAPPEX Helicity Tree");

    Int_t fHelicity_rep=0,fHelicity_act=0,fQRT=0,fMPS=0,fPairSync=0,fTimeStamp=0,fError=0;
    Int_t fSeed=0;
    Int_t fBCMup=0,fBCMdown=0,fTime=0;
    Int_t fBCMuph=0,fBCMdownh=0;
    Int_t fIRing=0,fIHappex=0;
    Int_t fEvNum=0;
    Int_t fT1=0,fT2=0;
    Int_t fL1A=0;
    Int_t N=0;
    Int_t temp[10];
    Int_t gEvNum=0,gEvNumMin=100000000,gEvNumMax=0;
    Int_t nentries;
    THaEvent *event=new THaEvent();

    t->SetBranchAddress("Event_Branch",&event);
    
    t1->Branch(Form("hel.%s.hel_rep",arm),&fHelicity_rep,"TIR Reported Helicity/I");
    t1->Branch(Form("hel.%s.hel_act",arm),&fHelicity_act,"TIR Actual Helicity/I");
    t1->Branch(Form("hel.%s.qrt",arm),&fQRT,"TIR QRT/I");
    t1->Branch(Form("hel.%s.mps",arm),&fMPS,"TIR MPS/I");
    t1->Branch(Form("hel.%s.pairsync",arm),&fPairSync,"TIR PairSync/I");
    t1->Branch(Form("hel.%s.timestamp",arm),&fTimeStamp,"TIR TimeStamp/I");
    t1->Branch(Form("hel.%s.seed",arm),&fSeed,"TIR Seed/I");   
    t1->Branch(Form("hel.%s.error",arm),&fError,"TIR Decode Error/I");
    t1->Branch(Form("hel.%s.numring",arm),&fIRing,"Event num in the ring/I");   
    t1->Branch(Form("hel.%s.numhappex",arm),&fIHappex,"Event num in HAPPEX/I");  
    t1->Branch(Form("hel.%s.bcmup",arm),&fBCMup,"Helicity Gated Upstream BCM/I");
    t1->Branch(Form("hel.%s.bcmdown",arm),&fBCMdown,"Helicity Gated Downstream BCM/I");

    if((arm[0]=='L')||(arm[0]=='R')){
        t1->Branch(Form("hel.%s.bcmuphappex",arm),&fBCMuph,"Helicity Gated Upstream BCM/I");
        t1->Branch(Form("hel.%s.bcmdownhappex",arm),&fBCMdownh,"Helicity Gated Downstream BCM/I");
    }
    t1->Branch(Form("hel.%s.time",arm),&fTime,"Helicity Gated Time/I");

    fscanf(fp1,"%d",&N);
    nentries=t->GetEntries();
    for(Int_t k=0;k<nentries;k++){
        t->GetEntry(k);
        gEvNum=Int_t(event->GetHeader()->GetEvtNum());
        printf("%d\n",gEvNum);
        if(gEvNum<gEvNumMin)gEvNumMin=gEvNum;
        if(gEvNum>gEvNumMax)gEvNumMax=gEvNum;
        do{
            if((arm[0]=='L')||(arm[0]=='R')){
                fscanf(fp1,"%d%d%d%d%d%d%d%x%d%d%d%d%d%d%d%d",&temp[0],&fHelicity_act,&fHelicity_rep,&fQRT,&fPairSync,&fMPS,&fTimeStamp,&fSeed,&fError,&fIRing,&fIHappex,&fBCMup,&fBCMdown,&fTime,&fBCMuph,&fBCMdownh);
            }
            else{
                fscanf(fp1,"%d%d%d%d%d%d%d%x%d%d%d%d%d%d",&temp[0],&fHelicity_act,&fHelicity_rep,&fQRT,&fPairSync,&fMPS,&fTimeStamp,&fSeed,&fError,&fIRing,&fIHappex,&fBCMup,&fBCMdown,&fTime);
            }
        }while((temp[0]!=gEvNum)&&(!feof(fp1)));
        if(temp[0]==gEvNum){
            t1->Fill();
        }
        else{
            fHelicity_rep=0;
            fHelicity_act=0;
            fQRT=0;
            fPairSync=0;
            fMPS=0;
            fTimeStamp=0;
            fSeed=0;
            fError=0;
            fIRing=0;
            fIHappex=0;
            fBCMup=0;
            fBCMdown=0;
            fTime=0;
            fBCMuph=0;
            fBCMdownh=0;
            t1->Fill();
        }
    }

    t1->Write();
    
    t2->Branch(Form("hel.%s.evnum",arm),&fEvNum,"Event Number/I");
    t2->Branch(Form("hel.%s.hel_rep",arm),&fHelicity_rep,"Ring Reported Helicity/I");
    t2->Branch(Form("hel.%s.hel_act",arm),&fHelicity_act,"Ring Actual Helicity/I");
    t2->Branch(Form("hel.%s.qrt",arm),&fQRT,"Ring QRT/I");
    t2->Branch(Form("hel.%s.time",arm),&fTime,"Ring Time/I");
    t2->Branch(Form("hel.%s.seed",arm),&fSeed,"Ring Seed/I");   
    t2->Branch(Form("hel.%s.error",arm),&fError,"Ring Decode Error/I");
    t2->Branch(Form("hel.%s.bcmup",arm),&fBCMup,"Helicity Gated Upstream BCM/I");
    t2->Branch(Form("hel.%s.bcmdown",arm),&fBCMdown,"Helicity Gated Downstream BCM/I");
    t2->Branch(Form("hel.%s.L1A",arm),&fL1A,"Helicity Gated L1A/I");

    if(arm[0]=='L'){
        t2->Branch(Form("hel.%s.T3",arm),&fT1,"Helicity Gated T3/I");
        t2->Branch(Form("hel.%s.T4",arm),&fT2,"Helicity Gated T4/I");
    }
    else if(arm[0]=='R'){
        t2->Branch(Form("hel.%s.T1",arm),&fT1,"Helicity Gated T1/I");
        t2->Branch(Form("hel.%s.T2",arm),&fT2,"Helicity Gated T2/I");
    }
    else if(arm[0]=='T'){
        t2->Branch(Form("hel.%s.T1",arm),&fT1,"Helicity Gated T1/I");
    }

    fscanf(fp2,"%d",&N);
    for(Int_t k=0;k<N;k++){
        if((arm[0]=='L')||(arm[0]=='R')){
            fscanf(fp2,"%d%d%d%d%x%d%d%d%d%d%d%d%d",&fEvNum,&fHelicity_act,&fHelicity_rep,&fQRT,&fSeed,&fError,&fTime,&fBCMup,&fBCMdown,&fL1A,&fT1,&fT2,&temp[0]);
        }
        else if(arm[0]=='T'){
            fscanf(fp2,"%d%d%d%d%x%d%d%d%d%d%d%d",&fEvNum,&fHelicity_act,&fHelicity_rep,&fQRT,&fSeed,&fError,&fTime,&fBCMup,&fBCMdown,&fL1A,&fT1,&temp[0]);
        }
        if((fEvNum>=gEvNumMin)&&(fEvNum<=gEvNumMax))t2->Fill();
    }

    t2->Write();

    if((arm[0]=='L')||(arm[0]=='R')){
        t3->Branch(Form("hel.%s.evnum",arm),&fEvNum,"Event Number/I");
        t3->Branch(Form("hel.%s.hel_rep",arm),&fHelicity_rep,"HAPPEX Reported Helicity/I");
        t3->Branch(Form("hel.%s.hel_act",arm),&fHelicity_act,"HAPPEX Actual Helicity/I");
        t3->Branch(Form("hel.%s.qrt",arm),&fQRT,"HAPPEX QRT/I");
        t3->Branch(Form("hel.%s.seed",arm),&fSeed,"HAPPEX Seed/I");
        t3->Branch(Form("hel.%s.error",arm),&fError,"HAPPEX Decode Error/I");
        t3->Branch(Form("hel.%s.bcmup",arm),&fBCMup,"Helicity Gated Upstream BCM/I");
        t3->Branch(Form("hel.%s.bcmdown",arm),&fBCMdown,"Helicity Gated Downstream BCM/I");

        fscanf(fp3,"%d",&N);
        for(Int_t k=0;k<N;k++){
            fscanf(fp3,"%d%d%d%d%x%d%d%d",&fEvNum,&fHelicity_act,&fHelicity_rep,&fQRT,&fSeed,&fError,&fBCMup,&fBCMdown);
            if((fEvNum>=gEvNumMin)&&(fEvNum<=gEvNumMax))t3->Fill();
        }

        t3->Write();
    }
    
    f->Close();

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    
    return 0;
}

Int_t construct(Int_t nrun)
{
    printf("Generating a new rootfile ...\n");

    fp1=fopen(Form("hel_%d.dat",nrun),"r");
    fp2=fopen(Form("helRIN_%d.dat",nrun),"r");
    fp3=fopen(Form("helHAP_%d.dat",nrun),"r");

    TFile *f=new TFile(Form("hel_%d.root",nrun),"RECREATE");

    TTree *t1=new TTree("hel","TIR Helicity Tree");
    TTree *t2=new TTree("hel_ring","Ring Helicity Tree");
    TTree *t3=new TTree("hel_happ","HAPPEX Helicity Tree");

    Int_t fHelicity_rep=0,fHelicity_act=0,fQRT=0,fMPS=0,fPairSync=0,fTimeStamp=0,fError=0;
    Int_t fSeed=0;
    Int_t fBCMup=0,fBCMdown=0,fTime=0;
    Int_t fBCMuph=0,fBCMdownh=0;
    Int_t fIRing=0,fIHappex=0;
    Int_t fEvNum=0;
    Int_t fT1=0,fT2=0;
    Int_t fL1A=0;
    Int_t N=0;
    Int_t temp[10];

    t1->Branch(Form("hel.%s.hel_rep",arm),&fHelicity_rep,"TIR Reported Helicity/I");
    t1->Branch(Form("hel.%s.hel_act",arm),&fHelicity_act,"TIR Actual Helicity/I");
    t1->Branch(Form("hel.%s.qrt",arm),&fQRT,"TIR QRT/I");
    t1->Branch(Form("hel.%s.mps",arm),&fMPS,"TIR MPS/I");
    t1->Branch(Form("hel.%s.pairsync",arm),&fPairSync,"TIR PairSync/I");
    t1->Branch(Form("hel.%s.timestamp",arm),&fTimeStamp,"TIR TimeStamp/I");
    t1->Branch(Form("hel.%s.seed",arm),&fSeed,"TIR Seed/I");   
    t1->Branch(Form("hel.%s.error",arm),&fError,"TIR Decode Error/I");
    t1->Branch(Form("hel.%s.numring",arm),&fIRing,"Event num in the ring/I");   
    t1->Branch(Form("hel.%s.numhappex",arm),&fIHappex,"Event num in HAPPEX/I");  
    t1->Branch(Form("hel.%s.bcmupring",arm),&fBCMup,"Helicity Gated Upstream BCM/I");
    t1->Branch(Form("hel.%s.bcmdownring",arm),&fBCMdown,"Helicity Gated Downstream BCM/I");
    if((arm[0]=='L')||(arm[0]=='R')){
        t1->Branch(Form("hel.%s.bcmuphappex",arm),&fBCMuph,"Helicity Gated Upstream BCM/I");
        t1->Branch(Form("hel.%s.bcmdownhappex",arm),&fBCMdownh,"Helicity Gated Downstream BCM/I");
    }
    t1->Branch(Form("hel.%s.time",arm),&fTime,"Helicity Gated Time/I");

    fscanf(fp1,"%d",&N);
    for(Int_t k=0;k<N;k++){
        if((arm[0]=='L')||(arm[0]=='R')){
            fscanf(fp1,"%d%d%d%d%d%d%d%x%d%d%d%d%d%d%d%d",&temp[0],&fHelicity_act,&fHelicity_rep,&fQRT,&fPairSync,&fMPS,&fTimeStamp,&fSeed,&fError,&fIRing,&fIHappex,&fBCMup,&fBCMdown,&fTime,&fBCMuph,&fBCMdownh);           
        }
        else{
            fscanf(fp1,"%d%d%d%d%d%d%d%x%d%d%d%d%d%d",&temp[0],&fHelicity_act,&fHelicity_rep,&fQRT,&fPairSync,&fMPS,&fTimeStamp,&fSeed,&fError,&fIRing,&fIHappex,&fBCMup,&fBCMdown,&fTime);
        }
        t1->Fill();
    }

    t1->Write();
    
    t2->Branch(Form("hel.%s.evnum",arm),&fEvNum,"Event Number/I");
    t2->Branch(Form("hel.%s.hel_rep",arm),&fHelicity_rep,"Ring Reported Helicity/I");
    t2->Branch(Form("hel.%s.hel_act",arm),&fHelicity_act,"Ring Actual Helicity/I");
    t2->Branch(Form("hel.%s.qrt",arm),&fQRT,"Ring QRT/I");
    t2->Branch(Form("hel.%s.time",arm),&fTime,"Ring Time/I");
    t2->Branch(Form("hel.%s.seed",arm),&fSeed,"Ring Seed/I");   
    t2->Branch(Form("hel.%s.error",arm),&fError,"Ring Decode Error/I");
    t2->Branch(Form("hel.%s.bcmup",arm),&fBCMup,"Helicity Gated Upstream BCM/I");
    t2->Branch(Form("hel.%s.bcmdown",arm),&fBCMdown,"Helicity Gated Downstream BCM/I");
    t2->Branch(Form("hel.%s.L1A",arm),&fL1A,"Helicity Gated L1A/I");

    if(arm[0]=='L'){
        t2->Branch(Form("hel.%s.T3",arm),&fT1,"Helicity Gated T3/I");
        t2->Branch(Form("hel.%s.T4",arm),&fT2,"Helicity Gated T4/I");
    }
    else if(arm[0]=='R'){
        t2->Branch(Form("hel.%s.T1",arm),&fT1,"Helicity Gated T1/I");
        t2->Branch(Form("hel.%s.T2",arm),&fT2,"Helicity Gated T2/I");
    }
    else if(arm[0]=='T'){
        t2->Branch(Form("hel.%s.T1",arm),&fT1,"Helicity Gated T1/I");
    }

    fscanf(fp2,"%d",&N);
    for(Int_t k=0;k<N;k++){
        if((arm[0]=='L')||(arm[0]=='R')){
            fscanf(fp2,"%d%d%d%d%x%d%d%d%d%d%d%d%d",&fEvNum,&fHelicity_act,&fHelicity_rep,&fQRT,&fSeed,&fError,&fTime,&fBCMup,&fBCMdown,&fL1A,&fT1,&fT2,&temp[0]);
        }
        else if(arm[0]=='T'){
            fscanf(fp2,"%d%d%d%d%x%d%d%d%d%d%d%d",&fEvNum,&fHelicity_act,&fHelicity_rep,&fQRT,&fSeed,&fError,&fTime,&fBCMup,&fBCMdown,&fL1A,&fT1,&temp[0]);
        }
        t2->Fill();
    }

    t2->Write();

    if((arm[0]=='L')||(arm[0]=='R')){
        t3->Branch(Form("hel.%s.evnum",arm),&fEvNum,"Event Number/I");
        t3->Branch(Form("hel.%s.hel_rep",arm),&fHelicity_rep,"HAPPEX Reported Helicity/I");
        t3->Branch(Form("hel.%s.hel_act",arm),&fHelicity_act,"HAPPEX Actual Helicity/I");
        t3->Branch(Form("hel.%s.qrt",arm),&fQRT,"HAPPEX QRT/I");
        t3->Branch(Form("hel.%s.seed",arm),&fSeed,"HAPPEX Seed/I");
        t3->Branch(Form("hel.%s.error",arm),&fError,"HAPPEX Decode Error/I");
        t3->Branch(Form("hel.%s.bcmup",arm),&fBCMup,"Helicity Gated Upstream BCM/I");
        t3->Branch(Form("hel.%s.bcmdown",arm),&fBCMdown,"Helicity Gated Downstream BCM/I");

        fscanf(fp3,"%d",&N);
        for(Int_t k=0;k<N;k++){
            fscanf(fp3,"%d%d%d%d%x%d%d%d",&fEvNum,&fHelicity_act,&fHelicity_rep,&fQRT,&fSeed,&fError,&fBCMup,&fBCMdown);
            t3->Fill();
        }

        t3->Write();
    }
    
    f->Close();

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    
    return 0;
}

void usage(int argc, char** argv)
{
    printf("usage: %s [options] RUN_NUMBER\n", argv[0]);
    printf("  -c, --cfgfile=config.cfg   Set configuration file name\n");
    printf("  -h, --help                 This small usage guide\n");
    printf("  -i, --infodir=.            Set helicity info directory\n");
    printf("  -r, --rootdir=.            Set rootfile directory\n");
}
