#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "THaEvent.h"

#include "hel.h"

FILE *fp1,*fp2,*fp3;

Int_t insert(Int_t nrun,Char_t* filename);
Int_t construct(Int_t nrun);
Bool_t isexist(Char_t* fname);

int main(int argc,char* argv[])
{
    Int_t nrun=atoi(argv[1]);

    if(nrun<20000){
        strcpy(arm,"L");
        nring=7;
        nhap=2;
    }    
    else if(nrun<40000){
        strcpy(arm,"R");
        nring=7;
        nhap=2;
    }
    else{
        strcpy(arm,"TA");
        nring=6;
    }

    Int_t select=atoi(argv[2]);
    Int_t filecount=1;
    Char_t filename[300];

    if(select==1){
        sprintf(filename,"g2p_%d.root",nrun);
        insert(nrun,filename);
        sprintf(filename,"g2p_%d_%d.root",nrun,filecount);
        while(isexist(filename)){
            insert(nrun,filename);
            filecount++;
            sprintf(filename,"g2p_%d_%d.root",nrun,filecount);
        }
    }
    else if(select==2)construct(nrun);
    
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

Bool_t isexist(Char_t* fname)
{
    FILE *temp;
    Bool_t isopen;

    if((temp=fopen(fname,"r"))==NULL){
        isopen=false;
    }
    else{
        isopen=true;
        fclose(temp);
    }

  return isopen;
}


