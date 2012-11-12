#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "hel.h"

#define LEN 15000000

FILE *fp,*fp1,*fp2;

//Global variables
Int_t gHelicity_act[LEN];
Int_t gSeed_rep[LEN];
Int_t gPhase[LEN];
Int_t gBCMu[LEN];
Int_t gBCMd[LEN];
Int_t gTime[LEN];
Int_t gError[LEN];
Int_t gN;

Int_t readin(Int_t nrun);
Int_t readinh(Int_t nrun);
Int_t align(Int_t nrun);
Int_t alignh(Int_t nrun);

int main(int argc,char* argv[])
{
    Int_t nrun=atoi(argv[1]);

    if(nrun<20000){
        strcpy(arm,"L");
        nring=7;
    }    
    else if(nrun<40000){
        strcpy(arm,"R");
        nring=7;
    }
    else{
        strcpy(arm,"TA");
        nring=6;
    }
    
    readin(nrun);
    align(nrun);

    if(nrun<40000){
        readinh(nrun);
        alignh(nrun);
    }
    
    return 0;
}

Int_t readin(Int_t nrun)
{
    fp=fopen(Form("helRIN_%d.dat",nrun),"r");

    Int_t fHelRing_rep,fQRTRing,fPhase=0,temp1;
    
    fscanf(fp,"%d",&gN);
    for(Int_t i=0;i<gN;i++){
        fscanf(fp,"%d%d%d%d%x%d%d%d%d",&temp1,&gHelicity_act[i],&fHelRing_rep,&fQRTRing,&gSeed_rep[i],&gError[i],&gTime[i],&gBCMu[i],&gBCMd[i]);
        if(fQRTRing==1)fPhase=0;
        else fPhase++;
        if(gError[i]>0)fPhase=0;
        gPhase[i]=fPhase;
        for(Int_t k=0;k<nring-3;k++)fscanf(fp,"%d",&temp1);
    }

    fclose(fp);
    
    return 0;
}

Int_t readinh(Int_t nrun)
{
    fp=fopen(Form("helHAP_%d.dat",nrun),"r");

    Int_t fHelRing_rep,fQRTRing,fPhase=0,temp1;

    for(Int_t i=0;i<LEN;i++){
        gHelicity_act[i]=0;
        gSeed_rep[i]=0;
        gError[i]=0;
        gTime[i]=0;
        gBCMu[i]=0;
        gBCMd[i]=0;
        gPhase[i]=0;
    }
    
    fscanf(fp,"%d",&gN);
    for(Int_t i=0;i<gN;i++){
        fscanf(fp,"%d%d%d%d%x%d%d%d",&temp1,&gHelicity_act[i],&fHelRing_rep,&fQRTRing,&gSeed_rep[i],&gError[i],&gBCMu[i],&gBCMd[i]);
        if(fQRTRing==1)fPhase=0;
        else fPhase++;
        if(gError[i]>0)fPhase=0;
        gPhase[i]=fPhase;
    }

    fclose(fp);
    
    return 0;
}

Int_t align(Int_t nrun)
{
    printf("Constructing new TIR helicity information ...\n");

    fp1=fopen(Form("helTIR_%d.dat",nrun),"r");
    if(nrun<40000){
        fp2=fopen(Form("hel_%d.tmp",nrun),"w");
    }
    else{
        fp2=fopen(Form("hel_%d.dat",nrun),"w");
    }
    
    Int_t fHelicity_rep=0,fHelicity_act=0,fQRT=0,fMPS=0,fPairSync=0,fTimeStamp=0,fError=0;
    Int_t temp[10];
    Int_t fSeedTIR_rep=0,fPhaseTIR=0,fPolarityTIR_rep=0;
    Int_t pLastp=0,pLastm=0;
    Int_t fBCMRingu,fBCMRingd,fTimeStampRing;
    Bool_t fNewFlag=kTRUE;
    Int_t N=0;

    fscanf(fp1,"%d",&N);
    fprintf(fp2,"%d\n",N);
    for(Int_t k=0;k<N;k++){
        fscanf(fp1,"%d%d%d%d%d%d%d%x%d%d%d",&temp[0],&fHelicity_act,&fHelicity_rep,&fQRT,&fPairSync,&fMPS,&fTimeStamp,&fSeedTIR_rep,&fError,&temp[1],&temp[2]); 
        fBCMRingu=0;
        fBCMRingd=0;
        fTimeStampRing=0;
        if(fError==0){
            fPolarityTIR_rep=fSeedTIR_rep&0x01;
            if(fQRT==1){
                fPhaseTIR=0;
            }
            else if(fHelicity_rep==fPolarityTIR_rep){
                fPhaseTIR=3;
            }
            else if(fPairSync==1){
                fPhaseTIR=2;
            }
            else{
                fPhaseTIR=1;
            }
            if(fNewFlag){
                for(Int_t i=(pLastp<pLastm)?pLastp:pLastm+1;i<gN;i++){
                    if((gPhase[i]==fPhaseTIR)&&(gSeed_rep[i]==fSeedTIR_rep)){
                        pLastp=i-1;
                        pLastm=i-1;
                        fNewFlag=kFALSE;
                        break;
                    }
                }
                if(!((gPhase[pLastp+1]==fPhaseTIR)&&(gSeed_rep[pLastp+1]==fSeedTIR_rep))){
                    fError=fError|0x20;
                    fNewFlag=kTRUE;
                    fprintf(fp2,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%08x\t%d\t%d\t%d\t%d\t%d\t%d\n",temp[0],fHelicity_act,fHelicity_rep,fQRT,fPairSync,fMPS,fTimeStamp,fSeedTIR_rep,fError,temp[1],temp[2],fBCMRingu,fBCMRingd,fTimeStampRing);
                    continue;
                }
            }
            if(fHelicity_act==1){
                if(!((gPhase[pLastp]==fPhaseTIR)&&(gSeed_rep[pLastp]==fSeedTIR_rep))){
                    for(Int_t i=pLastp+1;(i<pLastp+1000)&&(i<gN);i++){
                        if(gHelicity_act[i]==1){
                            fBCMRingu+=gBCMu[i];
                            fBCMRingd+=gBCMd[i];
                            fTimeStampRing+=gTime[i];
                        }
                        if((gPhase[i]==fPhaseTIR)&&(gSeed_rep[i]==fSeedTIR_rep)){
                            pLastp=i;
                            break;
                        }
                    }
                    if((gPhase[pLastp]==fPhaseTIR)&&(gSeed_rep[pLastp]==fSeedTIR_rep)){
                        fError=fError;
                    }
                    else{
                        fError=fError|0x10;
                        fNewFlag=kTRUE;
                    }
                }
            }
            else{
                if(!((gPhase[pLastm]==fPhaseTIR)&&(gSeed_rep[pLastm]==fSeedTIR_rep))){
                    for(Int_t i=pLastm+1;(i<pLastm+1000)&&(i<gN);i++){
                        if(gHelicity_act[i]==-1){
                            fBCMRingu+=gBCMu[i];
                            fBCMRingd+=gBCMd[i];
                            fTimeStampRing+=gTime[i];
                        }
                        if((gPhase[i]==fPhaseTIR)&&(gSeed_rep[i]==fSeedTIR_rep)){
                            pLastm=i;
                            break;
                        }
                    }
                    if((gPhase[pLastm]==fPhaseTIR)&&(gSeed_rep[pLastm]==fSeedTIR_rep)){
                        fError=fError;
                    }
                    else{
                        fError=fError|0x10;
                        fNewFlag=kTRUE;
                    }
                }
            }
        }
        else if(fError==8){
            fError=fError;
        }
        else{
            fError=fError;
            fNewFlag=kTRUE;
        }
        fprintf(fp2,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%08x\t%d\t%d\t%d\t%d\t%d\t%d\n",temp[0],fHelicity_act,fHelicity_rep,fQRT,fPairSync,fMPS,fTimeStamp,fSeedTIR_rep,fError,temp[1],temp[2],fBCMRingu,fBCMRingd,fTimeStampRing); 
    }

    fclose(fp1);
    fclose(fp2);
    
    return 0;
}

Int_t alignh(Int_t nrun)
{
    fp1=fopen(Form("hel_%d.tmp",nrun),"r");
    fp2=fopen(Form("hel_%d.dat",nrun),"w");

    Int_t fHelicity_rep=0,fHelicity_act=0,fQRT=0,fMPS=0,fPairSync=0,fTimeStamp=0,fError=0;
    Int_t temp[10];
    Int_t fSeedTIR_rep=0,fPhaseTIR=0,fPolarityTIR_rep=0;
    Int_t pLastp=0,pLastm=0;
    Int_t fBCMRingu,fBCMRingd;
    Bool_t fNewFlag=kTRUE;
    Int_t N=0;

    fscanf(fp1,"%d",&N);
    fprintf(fp2,"%d\n",N);
    for(Int_t k=0;k<N;k++){
        fscanf(fp1,"%d%d%d%d%d%d%d%x%d%d%d%d%d%d",&temp[0],&fHelicity_act,&fHelicity_rep,&fQRT,&fPairSync,&fMPS,&fTimeStamp,&fSeedTIR_rep,&fError,&temp[1],&temp[2],&temp[3],&temp[4],&temp[5]); 
        fBCMRingu=0;
        fBCMRingd=0;
        if(fError==0){
            fPolarityTIR_rep=fSeedTIR_rep&0x01;
            if(fQRT==1){
                fPhaseTIR=0;
            }
            else if(fHelicity_rep==fPolarityTIR_rep){
                fPhaseTIR=3;
            }
            else if(fPairSync==1){
                fPhaseTIR=2;
            }
            else{
                fPhaseTIR=1;
            }
            if(fNewFlag){
                for(Int_t i=(pLastp<pLastm)?pLastp:pLastm+1;i<gN;i++){
                    if((gPhase[i]==fPhaseTIR)&&(gSeed_rep[i]==fSeedTIR_rep)){
                        pLastp=i-1;
                        pLastm=i-1;
                        fNewFlag=kFALSE;
                        break;
                    }
                }
                if(!((gPhase[pLastp+1]==fPhaseTIR)&&(gSeed_rep[pLastp+1]==fSeedTIR_rep))){
                    fError=fError|0x200;
                    fNewFlag=kTRUE;
                    fprintf(fp2,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%08x\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",temp[0],fHelicity_act,fHelicity_rep,fQRT,fPairSync,fMPS,fTimeStamp,fSeedTIR_rep,fError,temp[1],temp[2],temp[3],temp[4],temp[5],fBCMRingu,fBCMRingd);
                    continue;
                }
            }
            if(fHelicity_act==1){
                if(!((gPhase[pLastp]==fPhaseTIR)&&(gSeed_rep[pLastp]==fSeedTIR_rep))){
                    for(Int_t i=pLastp+1;(i<pLastp+1000)&&(i<gN);i++){
                        if(gHelicity_act[i]==1){
                            fBCMRingu+=gBCMu[i];
                            fBCMRingd+=gBCMd[i];
                        }
                        if((gPhase[i]==fPhaseTIR)&&(gSeed_rep[i]==fSeedTIR_rep)){
                            pLastp=i;
                            break;
                        }
                    }
                    if((gPhase[pLastp]==fPhaseTIR)&&(gSeed_rep[pLastp]==fSeedTIR_rep)){
                        fError=fError;
                    }
                    else{
                        fError=fError|0x100;
                        fNewFlag=kTRUE;
                    }
                }
            }
            else{
                if(!((gPhase[pLastm]==fPhaseTIR)&&(gSeed_rep[pLastm]==fSeedTIR_rep))){
                    for(Int_t i=pLastm+1;(i<pLastm+1000)&&(i<gN);i++){
                        if(gHelicity_act[i]==-1){
                            fBCMRingu+=gBCMu[i];
                            fBCMRingd+=gBCMd[i];
                        }
                        if((gPhase[i]==fPhaseTIR)&&(gSeed_rep[i]==fSeedTIR_rep)){
                            pLastm=i;
                            break;
                        }
                    }
                    if((gPhase[pLastm]==fPhaseTIR)&&(gSeed_rep[pLastm]==fSeedTIR_rep)){
                        fError=fError;
                    }
                    else{
                        fError=fError|0x100;
                        fNewFlag=kTRUE;
                    }
                }
            }
        }
        else if(fError==8){
            fError=fError;
        }
        else{
            fError=fError;
            fNewFlag=kTRUE;
        }
        fprintf(fp2,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%08x\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",temp[0],fHelicity_act,fHelicity_rep,fQRT,fPairSync,fMPS,fTimeStamp,fSeedTIR_rep,fError,temp[1],temp[2],temp[3],temp[4],temp[5],fBCMRingu,fBCMRingd); 
    }

    fclose(fp1);
    fclose(fp2);
    
    return 0;
}


