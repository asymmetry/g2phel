#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "TROOT.h"

#include "hel.h"

#define LEN 15000000

FILE *fp,*fp1,*fp2;

//Global variables
Int_t gHelicity_rep[LEN];
Int_t gHelicity_act[LEN];
Int_t gQRT[LEN];
Int_t gSeed_rep[LEN];
Int_t gError[LEN];
Int_t gN;

Int_t readin(Int_t nrun,Char_t *name);
Int_t predictring(Int_t nrun);
Int_t printout(Int_t nrun,Char_t *name);
Int_t RanBit30(Int_t &runseed);
Int_t BitRan30(Int_t &runseed);

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
        nhap=0;
    }

    if(strcmp(argv[2],"HAP")==0)nring=nhap;
    
    readin(nrun,argv[2]);
    predictring(nrun);
    printout(nrun,argv[2]);
    
    return 0;
}

Int_t readin(Int_t nrun,Char_t *name)
{
    printf("Reading %s helicity information ...\n",name);

    fp=fopen(Form("hel%s_%d.tmp",name,nrun),"r");

    Int_t temp1;

    gN=0;

    fscanf(fp,"%d",&temp1);
    while(!feof(fp)){
        fscanf(fp,"%d%d",&gHelicity_rep[gN],&gQRT[gN]);
        for(Int_t i=0;i<nring;i++)
            fscanf(fp,"%d",&temp1);
        gN++;
        fscanf(fp,"%d",&temp1);
    }

    fclose(fp);
    
    return 0;
}

Int_t predictring(Int_t nrun){
    printf("Predicting helicity information ...\n");

    Int_t fPhaseRing_rep=0;
    Int_t fPolarityRing_rep=0,fPolarityRing_act=0;
    Int_t fSeedRing_rep=0,fSeedRing_act=0;
    Int_t fNSeedRing=0;
    
    for(Int_t i=0;i<gN;i++){
        gError[i]=0;
        gSeed_rep[i]=0;
        gHelicity_act[i]=0;
        
        if(gQRT[i]==1){
            fPhaseRing_rep=0;
        }
        else{
            fPhaseRing_rep+=1;
        }
        if(fPhaseRing_rep>=4){
            fNSeedRing=0;
            gError[i]|=0x04;
        }
        
        if((fNSeedRing==MAXBIT)&&gQRT[i]==1){
            fPolarityRing_rep=RanBit30(fSeedRing_rep);
            fPolarityRing_act=RanBit30(fSeedRing_act);
            if(fPolarityRing_rep!=gHelicity_rep[i]){
                fNSeedRing=0;
                gError[i]|=0x02;
            }
        }

        if((fNSeedRing<MAXBIT)&&gQRT[i]==1){
            fSeedRing_rep=((fSeedRing_rep<<1)&0x3FFFFFFF)|gHelicity_rep[i];
            fNSeedRing+=1;
            if(fNSeedRing==MAXBIT){
                fSeedRing_act=fSeedRing_rep;
                for(Int_t j=0;j<NDELAY;j++)fPolarityRing_act=RanBit30(fSeedRing_act);
                gError[i]=0;
            }
        }

        if(fNSeedRing==MAXBIT){
            if(fPolarityRing_act==1){
                if(fPhaseRing_rep==0||fPhaseRing_rep==3)
                    gHelicity_act[i]=1;
                else
                    gHelicity_act[i]=-1;
            }
            else{
                if(fPhaseRing_rep==0||fPhaseRing_rep==3)
                    gHelicity_act[i]=-1;
                else
                    gHelicity_act[i]=1;
            }
            gSeed_rep[i]=fSeedRing_rep;
        }
        else
        {
            gError[i]|=0x01;
            gSeed_rep[i]=0;
            gHelicity_act[i]=0;
        }
    }

    for(Int_t i=gN-2;i>=0;i--){
        if(gQRT[i+1]==1){
            fPhaseRing_rep=3;
        }
        else{
            fPhaseRing_rep-=1;
        }
        if(fPhaseRing_rep>=0){
            if(gError[i]==0){
                if(gQRT[i]==1){
                    fSeedRing_rep=gSeed_rep[i];
                }
            }
            else{
                if(gQRT[i]==1){
                    if(fSeedRing_rep!=0){
                        fPolarityRing_rep=BitRan30(fSeedRing_rep);
                        if(gHelicity_rep[i]==fPolarityRing_rep){
                            fSeedRing_act=fSeedRing_rep;
                            for(Int_t j=0;j<NDELAY;j++)fPolarityRing_act=RanBit30(fSeedRing_act);
                            if(fPolarityRing_act==1){
                                gHelicity_act[i]=gHelicity_act[i+3]=1;
                                gHelicity_act[i+1]=gHelicity_act[i+2]=-1;
                                for(Int_t j=0;j<4;j++){
                                    gError[i+j]=0;
                                    gSeed_rep[i+j]=fSeedRing_rep;
                                }
                            }
                            else{
                                gHelicity_act[i]=gHelicity_act[i+3]=-1;
                                gHelicity_act[i+1]=gHelicity_act[i+2]=1;
                                for(Int_t j=0;j<4;j++){
                                    gError[i+j]=0;
                                    gSeed_rep[i+j]=fSeedRing_rep;
                                }
                            }
                        }
                        else{
                            fSeedRing_rep=0;
                        }                        
                    }
                }
            }
        }
        else{
            fSeedRing_rep=0;
        }
    }

    return 0;
}

Int_t printout(Int_t nrun, Char_t* name)
{
    fp1=fopen(Form("hel%s_%d.tmp",name,nrun),"r");
    fp2=fopen(Form("hel%s_%d.dat",name,nrun),"w");

    Int_t temp[10];

    fprintf(fp2,"%d\n",gN);
    for(Int_t i=0;i<gN;i++){
        fscanf(fp1,"%d%d%d",&temp[0],&temp[1],&temp[2]);
        fprintf(fp2,"%d\t%d\t%d\t%d\t%08x\t%d",temp[0],gHelicity_act[i],gHelicity_rep[i],gQRT[i],gSeed_rep[i],gError[i]);
        for(Int_t k=0;k<nring;k++){
            fscanf(fp1,"%d",&temp[3]);
            fprintf(fp2,"\t%d",temp[3]);
        }
        fprintf(fp2,"\n");
    }

    fclose(fp1);
    fclose(fp2);

    return 0;
}

Int_t RanBit30(Int_t& ranseed)
{
  // Take 7,28,29,30 bit of ranseed out
  UInt_t bit7    = ((ranseed&0x00000040)!=0);
  UInt_t bit28   = ((ranseed&0x08000000)!=0);
  UInt_t bit29   = ((ranseed&0x10000000)!=0);
  UInt_t bit30   = ((ranseed&0x20000000)!=0);

  UInt_t newbit  = (bit30^bit29^bit28^bit7)&0x1;

  if(ranseed<=0){
    newbit=0;
  }
  
  ranseed=((ranseed<<1)|newbit)&0x3FFFFFFF;

  return newbit;
}

Int_t BitRan30(Int_t& ranseed)
{
  // Take 7,28,29,30 bit of ranseed out
  UInt_t bit1    = ((ranseed&0x00000001)!=0);
  UInt_t bit8    = ((ranseed&0x00000080)!=0);
  UInt_t bit29   = ((ranseed&0x10000000)!=0);
  UInt_t bit30   = ((ranseed&0x20000000)!=0);

  UInt_t newbit  = (bit30^bit29^bit8^bit1)&0x1;

  if(ranseed<=0){
    newbit=0;
  }
  
  ranseed=((ranseed>>1)|(newbit<<29))&0x3FFFFFFF;

  newbit=((ranseed&0x1)!=0);

  return newbit;
}
