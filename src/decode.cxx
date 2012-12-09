#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <libconfig.h>

#include "TROOT.h"
//#include "TFile.h"
//#include "TTree.h"
//#include "TMath.h"
//#include "TChain.h"
#include "THaCodaData.h"
#include "THaCodaFile.h"

#include "hel.h"
#include "decode.h"

FILE *fp1, *fp2, *fp3;

Int_t extract(Int_t nrun, Int_t total);
Int_t decode_hel(Int_t* data);
Int_t findword(Int_t* data, struct rocinfo info);
Int_t adc18_decode_data(Int_t data, Int_t adcnum, Int_t &num, Int_t &val);
void usage(int argc, char** argv);

#include "isexist.h"

Bool_t USEBIN = kFALSE;
Int_t EVTLIMIT = -1;
Char_t CFGFILE[] = "./config.cfg";
Char_t RAWPATH[] = ".";
Char_t OUTPATH[] = ".";

int main(int argc, char** argv)
{
    int c;

    while (1) {
        static struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"bin", no_argument, 0, 'b'},
            {"event", required_argument, 0, 'e'},
            {"cfgfile", required_argument, 0, 'c'},
            {"rawpath", required_argument, 0, 'r'},
            {"outpath", required_argument, 0, 'o'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "bc:e:ho:r:", long_options, &option_index);

        if (c==-1) break;

        switch (c) {
        case 'b':
            USEBIN = kTRUE;
            break;
        case 'c':
            strcpy(CFGFILE,optarg);
            break;
        case 'e':
            EVTLIMIT = atoi(optarg);
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'o':
            strcpy(OUTPATH,optarg);
            break;
        case 'r':
            strcpy(RAWPATH,optarg);
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
    else{
        usage(argc, argv);
        exit(-1);
    }

    printf("%d\n%d\n%s\n%s\n%s\n%d\n",USEBIN,EVTLIMIT,CFGFILE,RAWPATH,OUTPATH,nrun);

    exit(0);

    if(nrun<20000){
        strcpy(arm,"L");
        nring=7;
        nhap=2;
        info_HEL.roc=11;
        info_RIN.roc=11;
        info_HAP.roc=25;
        info_TIM.roc=11;
        info_HEL.header=0;
        info_RIN.header=0xfb1b0000;
        info_HAP.header=0xbf1ff000;
        info_TIM.header=0;
        info_HEL.index=3;
        info_RIN.index=0;
        info_HAP.index=0;
        info_TIM.index=4;        
    }    
    else if(nrun<40000){
        strcpy(arm,"R");
        nring=7;
        nhap=2;
        info_HEL.roc=10;
        info_RIN.roc=10;
        info_HAP.roc=26;
        info_TIM.roc=10;
        info_HEL.header=0;
        info_RIN.header=0xfb1b0000;
        info_HAP.header=0;
        info_TIM.header=0;
        info_HEL.index=3;
        info_RIN.index=0;
        info_HAP.index=0;
        info_TIM.index=4;
    }
    else{
        strcpy(arm,"TA");
        nring=6;
        nhap=0;
        info_HEL.roc=12;
        info_RIN.roc=12;
        info_HAP.roc=0;
        info_TIM.roc=12;
        info_HEL.header=0;
        info_RIN.header=0xfb1b0000;
        info_HAP.header=0;
        info_TIM.header=0;
        info_HEL.index=3;
        info_RIN.index=0;
        info_HAP.index=0;
        info_TIM.index=4; 
    }

    extract(nrun, EVTLIMIT);

    return 0;
}

Int_t extract(Int_t nrun, Int_t total)
{
    printf("Extracting helicity information from run %d ...\n",nrun);

    Char_t filename[300];
    
    Int_t filecount=0;
    sprintf(filename,"g2p_%d.dat.%d",nrun,filecount);

    THaCodaData *coda;
    
    coda = new THaCodaFile();

    fp1=fopen(Form("helTIR_%d.tmp",nrun),"w");
    fp2=fopen(Form("helRIN_%d.tmp",nrun),"w");
    fp3=fopen(Form("helHAP_%d.tmp",nrun),"w");

    Int_t status,*data;
    Int_t evcount=0;

    while(isexist(filename)){
        if(coda->codaOpen(filename)==0){
            printf("Adding %s ...\n",filename);
        }
        else{
            break;
        }
        status=coda->codaRead();
        while(status==0){
            if((total!=-1)&&(evnum>=total))break;
            evcount++;
            if(evcount/1000*1000==evcount)printf("%d\n",evcount);
            data=coda->getEvBuffer();
            evlen=data[0]+1;
            evtype=data[1]>>16;
            switch(evtype){
            case 1:
            case 2:
            case 3:
            case 4:
            case 8:
                evnum=data[4];
                decode_hel(data);
                break;
            }
            status=coda->codaRead();
        }
        filecount++;
        sprintf(filename,"g2p_%d.dat.%d",nrun,filecount);
    }
    
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    
    return 0;
}

Int_t decode_hel(Int_t* data)
{
    Int_t pos,nroc;
    Int_t len,iroc;
    Int_t index;
    Int_t d,id=0,ida,idaw;
    Int_t fHelicityTIR=0,fQRTTIR=0,fPairSyncTIR=0,fMPSTIR=0;
    Int_t fTimeStampTIR,fIRing,fIHappex;
    Int_t fHelicityRing,fQRTRing;
    Int_t fDRing[8];
    Int_t fHelicityHappex[50],fQRTHappex[50];
    Int_t fDHappex[8][50];
    Int_t num,val;

    //get ROC table
    pos=data[2]+3;
    nroc=0;
    while((pos+1<evlen)&&(nroc<MAXROC)){
        len=data[pos];
        iroc=(data[pos+1]&0xff0000)>>16;
        rocpos[iroc]=pos;
        roclen[iroc]=len;
        irn[nroc++]=iroc;
        pos+=len+1;
    }

    //decode helicity
    //modify it with experiment setting
    if((index=findword(data,info_HEL))!=-1){
        d=data[index];
        if(info_HEL.roc==10){
            fHelicityTIR=(( d)&0x20)>>5;
            fQRTTIR     =(( d)&0x10)>>4;
            fPairSyncTIR=((~d)&0x40)>>6;
            fMPSTIR     =((~d)&0x80)>>7;
        }
        else if(info_HEL.roc==11){
            fHelicityTIR=(( d)&0x10)>>4;
            fQRTTIR     =(( d)&0x20)>>5;
            fPairSyncTIR=((~d)&0x40)>>6;
            fMPSTIR     =((~d)&0x80)>>7;
        }
        else if(info_HEL.roc==12){
            fHelicityTIR=((~d)&0x20)>>5;
            fQRTTIR     =((~d)&0x10)>>4;
            fPairSyncTIR=((~d)&0x40)>>6;
            fMPSTIR     =((~d)&0x80)>>7;
        }
        fprintf(fp1,"%8d\t%d\t%d\t%d\t%d\t",evnum,fHelicityTIR,fQRTTIR,fPairSyncTIR,fMPSTIR);
    }
    else{
        return -1;
    }

    if((index=findword(data,info_TIM))!=-1){
        fTimeStampTIR=data[index];
        fprintf(fp1,"%10d\t",fTimeStampTIR);
    }
    else{
        return -1;
    }

    if((index=findword(data,info_RIN))!=-1){
        fIRing=data[index++]&0x3ff;
        for(Int_t i=0;i<fIRing;i++){
            fDRing[0]=data[index++];
            d=data[index++];
            fHelicityRing=(d&0x01);
            fQRTRing     =(d&0x10)>>4;
            fprintf(fp2,"%8d\t%d\t%d\t%3d",evnum,fHelicityRing,fQRTRing,fDRing[0]);
            for(Int_t k=1;k<nring;k++){
                fDRing[k]=data[index++];
                fprintf(fp2,"\t%3d",fDRing[k]);
            }
            fprintf(fp2,"\n");
        }
    }
    else{
        fIRing=0;
    }
    fprintf(fp1,"%2d\t",fIRing);
    
    if(info_HAP.roc!=0){
        if((index=findword(data,info_HAP))!=-1){
            d=data[index++];
            if(((d&0xbf1ff000)==0xbf1ff000)&&((d&0xffffffff)!=0xffffffff)){
                fIHappex=d&0xfff;
            }
            else{
                fIHappex=0;
            }
            for(Int_t i=index;i<rocpos[info_HAP.roc]+roclen[info_HAP.roc];i++){
                d=data[i];
                if(((d&0xbfead000)==0xbfead000)&&((d&0xffffffff)!=0xffffffff)){
                    id=d&0xfff;
                }
                if(((d&0xbffec000)==0xbffec000)&&((d&0xffffffff)!=0xffffffff)){
                    d=data[i+1];
                    fHelicityHappex[id]=(d&2)>>1;
                    fQRTHappex[id]=(d&4)>>2;
                }
                if(((d&0xbfadc000)==0xbfadc000)&&((d&0xffffffff)!=0xffffffff)){
                    ida=d&0x0f;
                    idaw=(d&0xf0)>>4;
                    for(Int_t j=0;j<ida;j++){
                        for(Int_t k=i+2+j*idaw;k<i+2+(j+1)*idaw;k++){
                            d=data[k];
                            adc18_decode_data(d,j,num,val);
                            if(num==4)fDHappex[0][id]=val;
                            if(num==5)fDHappex[1][id]=val;
                        }
                    }
                }
            }
            for(Int_t i=0;i<fIHappex;i++){
                fprintf(fp3,"%8d\t%d\t%d",evnum,fHelicityHappex[i],fQRTHappex[i]);
                for(Int_t k=0;k<nhap;k++){
                    fprintf(fp3,"\t%6d",fDHappex[k][i]);
                }
                fprintf(fp3,"\n");
            }
        }
        else{
            fIHappex=0;
        }
    }
    else{
        fIHappex=0;
    }
    fprintf(fp1,"%2d\n",fIHappex);

    return 0;
}

Int_t findword(Int_t* data, struct rocinfo info)
{
    Int_t i;

    if(info.header==0)
        i=info.index+rocpos[info.roc];
    else{
        for(i=rocpos[info.roc];(i<rocpos[info.roc]+roclen[info.roc]-4)&&((data[i]&0xfffff000)!=info.header);i++);
        i+=info.index;
    }
    
    return (i<rocpos[info.roc]+roclen[info.roc]-4)?i:-1;
}

Int_t adc18_decode_data(Int_t data,Int_t adcnum,Int_t &num,Int_t &val)
{
    Int_t module_id,event_number;
    Int_t ch_number,divider,div_n,data_type;
    Int_t diff_value;
    Int_t sign,difference;
    Int_t ii;
    Double_t diff_avg;

    if(data&0x80000000){
        module_id=(0x1F)&(data>>26);
        event_number=data&0x3FFFFFF;
    }
    else{
        ch_number=(0x3)&(data>>29);
        div_n=((0x3)&(data>>25));
        divider=1;
        for(ii=0;ii<div_n;ii++)divider=divider*2;
        data_type=(0x7)&(data>>22);
        if(data_type==0){
            diff_value=(0x1FFFFF)&data;
            if(data&0x200000){
                sign=-1;
                difference=sign*((~diff_value&0x1FFFFF)+1);
            }
            else{
                sign=1;
                difference=diff_value;
            }
            diff_avg=((Float_t)difference)/((Float_t)divider);
            if(ch_number>=0&&ch_number<4){
                num=ch_number+adcnum*4;
                val=(Int_t)diff_avg;
            }
        }
    }
    return 0;
}

void usage(int argc, char** argv)
{
    printf("usage: %s [options] RUN_NUMBER\n", argv[0]);
    printf("  -b, --bin                  Set output binary file\n");
    printf("  -c, --cfgfile=config.cfg   Set configuration file name\n");
    printf("  -e, --event=-1             Set event limit\n");
    printf("  -h, --help                 This small usage guide\n");
    printf("  -o, --outpath=.            Set output path\n");
    printf("  -r, --rawpath=.            Set rawdata path\n");
}
