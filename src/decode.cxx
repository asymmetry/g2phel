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

Int_t extract(Int_t nrun);
Int_t decode_hel(Int_t* data);
Int_t findword(Int_t* data, struct rocinfo info);
Int_t adc18_decode_data(Int_t data, Int_t adcnum, Int_t &num, Int_t &val);
void usage(int argc, char** argv);

#include "isexist.h"

Int_t EVTLIMIT = -1;
Char_t CFGFILE[300] = "./config.cfg";
Char_t RAWPATH[300] = ".";
Char_t OUTPATH[300] = ".";

int main(int argc, char** argv)
{
    int c;

    while (1) {
        static struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"event", required_argument, 0, 'e'},
            {"cfgfile", required_argument, 0, 'c'},
            {"rawpath", required_argument, 0, 'r'},
            {"outpath", required_argument, 0, 'o'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "c:e:ho:r:", long_options, &option_index);

        if (c==-1) break;
        
        switch (c) {
        case 'c':
            strcpy(CFGFILE, optarg);
            break;
        case 'e':
            EVTLIMIT = atoi(optarg);
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'o':
            strcpy(OUTPATH, optarg);
            break;
        case 'r':
            strcpy(RAWPATH, optarg);
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

    config_init(&cfg);

    if (!config_read_file(&cfg, CFGFILE)) {
        fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg), config_error_line(&cfg), config_error_text(&cfg));
        config_destroy(&cfg);
        exit(EXIT_FAILURE);
    }

    Bool_t configerror = kFALSE;
    
    setting = config_lookup(&cfg, "rocinfo.hel");
    if (setting != NULL) {
        config_setting_lookup_int(setting, "roc", &info_HEL.roc);
        int temp;
        config_setting_lookup_int(setting, "header", &temp);
        info_HEL.header=(unsigned)temp;
        config_setting_lookup_int(setting, "index", &info_HEL.index);
    }
    else configerror = kTRUE;
    setting = config_lookup(&cfg, "rocinfo.ring");
    if (setting != NULL) {
        config_setting_lookup_int(setting, "roc", &info_RIN.roc);
        int temp;
        config_setting_lookup_int(setting, "header", &temp);
        info_RIN.header=(unsigned)temp;
        config_setting_lookup_int(setting, "index", &info_RIN.index);
    }
    else configerror = kTRUE;
    setting = config_lookup(&cfg, "rocinfo.time");
    if (setting != NULL) {
        config_setting_lookup_int(setting, "roc", &info_TIM.roc);
        int temp;
        config_setting_lookup_int(setting, "header", &temp);
        info_TIM.header=(unsigned)temp;
        config_setting_lookup_int(setting, "index", &info_TIM.index);
    }
    else configerror = kTRUE;
    setting = config_lookup(&cfg, "rocinfo.happex");
    if (setting != NULL) {
        USEHAPPEX=kTRUE;
        config_setting_lookup_int(setting, "roc", &info_HAP.roc);
        int temp;
        config_setting_lookup_int(setting, "header", &temp);
        info_HAP.header=(unsigned)temp;
        config_setting_lookup_int(setting, "index", &info_HAP.index);
    }
    else {
        USEHAPPEX=kFALSE;
    }

    setting = config_lookup(&cfg, "ringinfo.data");
    if (setting != NULL) {
        NRING=config_setting_length(setting);
    }
    else configerror = kTRUE;
    if (USEHAPPEX) {
        setting = config_lookup(&cfg, "happexinfo.data");
        if (setting != NULL) {
        NHAPPEX=config_setting_length(setting);
        }
        else configerror = kTRUE;
    }

    if (configerror) {
        printf("Invalid cfg file !\n");
        exit(-1);
    }
    
    extract(nrun);

    return 0;
}

Int_t extract(Int_t nrun)
{
    printf("Extracting helicity information from run %d ...\n",nrun);

    Char_t filename[300];
    
    Int_t filecount=0;
    sprintf(filename, "%s/g2p_%d.dat.%d", RAWPATH, nrun, filecount);

    THaCodaData *coda;
    
    coda = new THaCodaFile();

    fp1=fopen(Form("%s/helTIR_%d.tmp", OUTPATH, nrun),"w");
    fp2=fopen(Form("%s/helRIN_%d.tmp", OUTPATH, nrun),"w");
    if (USEHAPPEX) fp3=fopen(Form("%s/helHAP_%d.tmp", OUTPATH, nrun),"w");

    Int_t status,*data;
    Int_t evcount=0;

    while (isexist(filename)) {
        if(coda->codaOpen(filename)==0){
            printf("Adding %s ...\n",filename);
        }
        else{
            break;
        }
        status=coda->codaRead();
        while(status==0){
            if((EVTLIMIT!=-1)&&(evnum>=EVTLIMIT))break;
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
        sprintf(filename, "%s/g2p_%d.dat.%d", RAWPATH, nrun, filecount);
    }
    
    fclose(fp1);
    fclose(fp2);
    if (USEHAPPEX) fclose(fp3);
    
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
            for(Int_t k=1;k<NRING;k++){
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
    
    if (USEHAPPEX) {
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
                for(Int_t k=0;k<NHAPPEX;k++){
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
    printf("  -c, --cfgfile=config.cfg   Set configuration file name\n");
    printf("  -e, --event=-1             Set event limit\n");
    printf("  -h, --help                 This small usage guide\n");
    printf("  -o, --outpath=.            Set output path\n");
    printf("  -r, --rawpath=.            Set rawdata path\n");
}
