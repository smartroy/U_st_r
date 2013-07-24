#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "temppred.h"
#include "sorttask.h"
//#include "vardeclare.h"
extern double tmpr_array[4][6][6];
//extern double ratio[2];
extern double ratio_v;
extern double ratio_h;
extern FILE *fp[2];
extern double power_standby[40];
extern double pwr[2];
extern int end0;
extern int end1;
extern int heat_0,heat_1;
extern struct task *sorted_queue;
extern double blk_rel_ratio[17];

struct task *send_to_exe_queue(struct task *exe_queue){
  struct task *q_temp;
  
  if(exe_queue==NULL){
    q_temp=sorted_queue;
    sorted_queue=sorted_queue->next;
    exe_queue=q_temp;
    exe_queue->next=NULL;
  }else{
    q_temp=sorted_queue;//save next in sorted q since next function will change the struct pointed by sorted_queue
    sorted_queue=sorted_queue->next;
    exe_queue=sort_deadline(q_temp,exe_queue);
    
  }
  return exe_queue;
}

float exp_predict(float T_ss, float time_past, float T_init){

   //float tss;
   //float t_current;
   //float time_past;
   float t_predict;
   //float t_init;
   
   //tss=(t_current-t_init*exp(-0.015*time_past))/(1-exp(-0.015*time_past));  //0.015 is the chosen parameter  a in exponential function
   t_predict=T_ss-(T_ss-T_init)*exp(-0.017*(time_past+0.3));
   
   return t_predict;
}

void tmpr_table(){
   
   FILE *fp;
   int i,j,k;
   fp=fopen("5","r");
   for(i=0;i<6;i++){
      fscanf(fp,"%lf %lf %lf %lf %lf %lf",&tmpr_array[0][i][0],&tmpr_array[0][i][1],&tmpr_array[0][i][2],&tmpr_array[0][i][3],&tmpr_array[0][i][4],&tmpr_array[0][i][5]);
   }
   fclose(fp);
   
   fp=fopen("15","r");
   for(i=0;i<6;i++){
      fscanf(fp,"%lf %lf %lf %lf %lf %lf",&tmpr_array[1][i][0],&tmpr_array[1][i][1],&tmpr_array[1][i][2],&tmpr_array[1][i][3],&tmpr_array[1][i][4],&tmpr_array[1][i][5]);
   }
   fclose(fp);
   
   fp=fopen("25","r");
   for(i=0;i<6;i++){
      fscanf(fp,"%lf %lf %lf %lf %lf %lf",&tmpr_array[2][i][0],&tmpr_array[2][i][1],&tmpr_array[2][i][2],&tmpr_array[2][i][3],&tmpr_array[2][i][4],&tmpr_array[2][i][5]);
   }
   fclose(fp);
   
   fp=fopen("35","r");
   for(i=0;i<6;i++){
      fscanf(fp,"%lf %lf %lf %lf %lf %lf",&tmpr_array[3][i][0],&tmpr_array[3][i][1],&tmpr_array[3][i][2],&tmpr_array[3][i][3],&tmpr_array[3][i][4],&tmpr_array[3][i][5]);
   }
   fclose(fp);
   
   for(i=0;i<4;i++){
      for(j=0;j<6;j++){
         for(k=0;k<6;k++){
            printf("%lf ",tmpr_array[i][j][k]);
         }
         printf("\n");
      }
      printf("\n");
   }
  
}



void set_ratio(double t0,double t1,double base,int index){
   //ouble ratio_h,ratio_v;
   double diff0,diff1,diff;
   int i,j;
   int flag=0;
   diff0=t0-base;diff1=t1-base;
   if(diff0<diff1){diff=diff1;}
   else {diff=diff0;}
    
   for(i=0;i<6;i++){
      for(j=i;j<6;j++){
         if((tmpr_array[index][i][j]+diff)<90+0.5){ratio_h=i*0.1+0.5;ratio_v=j*0.1+0.5;flag=1;
         printf("set ratio debug: input: %f %f pred: %f   diff:%f  index:%d\n",t0,t1,tmpr_array[index][i][j]+diff,diff,index);
       }
       //  printf("entering");
         //else continue;
      }
   }
   
   if(flag==0) {ratio_h=0.5;ratio_v=0.5;}
   //ratio[0]=((ratio_v+0.5)/1.5)*((ratio_v+0.5)/1.5)*ratio_v;
   //ratio[1]=((ratio_h+0.5)/1.5)*((ratio_h+0.5)/1.5)*ratio_h;
   printf("%f  %f  %f\n",ratio_h,ratio_v,diff);
}

void init_standby_power(){
   int i;
   power_standby[0]=0.002;
   power_standby[1]=0.046;
   power_standby[2]=0;
   power_standby[3]=0.001;
   power_standby[4]=0.001;
   for(i=5;i<9;i++){
      power_standby[i]=0.025;
   }
   power_standby[9]=0.048;
   power_standby[10]=0.048;
   power_standby[11]=0.001;
   power_standby[12]=0.001;
   power_standby[13]=0;
   power_standby[14]=0.03;
   power_standby[15]=0;
   power_standby[16]=0.03;
   power_standby[17]=0.448/7;
   power_standby[18]=0.448/7*5;
   power_standby[19]=0.448/7;
   
   for(i=0;i<20;i++){
      power_standby[i]=power_standby[i]*2;
      power_standby[i+20]=power_standby[i];
   }
}

void readPower1(struct task *exe_queue){
  double pwr_tmp;
  int i;
  FILE *fp;
  //int speed[2];
  double pwr_in[20];
  //printf("ok in read1\n");
  fp=exe_queue->powerfile;
  
  if(feof(fp)!=1){
    fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&pwr_in[0],&pwr_in[1],&pwr_in[2],&pwr_in[3],&pwr_in[4],&pwr_in[5],&pwr_in[6],&pwr_in[7],&pwr_in[8],&pwr_in[9],&pwr_in[10],&pwr_in[11],&pwr_in[12],&pwr_in[13],&pwr_in[14],&pwr_in[15],&pwr_in[16],&pwr_in[17]);
    pwr_tmp=pwr_in[17];
    pwr_in[17]=pwr_tmp* 1.0 / 7.0;
    pwr_in[18]=pwr_tmp* 5.0 / 7.0;
    pwr_in[19]=pwr_tmp* 1.0 / 7.0;
    //printf("ok in read2\n");
    for(i=0;i<20;i++){
      pwr[i]+=pwr_in[i];
    }
  }
  exe_queue->powerfile=fp;
}

void readPower2(struct task *exe_queue){
  double pwr_tmp;
  int i;
  FILE *fp;
  //int speed[2];
  double pwr_in[20];

  fp=exe_queue->powerfile;
  if(feof(fp)!=1){
    fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&pwr_in[0],&pwr_in[1],&pwr_in[2],&pwr_in[3],&pwr_in[4],&pwr_in[5],&pwr_in[6],&pwr_in[7],&pwr_in[8],&pwr_in[9],&pwr_in[10],&pwr_in[11],&pwr_in[12],&pwr_in[13],&pwr_in[14],&pwr_in[15],&pwr_in[16],&pwr_in[17]);
    pwr_tmp=pwr_in[17];
    pwr_in[17]=pwr_tmp* 1.0 / 7.0;
    pwr_in[18]=pwr_tmp* 5.0 / 7.0;
    pwr_in[19]=pwr_tmp* 1.0 / 7.0;
    for(i=20;i<40;i++){
      pwr[i]+=pwr_in[i-20];
    }
  }
  exe_queue->powerfile=fp;
}
void init_rel_ratio(){
  FILE *fp;
  int i;
  printf("reading ratio\n");
  fp=fopen("rel_ratio_core_new","r");
  while(!feof(fp)){
    fscanf(fp,"%lf",&blk_rel_ratio[i]);
    i++;
    printf("%d\n",i);
  }
  printf("finish reading raito\n");
}
    
void readPower_core1(struct task *exe_queue){
  double pwr_tmp;
  int i;
  FILE *fp;
  //int speed[2];
  double pwr_in[20];
  //printf("ok in read1\n");
  fp=exe_queue->powerfile;
  
  if(feof(fp)!=1){
    fscanf(fp,"%lf ",&pwr_in[0]);
    
    //printf("ok in read2\n");
    
    pwr[0]+=pwr_in[0];
    
  }
  exe_queue->powerfile=fp;
}
void readPower_core2(struct task *exe_queue){
  double pwr_tmp;
  int i;
  FILE *fp;
  //int speed[2];
  double pwr_in[20];
  //printf("ok in read1\n");
  fp=exe_queue->powerfile;
  
  if(feof(fp)!=1){
    fscanf(fp,"%lf ",&pwr_in[0]);
    
    //printf("ok in read2\n");
    
    pwr[1]+=pwr_in[0];
    
  }
  exe_queue->powerfile=fp;
}
void power_input1(struct task *queue ,int speed){
  int i;
  if(queue->powerfile!=NULL){
    //printf("power file exist\n");
    if(queue->sim_exe>=speed){
        for(i=0;i<speed;i++){
            readPower_core1(queue);
        }
    }
    else{
      for(i=0;i<queue->sim_exe;i++){
        readPower_core1(queue);
      }
      if(queue->next!=NULL){
        power_input1(queue->next,speed-queue->sim_exe);
      }
    }
  }
  else{
    queue->powerfile=fopen(queue->name,"r");
    for(i=0;i<speed;i++){
      readPower_core1(queue);
    }
  }
}
void power_input2(struct task *queue ,int speed){
  int i;
  if(queue->powerfile!=NULL){
    //printf("power file exist\n");
    if(queue->sim_exe>=speed){
        for(i=0;i<speed;i++){
            readPower_core2(queue);
        }
    }
    else{
      for(i=0;i<queue->sim_exe;i++){
        readPower_core2(queue);
      }
      if(queue->next!=NULL){
        power_input2(queue->next,speed-queue->sim_exe);
      }
    }
  }
  else{
    queue->powerfile=fopen(queue->name,"r");
    for(i=0;i<speed;i++){
      readPower_core2(queue);
    }
  }
}

void power_input1_simple(struct task *queue ,int speed){
  int i;
  //if(queue->powerfile!=NULL){
    //printf("power file exist\n");
    if(queue->sim_exe>=speed){
        pwr[0]+=queue->avg_pwr*speed;
    }
    else{
      for(i=0;i<queue->sim_exe;i++){
        pwr[0]+=queue->avg_pwr*queue->sim_exe;
        //readPower_core1(queue);
      }
      if(queue->next!=NULL){
        power_input1_simple(queue->next,speed-queue->sim_exe);
      }
    }
  //}
  /*
  else{
    queue->powerfile=fopen(queue->name,"r");
    for(i=0;i<speed;i++){
      readPower_core1(queue);
    }
  }*/
}
void power_input2_simple(struct task *queue ,int speed){
  int i;
  //if(queue->powerfile!=NULL){
    //printf("power file exist\n");
    if(queue->sim_exe>=speed){
        pwr[1]+=queue->avg_pwr*speed;
    }
    else{
      for(i=0;i<queue->sim_exe;i++){
        pwr[1]+=queue->avg_pwr*queue->sim_exe;
        //readPower_core1(queue);
      }
      if(queue->next!=NULL){
        power_input2_simple(queue->next,speed-queue->sim_exe);
      }
    }
  //}
  /*
  else{
    queue->powerfile=fopen(queue->name,"r");
    for(i=0;i<speed;i++){
      readPower_core1(queue);
    }
  }*/
}