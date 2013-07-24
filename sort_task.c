#include <stdlib.h>
#include <stdio.h>
#include "sorttask.h"
#include <string.h>
#include "flp.h"
#include "RC.h"
#include <math.h>

#define LEN sizeof(struct task)

extern long double rel_bd[2]; //reliability due to oxide breakdown
extern long double rel_bd_pre[2];
//extern double temp_lamda_pre[20],temp_lamda2_pre[20];
extern long double rel_bd_blk[20],rel_bd_blk2[20];
extern long double rel_bd_blk_pre[20],rel_bd_blk2_pre[20];
//extern int n;
extern long double lamda_bd[2];
extern struct task tasklist[];
extern struct core real_c[2];
extern FILE *fp_lamda_bd_debug;
extern int global_out_count;
extern double blk_rel_ratio[17];
extern double eff_time[20];
extern double eff_time2[20];
#define MTTF 3e8

struct task *sort_deadline(struct task *sorted_queue,struct task *exe_queue){
  struct task *temp,*temp2,*previous;
  temp2=sorted_queue;
  //temp2->next=NULL;
  temp=exe_queue;
 
  while((temp->deadline<=temp2->deadline)&&(temp->next!=NULL)){
    previous=temp;
    temp=temp->next;
  }
  //printf("sort deadline\n");
  if(temp->deadline<=temp2->deadline){
    //printf("sort deadline2\n");
    temp2->next=temp->next;
    temp->next=temp2;
  }
 
  else{
    if(temp==exe_queue){
      temp2->next=temp;
      exe_queue=temp2;
    }
    else{
      temp2->next=temp;
      previous->next=temp2;
    }
  }
  //printf("sort deadline3\n");
  return exe_queue;
}


struct task *sort_task(struct task *initial_queue){

  struct task *head,*temp,*temp2,*mv_pointer,*previous;
  head=initial_queue;
  temp=initial_queue->next;
  head->next=NULL;
  while(temp!=NULL){
    mv_pointer=head;
    while((temp->arrival>=mv_pointer->arrival)&&(mv_pointer->next!=NULL)){
      previous=mv_pointer;
      mv_pointer=mv_pointer->next;
    }
    if(temp->arrival>=mv_pointer->arrival){
      temp2=temp;
      temp=temp->next;     
      temp2->next=mv_pointer->next; 
      mv_pointer->next=temp2;
    }
    else{
      temp2=temp;
      temp=temp->next;
      if(mv_pointer!=head){
	temp2->next=previous->next;
	previous->next=temp2;
      }
      else{
	temp2->next=head;
	head=temp2;
      }

    }
  }
  

  return head;
}
void change_assign(struct task *queue){
  struct task *temp;
  if(queue->assign==1){
    queue->assign=2;
    temp=real_c[0].current_task;
    real_c[0].current_task=real_c[1].current_task;
    real_c[1].current_task=temp;
  }
  else{
    queue->assign=1;
    temp=real_c[1].current_task;
    real_c[1].current_task=real_c[0].current_task;
    real_c[0].current_task=temp;
  }
}


int lowest_speed(struct task *queue,int sim_time){
  struct task *temp;
  int exe_remain=0;
  double speed_limit=0;
  double min_speed=0;
  //double speed_temp=0;
  int speed=0;
  temp=queue;
  while(temp!=NULL){
    exe_remain+=temp->sim_exe;
    speed_limit=(double) exe_remain/(temp->deadline-sim_time);
    if(min_speed<=speed_limit){
      min_speed=speed_limit;
    }
    if(temp->next==NULL){
      temp=NULL;
    }else{
      temp=temp->next;
    }
  }
  //printf("min_speed is %f",min_speed);
  while(speed<10*min_speed){
    speed+=1;
  }
  //printf("speed is %f",speed);
  return speed;
}

double lamda_cal(double temp[40+EXTRA],double speed,flp_t *flp,int core_n,long double lamda){
  double each_unit=0;
  int i;
  if(speed==0)speed=1;
  for(i=0;i<20;i++){
    each_unit+=pow(flp->units[i+core_n*20].width*flp->units[i+core_n*20].height,2)*expl(-0.95/(8.62e-5*temp[i+core_n*20]))/temp[i+core_n*20];
  }
  lamda+=each_unit*pow(speed,2)*1e8;//0.95 and 2 are for Al, according to JEDEC
  //printf("lamda:  %lf %lf %lf %lf %lf\n",lamda,flp->units[2+core_n*20].width,flp->units[2+core_n*20].height,exp(-100000/(8.62*temp[2+core_n*20]))*1e12,temp[2+core_n*20]);
  return lamda;
}

void lamda_cal_bd(double temp_lamda[20],double temp_lamda2[20],double temp_lamda_pre[20],double temp_lamda2_pre[20],double interval,double sim_time){
  int i;
  long double alpha[2],alpha_pre[2];
  long double temp;
  
  rel_bd_pre[0]=rel_bd[0];rel_bd_pre[1]=rel_bd[1];
  rel_bd[0]=1;rel_bd[1]=1;
  for(i=0;i<1;i++){
    //compute alpha
    /*
    alpha[0]=15*exp(-2.3*2.75)*exp(0.7/(temp_lamda[i]*8.62e-5));
    alpha[1]=15*exp(-2.3*2.75)*exp(0.7/(temp_lamda2[i]*8.62e-5));
    alpha_pre[0]=15*exp(-2.3*2.75)*exp(0.7/(temp_lamda_pre[i]*8.62e-5));
    alpha_pre[1]=15*exp(-2.3*2.75)*exp(0.7/(temp_lamda2_pre[i]*8.62e-5));//the only part matters is the latter, the former should be unchanged, only temperature changes!!!
    */
    alpha[0]=16e7*pow(1/1.1,78+0.081*temp_lamda[i])*exp((0.759-66.8/temp_lamda[i]-8.37e-4*temp_lamda[i])/(8.62*temp_lamda[i]*1e-5));
    alpha[1]=16e7*pow(1/1.1,78+0.081*temp_lamda2[i])*exp((0.759-66.8/temp_lamda2[i]-8.37e-4*temp_lamda2[i])/(8.62*temp_lamda2[i]*1e-5));
    alpha_pre[0]=16e7*pow(1/1.1,78+0.081*temp_lamda_pre[i])*exp((0.759-66.8/temp_lamda_pre[i]-8.37e-4*temp_lamda_pre[i])/(8.62*temp_lamda_pre[i]*1e-5));
    alpha_pre[1]=16e7*pow(1/1.1,78+0.081*temp_lamda2_pre[i])*exp((0.759-66.8/temp_lamda2_pre[i]-8.37e-4*temp_lamda2_pre[i])/(8.62*temp_lamda2_pre[i]*1e-5));


    //compute block reliability
    /*
    if(rel_bd_blk_pre[i]-exp(-pow((sim_time)*blk_rel_ratio[i]/alpha[0],2))<0){
      temp=0;
    }
    else{
      temp=rel_bd_blk_pre[i]-exp(-pow((sim_time)*blk_rel_ratio[i]/alpha[0],2));
    }
    temp=exp(-pow((sim_time-interval)*blk_rel_ratio[i]/alpha[0],2))-exp(-pow((sim_time)*blk_rel_ratio[i]/alpha[0],2));
    rel_bd_blk[i]=rel_bd_blk[i]*(1-temp);
    
    if(rel_bd_blk2_pre[i]-exp(-pow((sim_time)*blk_rel_ratio[i]/alpha[1],2))<0){
      temp=0;
    }
    else{
      temp=rel_bd_blk2_pre[i]-exp(-pow((sim_time)*blk_rel_ratio[i]/alpha[1],2));
    }
    temp=exp(-pow((sim_time-interval)*blk_rel_ratio[i]/alpha[1],2))-exp(-pow((sim_time)*blk_rel_ratio[i]/alpha[1],2));
    rel_bd_blk2[i]=rel_bd_blk2[i]*(1-temp);
    */
    /*
    if((exp(-pow((sim_time-interval)*blk_rel_ratio[i]/alpha_pre[0],2))-exp(-pow(sim_time*blk_rel_ratio[i]/alpha[0],2)))<0){
      temp=0;
    }
    else{
      temp=exp(-pow((sim_time-interval)*blk_rel_ratio[i]/alpha_pre[0],2))-exp(-pow(sim_time*blk_rel_ratio[i]/alpha[0],2));
    }
    rel_bd_blk[i]=rel_bd_blk[i]*(1-temp);
    
    if((exp(-pow((sim_time-interval)*blk_rel_ratio[i]/alpha_pre[1],2))-exp(-pow(sim_time*blk_rel_ratio[i]/alpha[1],2)))<0){
      temp=0;
    }
    else{
      temp=exp(-pow((sim_time-interval)*blk_rel_ratio[i]/alpha_pre[1],2))-exp(-pow(sim_time*blk_rel_ratio[i]/alpha[1],2));
    }
    rel_bd_blk2[i]=rel_bd_blk2[i]*(1-temp);
    */
    //temp=exp(-pow((sim_time+1e7-interval)*blk_rel_ratio[i]/alpha[0],2))-exp(-pow((sim_time+1e7)*blk_rel_ratio[i]/alpha[0],2));
    //temp=exp(-pow((sim_time-interval+1.2e7)/alpha[0],2))-exp(-pow((sim_time+1.2e7)/alpha[0],2));
    
    

    temp=exp(-pow((sim_time-interval)/alpha[0],2))-exp(-pow(sim_time/alpha[0],2));
    rel_bd_blk[i]=rel_bd_blk[i]*(1-temp);
    temp=exp(-pow((sim_time-interval)/alpha[1],2))-exp(-pow(sim_time/alpha[1],2));
    rel_bd_blk2[i]=rel_bd_blk2[i]*(1-temp);
    
    /*
    temp=exp(-(sim_time-interval)/alpha[0])-exp(-sim_time/alpha[0]);
    rel_bd_blk[i]=rel_bd_blk[i]*(1-temp);
    temp=exp(-(sim_time-interval)/alpha[1])-exp(-sim_time/alpha[1]);
    rel_bd_blk2[i]=rel_bd_blk2[i]*(1-temp);
    */
  //compute chip reliability
    /*
    eff_time[i]+=MTTF/alpha[0]*interval;
    eff_time2[i]+=MTTF/alpha[1]*interval;
    rel_bd_blk[i]=exp(-pow(eff_time[i]/MTTF,2));
    rel_bd_blk2[i]=exp(-pow(eff_time2[i]/MTTF,2));
    */
    //rel_bd_blk[i]=exp(-eff_time[i]/MTTF);
    //rel_bd_blk2[i]=exp(-eff_time2[i]/MTTF);  
    
    //rel_bd_blk_pre[i]=rel_bd_blk[i];
    //rel_bd_blk2_pre[i]=rel_bd_blk2[i];
    rel_bd[0]*=rel_bd_blk[i];
    rel_bd[1]*=rel_bd_blk2[i];

    
    //printf("alpha[0] %f rel_bd_blk %Le    ",alpha[0],rel_bd_blk[i]);
    
  }
  //printf("rel_bd[0] %Le rel_bd_pre[0] %Le  lamda_bd  %Le   rel_pre[0] %Le rel_pre[1] %Le sample blk_rel:%Le sample temp: %f\n",rel_bd[0],rel_bd_pre[0],(rel_bd_pre[0]-rel_bd[0])/interval/rel_bd[0],rel_bd[1],rel_bd_pre[1],rel_bd_blk[2],temp_lamda[0]);
  //compute lamda
  //lamda_bd[0]=(rel_bd_pre[0]-rel_bd[0])/interval/rel_bd[0];
  //lamda_bd[1]=(rel_bd_pre[1]-rel_bd[1])/interval/rel_bd[1];
 
   //fprintf(fp,"%f %f simtime:  %f  rel_bd[0]: %Le  pre: %Le  rel_bd[1]: %Le  pre:%Le total:  %Le\n",temp_lamda[0],temp_lamda2[0],sim_time,rel_bd[0],rel_bd_pre[0],rel_bd[1],rel_bd_pre[1],rel_bd[0]*rel_bd[1]);
   
 //return fp;
}

void q_sort(struct task array[],int size){
  quick_sort(array,0,size);
}
void quick_sort(struct task array[],int left,int right){
  int pivot;
  if(left<right){
    pivot=partition(array,left,right);
    quick_sort(array,left,pivot);
    quick_sort(array,pivot+1,right);
  }
}
int partition(struct task array[],int left,int right){
  int i,j;
  double p;
  i=left;
  p=array[left].utilization;
  for(j=left+1;j<right;j++){
    if(array[j].utilization>p){
      i++;
      swap(array,i,j);
    }
  }
  swap(array,left,i);
  return i;
}
void swap(struct task array[],int i, int j){
  struct task temp;
  temp=array[i];
  array[i]=array[j];
  array[j]=temp;
}





void q_sort_c(struct core array[],int size){
  quick_sort_c(array,0,size);
}
void quick_sort_c(struct core array[],int left,int right){
  int pivot;
  if(left<right){
    pivot=partition_c(array,left,right);
    quick_sort_c(array,left,pivot);
    quick_sort_c(array,pivot+1,right);
  }
}
int partition_c(struct core array[],int left,int right){
  int i,j;
  double p;
  i=left;
  p=array[left].utilization;
  for(j=left+1;j<right;j++){
    if(array[j].utilization<p){
      i++;
      swap_c(array,i,j);
    }
  }
  swap_c(array,left,i);
  return i;
}
void swap_c(struct core array[],int i, int j){
  struct core temp;
  temp=array[i];
  array[i]=array[j];
  array[j]=temp;
}