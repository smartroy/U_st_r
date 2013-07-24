#include <stdlib.h>
#include <stdio.h>
#include "sorttask.h"
#include <string.h>
#include "RC.h"
#include "flp.h"
#include "temppred.h"
#include <math.h>


#define OUT_COUNT 100
#define TOTAL_MODULES 18
#define MODULE_DECODE 0
#define MODULE_BRANCH 1
#define MODULE_RAT 2
#define MODULE_RUU 3
#define MODULE_LSQ 4
#define MODULE_IALU1 5
#define MODULE_IALU2 6
#define MODULE_IALU3 7
#define MODULE_IALU4 8
#define MODULE_FALU1 9
#define MODULE_FALU2 10
#define MODULE_INTREG 11
#define MODULE_FPREG 12
#define MODULE_ITLB 13
#define MODULE_IL1 14
#define MODULE_DTLB 15
#define MODULE_DL1 16
#define MODULE_L2 17
#define MODULE_L2_LEFT MODULE_L2
#define MODULE_L2_BOTTOM MODULE_L2 + 1
#define MODULE_L2_RIGHT MODULE_L2 + 2
#define MODULE_L2_EXTRA 2

#define TIME_SLICE 1e-1
#define SCH_COUNT 30000
#define INTERVAL 10000
#define U 1e-2
#define T_MAX 90

#define R_H_M 0.78
#define R_M_L 0.765
#define GAP_T 20
#define LEN sizeof(struct task)
#define FINE 1e-5
#define COARSE 3
#define DEV 5
struct task *sorted_queue,*sorted_queue_base;
struct core real_c[2];
flp_t *flp;
double pwr[2]; // +2 because L2 include 3 parts: left, right and bottom
double tmpr[2+EXTRA]; // +2 because L2 include 3 parts: left, right and bottom 
double pwr_vec[2+EXTRA];
double power_standby[40];

double tmpr_array[4][6][6];
double ratio[2];//power ratio
double ratio_v,ratio_h;// frequency and voltage ratio
FILE *fp_out,*fp_out2,*fp_out4,*fp_standby,*fp_previous,*fp_out5,*fp_lamda,*fp_lamda_sum,*fp_lamda_bd,*fp_lamda_bd_debug,*fp_schedule;
int end0,end1,one_task;
int heat_0=0;
int heat_1=0;
long double rel_bd[2]; //reliability due to oxide breakdown
long double rel_bd_pre[2];
long double rel_bd_blk[20],rel_bd_blk2[20];
long double rel_bd_blk_pre[20],rel_bd_blk2_pre[20];
double temp_lamda[20],temp_lamda2[20],temp_lamda_pre[20],temp_lamda2_pre[20];
long double lamda_bd[2];
int global_out_count;
double A[4096],B[4096];
double blk_rel_ratio[17];
double eff_time[20],eff_time2[20];
double core_U[2];
struct queue_util{
  unsigned int deadline;
  double utilization;
  int temporary;
  struct queue_util *next;
}queue_util;
int exe_mean[4];
int sigma;

struct task tasklist[]={
  {"mpeg2enc_core",0,80,200,0,101,0,80,5000,0,0,0,0,0.09403519,NULL,NULL},
  {"h263enc_core",0,100,200,0,85,0,100,10000,0,1,0,0,0.09420586,NULL,NULL},
  //{"susan_relax_scale",0,40,100,0,101,0,40,5000,0,2,0,0,NULL,NULL},
  {"basicmath_core",0,100,200,0,101,0,100,5000,0,2,0,0,0.07972,NULL,NULL},
  {"patricia_core",0,50,100,0,101,0,50,5000,0,3,0,0,0.0657524,NULL,NULL},
  // {"cc1_64kl13g",0,0.2,1342,250,0,4,NULL,NULL},
  //{"di__cbc_e",0,0.12,0.000039,0,103.5,0,0.12,1,0,2,NULL,NULL},
  // {"mpeg_d_15_20",0,0.2,1342,30,0,5,NULL,NULL},
  // {"bzip_cbc",5,50.0,21.47484,0,101,0,50.0,2,0,5,NULL,NULL},
   //{"bzip_cbc",0.2,0.5,0.2,0,101,0,0.5,2,0,5,NULL,NULL},
};

int n;
struct task *add_task(int period_lcm);
void init();
int get_lcm(int sub_lcm, int input);
struct task *form_list(struct task *sorted,double lcm);
struct task *form_list(struct task *sorted,double lcm);
struct queue_util *update_util(unsigned int sim_time,struct queue_util *queue,int number);
struct queue_util *add_u(unsigned long deadline,double utilization,struct queue_util *queue,int temporary);
double cal_u(struct queue_util *queue);
void print_matrix(double *a,FILE *fp,int m,int n);

main(int argc,char *argv[]){
  struct task *initial_queue,*exe_queue,*exe_queue2,*task_1,*task_2,*q_temp;
  struct task *q_pointer,*q_pointer2,*q_free;//for debug & schedule test
  struct queue_util *q_u1,*q_u2; //queue utilization list
  unsigned int sim_time=0;
  int all_task_in=0;
  int available_finished=0;
  int offset;
  int speed[2];
  unsigned int output_count=0;
  int i,omit_lateral;
  int last_scaled=0;//indicator for limitting scaled too often
  int speeding_up=0;
  int freq_scale=0;//to indicate the freq is scaling,bring in penalty
  double THRESH=70;
  //double speed_limit[2];//for schedule test
  double min_speed;//for schedule test
  int no_of_task=0;//for schedule test
  double exe_remain[2];//for schedule test
  double speed_temp=0;
  int speed_limit[2];
  speed_limit[0]=0;speed_limit[1]=0;

  float temp_pred[2];
  double tmax[2];
  double tmax_core;
  int schedule_count=0;
  long double lamda[2],rel[2];//for lamda and reliability calculating
  int lamda_count,rel_count;
  int lamda_bd_count=0;
  int lamda_bd_count_out=0;
  int task_num=sizeof(tasklist)/sizeof(struct task);
  int lcm=0;
  unsigned long round=2500000;
  int new_task[2];//indicating if new iteration of task enter queue;
  int debug_count=0;
  int debug_count2=0;
  new_task[0]=0;new_task[1]=0;
  int schedule_switch=0;
  ratio[0]=1;ratio[1]=1;ratio_v=1;ratio_h=1;
  int swap_f[4];
  int swap_counter[4];
  for(i=0;i<4;i++){
   swap_f[i]=0;
   swap_counter[i]=0;
  }
  unsigned int swap_total=0;
  double tmpr_acc=0;
  double swap_threshold=FINE;//50;

  FILE *fp_input_ratio;
  int input_u_ratio;
  FILE *fp_input_mean;
  FILE *fp_input_sigma;

//tilts
  int tilts_count=0;
  FILE *fp_A,*fp_B;
  fp_A=fopen("matrix_A","w");
  fp_B=fopen("matrix_B","w");
  if(argc!=4){
    printf("using default U ratio\n");
    input_u_ratio=10;
    for(i=0;i<sizeof(tasklist)/sizeof(struct task);i++){
      exe_mean[i]=5;
    }
    sigma=DEV;
  }
  else{
    for(i=0;i<argc;i++){
      printf("%s\n",argv[i]);
    }
    fp_input_ratio=fopen(argv[1],"r");
    fscanf(fp_input_ratio,"%d",&input_u_ratio);
    printf("U ratio %d\n",input_u_ratio);
    fp_input_mean=fopen(argv[2],"r");
    for(i=0;i<sizeof(tasklist)/sizeof(struct task);i++){
      fscanf(fp_input_mean,"%d",&exe_mean[i]);
      printf("mean exe %d\n",exe_mean[i]);
    }
    fp_input_sigma=fopen(argv[3],"r");
    fscanf(fp_input_sigma,"%d",&sigma);
    printf("sigma %d\n",sigma);
  }
  for(i=0;i<sizeof(tasklist)/sizeof(struct task);i++){
    tasklist[i].deadline=tasklist[i].deadline/10*input_u_ratio;
    tasklist[i].period=tasklist[i].period/10*input_u_ratio;
  }
  for(i=0;i<sizeof(tasklist)/sizeof(struct task);i++){
    printf("deadline  %d  period  %d\n",tasklist[i].deadline,tasklist[i].period);
  }


  exe_queue=NULL;
  exe_queue2=NULL;
  q_u1=NULL;
  q_u2=NULL; 


  tmpr_table();
  flp=read_flp("dual_core_only.flp");
  fp_out=fopen("output_p","w");
  create_RC_matrices(flp,omit_lateral);
  init();
  printf("opening standby_power\n");
  init_standby_power();
  init_rel_ratio();
 
  
  //use QSORT to order task and assign task
  q_sort(tasklist,sizeof(tasklist)/sizeof(struct task));
  for(i=0;i<sizeof(tasklist)/sizeof(struct task);i++){
    printf("%f   %d\n",tasklist[i].utilization,tasklist[i].number);
  }
  for(i=0;i<sizeof(tasklist)/sizeof(struct task);i++){
    q_sort_c(real_c,sizeof(real_c)/sizeof(struct core));
    swap_f[tasklist[i].number]=real_c[0].number;
    real_c[0].utilization+=tasklist[i].utilization;
    
  }
  for(i=0;i<sizeof(tasklist)/sizeof(struct task);i++){
    printf("%d\n",swap_f[tasklist[i].number]);
  }
   
  for(i=0;i<2;i++){
    core_U[real_c[i].number-1]=real_c[i].utilization;
    real_c[i].utilization=0;
  }
  for(i=0;i<2;i++){
    real_c[i].number=i+1;
    printf("%f   ",core_U[i]);
  }
  printf("\n");

  //return 0;
	 
   

//tilts
  /*
  generate_AB(51,40,1e-6);
  double base_intvl = 1e-6;
  double times = 1e-1 / base_intvl;
  int exp2 = 2; int num_exp2 = 1;
  while (exp2 < times) {
    exp2 *= 2; num_exp2 ++;
  }
  //generate_AB(51,40,1e-1);
  for(i=0;i<num_exp2-1;i++){
    calculate_2timeAB(A,B,51,40);
  }
  print_matrix(A,fp_A,51,51);
  print_matrix(B,fp_B,51,40);
  init();*/
  

  tmax[1]=T_INIT;tmax[0]=T_INIT;
  real_c[0].speed=0;
  real_c[1].speed=0;
  speed[0]=0;speed[1]=0;
  fp_out2=fopen("output_t","w");
  fp_out4=fopen("output_debug","w");
  fp_out5=fopen("out_error","w");
  fp_lamda=fopen("out_lamda","w");
  fp_lamda_sum=fopen("out_lamda_sum","w");
  fp_lamda_bd=fopen("lamda_bd","w");
  fp_lamda_bd_debug=fopen("lamda_bd_debug","w");
  fp_schedule=fopen("schedule_point","w");

  lamda[0]=0;lamda[1]=0;rel[0]=0;rel[1]=0;
  lamda_count=0;rel_count=0;
  //get lcm of periods
  for(i=0;i<task_num;i++){
    printf("period: %d\n",tasklist[i].period);
    lcm=get_lcm(lcm,tasklist[i].period);
    printf("lcm: %d\n",lcm);
  }
  initial_queue=add_task(lcm);
  q_pointer=initial_queue;
  while(q_pointer!=NULL){
    printf("%d  ",q_pointer->arrival);
    q_pointer=q_pointer->next;
  }
  printf("\n");
  srand(time(NULL));
  sorted_queue_base=sort_task(initial_queue);
  q_pointer=sorted_queue_base;
  while(q_pointer!=NULL){
    printf("%d  ",q_pointer->number);
    q_pointer=q_pointer->next;
  }
  printf("\n");
  sorted_queue=form_list(sorted_queue_base,0);
  //return 0;
   
  while(!(all_task_in&&available_finished)){
   
    //add task into exe_q according to arrival time and sort according to deadline
    if(all_task_in==0){
      while(sorted_queue!=NULL && (sorted_queue->arrival==sim_time)&&(sorted_queue->next!=NULL)){
	//printf("insert task\n");
	if(sorted_queue!=NULL){
	  //printf("first queue  %d \n",sorted_queue->number);
	}else{
	  //printf("sort queue null\n");
	}
	//printf("no problem2 %d %d\n",sorted_queue->arrival,sim_time);
	q_u1=update_util(sim_time,q_u1,0);
	q_u2=update_util(sim_time,q_u2,1);
	real_c[0].utilization=cal_u(q_u1);
	real_c[1].utilization=cal_u(q_u2);
	
	if(schedule_switch==1){
	  
	  if((rel_bd[0]-rel_bd[1])>0){   //task allocation when entering system
	    //printf("state 1\n");
	    if(swap_f[sorted_queue->number]==1){
	      real_c[0].utilization+=sorted_queue->utilization;
	      q_u1=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u1,0);
	      //printf("state 1 ok1\n");
	      exe_queue=send_to_exe_queue(exe_queue);
	    }
	    else{
	      if(core_U[0]+sorted_queue->utilization<=1){
		core_U[0]+=sorted_queue->utilization;
		sorted_queue->temporary=1;
		real_c[0].utilization+=sorted_queue->utilization;
		q_u1=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u1,1);
		//printf("state 1 ok1\n");
		exe_queue=send_to_exe_queue(exe_queue);
	      }
	      else{
		real_c[1].utilization+=sorted_queue->utilization;
		q_u2=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u2,0);
		// printf("ok 1\n");
		exe_queue2=send_to_exe_queue(exe_queue2);
	      }
	    }
	  }
	 
	  else {
	    //printf("state 2\n");
	    if(swap_f[sorted_queue->number]==2){
	      real_c[1].utilization+=sorted_queue->utilization;
	      q_u2=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u2,0);
	      // printf("ok 1\n");
	      exe_queue2=send_to_exe_queue(exe_queue2);
	    }
	    else{
	      if(core_U[1]+sorted_queue->utilization<=1){
		core_U[1]+=sorted_queue->utilization;
		sorted_queue->temporary=1;
		real_c[1].utilization+=sorted_queue->utilization;
		q_u2=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u2,1);
		exe_queue2=send_to_exe_queue(exe_queue2);
	      }
	      else{
		real_c[0].utilization+=sorted_queue->utilization;
		q_u1=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u1,0);
		// printf("state 1 ok1\n");
		exe_queue=send_to_exe_queue(exe_queue);
	      }
	    }
		
	  }
	    
	    //  printf("ok 2\n") 
	  
	  
	}
	else{
	  if(swap_f[sorted_queue->number]==1){
	    real_c[0].utilization+=sorted_queue->utilization;
	    q_u1=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u1,0);
	    exe_queue=send_to_exe_queue(exe_queue);
	  }
	  else{
	    real_c[1].utilization+=sorted_queue->utilization;
	    q_u2=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u2,0);
	    exe_queue2=send_to_exe_queue(exe_queue2);
	  }
	}
	
	//fprintf(fp_out5,"send to core1\n");
	/*
	if(exe_queue==NULL){speed[0]=0;}
	else{
	  speed[0]=1;
	  while(speed[0]<real_c[0].utilization*10){
	    speed[0]++;
	  }
	}
	if(exe_queue2==NULL){speed[1]=0;}
	else{
	  speed[1]=1;
	  while(speed[1]<real_c[1].utilization*10){
	    speed[1]++;
	  }
	}
   */
	if(exe_queue==NULL){speed[0]=0;}
	else{speed[0]=10;}
	if(exe_queue2==NULL){speed[1]=0;}
	else{speed[1]=10;}
	//fprintf(fp_out4,"simtime: %d  speed[0]: %d %f %d  speed[1]: %d %f %d   %d\n",sim_time,speed[0],real_c[0].utilization,speed_limit[0],speed[1],real_c[1].utilization,speed_limit[1],swap_total);
	
	
      }
      
      // printf("here ok upper6\n");//dealing with the last task in sorted_queue(sort according to arrival time)
      if(sorted_queue!=NULL && sorted_queue->arrival<=sim_time){
	q_u1=update_util(sim_time,q_u1,0);
	q_u2=update_util(sim_time,q_u2,1);
	real_c[0].utilization=cal_u(q_u1);
	real_c[1].utilization=cal_u(q_u2);
	if(schedule_switch==1){
	  if((rel_bd[0]-rel_bd[1])>0){   //task allocation when entering system
	    //printf("state 1\n");
	    if(swap_f[sorted_queue->number]==1){
	      real_c[0].utilization+=sorted_queue->utilization;
	      q_u1=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u1,0);
	      //printf("state 1 ok1\n");
	      exe_queue=send_to_exe_queue(exe_queue);
	    }
	    else{
	      if(core_U[0]+sorted_queue->utilization<=1){
		core_U[0]+=sorted_queue->utilization;
		sorted_queue->temporary=1;
		real_c[0].utilization+=sorted_queue->utilization;
		q_u1=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u1,1);
		//printf("state 1 ok1\n");
		exe_queue=send_to_exe_queue(exe_queue);
	      }
	      else{
		real_c[1].utilization+=sorted_queue->utilization;
		q_u2=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u2,0);
		// printf("ok 1\n");
		exe_queue2=send_to_exe_queue(exe_queue2);
	      }
	    }
	  }
	  else {
	    //printf("state 2\n");
	    if(swap_f[sorted_queue->number]==2){
	      real_c[1].utilization+=sorted_queue->utilization;
	      q_u2=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u2,0);
	      // printf("ok 1\n");
	      exe_queue2=send_to_exe_queue(exe_queue2);
	    }
	    else{
	      if(core_U[1]+sorted_queue->utilization<=1){
		core_U[1]+=sorted_queue->utilization;
		sorted_queue->temporary=1;
		real_c[1].utilization+=sorted_queue->utilization;
		q_u2=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u2,1);
		exe_queue2=send_to_exe_queue(exe_queue2);
	      }
	      else{
		real_c[0].utilization+=sorted_queue->utilization;
		q_u1=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u1,0);
		// printf("state 1 ok1\n");
		exe_queue=send_to_exe_queue(exe_queue);
	      }
	    }
		
	  }
	  
	  
	}
	else{
	  if(swap_f[sorted_queue->number]==1){
	    real_c[0].utilization+=sorted_queue->utilization;
	    q_u1=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u1,0);
	    exe_queue=send_to_exe_queue(exe_queue);
	  }
	  else{
	    real_c[1].utilization+=sorted_queue->utilization;
	    q_u2=add_u(sorted_queue->deadline,sorted_queue->utilization,q_u2,0);
	    exe_queue2=send_to_exe_queue(exe_queue2);
	  }
	  
	}
	/*
	if(exe_queue==NULL){speed[0]=0;}
	else{speed[0]=10;}
	if(exe_queue2==NULL){speed[1]=0;}
	else{speed[1]=10;}
	*/
  /*
	if(exe_queue==NULL){speed[0]=0;}
	else{
	  speed[0]=1;
	  while(speed[0]<real_c[0].utilization*10){
	    speed[0]++;
	  }
	}
	if(exe_queue2==NULL){speed[1]=0;}
	else{
	  speed[1]=1;
	  while(speed[1]<real_c[1].utilization*10){
	    speed[1]++;
	  }
	}
   */
	if(exe_queue==NULL){speed[0]=0;}
	else{speed[0]=10;}
	if(exe_queue2==NULL){speed[1]=0;}
	else{speed[1]=10;}
	
	//if(exe_queue2!=NULL){
	//fprintf(fp_out4,"simtime: %d  speed[0]: %d %f %d  speed[1]: %d %f %d  %d\n",sim_time,speed[0],real_c[0].utilization,speed_limit[0],speed[1],real_c[1].utilization,speed_limit[1],swap_total);
	  //}
	  //else{fprintf(fp_out4,"core2 idle\n");}
      }
      //fprintf(fp_out4,"simtime: %d  %f  %f %d\n",sim_time,real_c[0].utilization,real_c[1].utilization,schedule_switch);

    }

    if(sorted_queue==NULL && round!=0){
      //printf("add new lcm of tasks\n");
      sorted_queue=form_list(sorted_queue_base,lcm);
      q_pointer=sorted_queue_base;
      while(q_pointer!=NULL){
	q_pointer->arrival+=lcm;
	q_pointer->deadline+=lcm;
	q_pointer=q_pointer->next;
      }
      round--;
    }
    
    //printf("cycle start\n");
      //fprintf(fp_out5,"1speed0,%d speed1 %d\n",speed[0],speed[1]);
    //set speed
    real_c[0].speed=speed[0];
    real_c[1].speed=speed[1];
    new_task[0]=0;
    new_task[1]=0;
    // fprintf(fp_out5,"2speed0,%d speed1 %d\n",speed[0],speed[1]);
    //printf("ok 3\n");
    //read power profile
    if(exe_queue!=NULL){
      exe_queue->assign=1;
    }
    if(exe_queue2!=NULL){
      exe_queue2->assign=1;
    }
   // printf("read power start\n");
   
    if(exe_queue!=NULL){
      power_input1(exe_queue,speed[0]);
    }
    else{
      //for(i=0;i<20;i++){
      pwr[0]=0;//power_standby[0]/100;
      //}
    }
    //printf("read power 2\n");
    if(exe_queue2!=NULL){
      power_input2(exe_queue2,speed[1]);
    }
    else{
      //for(i=0;i<20;i++){
      pwr[1]=0;//power_standby[0]/100;
      //}
    }
    //printf("read power finished\n");
    //fprintf(fp_out5,"5speed0,%d speed1 %d\n",speed[0],speed[1]);
    
    
    
    //dealing with unit issue
    //   printf("ok 4\n");
    for(i=0;i<2;i++){
      pwr[i]=pwr[i]/1e-1;//uJ, time_slice=1e-5
      //pwr[i+20]=pwr[i+20]/1e-1;
    }
    for(i=0;i<2;i++){
      pwr_vec[i]=pwr[i];
      //pwr_vec[i+20]=pwr[i+20];
    }		    
    /*
    for(i=0;i<20;i++){
      pwr_vec[i]=0;
      pwr_vec[i+20]=0;
    }*/

    //send power value into hotspot power vector
    
    
   // printf("no problem4\n");
    //printf("temperature: %f %f  %d  %d\n",tmpr[0],tmpr[1],speed[0],speed[1]);
    //printf("%f   %f\n",pwr_vec[0],pwr_vec[1]);
    compute_temp_simple(pwr_vec,tmpr,1e-1);
    //printf("temperature: %f %f %d\n",tmpr[0],tmpr[1],sim_time+1);
    //return 0;
    /*if(exe_queue!=NULL){
    	printf("%d",exe_queue->number);
    }
    else{
    	printf("exe_queue empty");
    }
    if(sim_time>10000){
    	return 0;
    }*/
    

    tmax[0]=tmpr[0];
    tmax[1]=tmpr[1];
    /*
    for(i=1;i<20;i++){
      if(tmpr[i]>tmax[0]){tmax[0]=tmpr[i];}
      if(tmpr[20+i]>tmax[1]){tmax[1]=tmpr[i+20];}
    }*/
    if(tmax[0]>=tmax[1]){tmax_core=tmax[0];}
    else{tmax_core=tmax[1];}
    real_c[0].tmpr=tmax[0]-273.15;
    real_c[1].tmpr=tmax[1]-273.15;
    

    
    if(lamda_bd_count<2){
      for(i=0;i<1;i++){
	temp_lamda[i]+=tmpr[0];
	temp_lamda2[i]+=tmpr[1];
      }
      lamda_bd_count++;
      //printf("%d  %f   %f\n",lamda_bd_count,temp_lamda[1],temp_lamda2[1]);
    }
    else{
      
      for(i=0;i<1;i++){
	temp_lamda[i]/=lamda_bd_count;
	temp_lamda2[i]/=lamda_bd_count;
      }
      rel_bd_pre[0]=rel_bd[0];
      rel_bd_pre[1]=rel_bd[1];
      //return 0;
      lamda_cal_bd(temp_lamda,temp_lamda2,temp_lamda_pre,temp_lamda2_pre,0.1*lamda_bd_count,(double)sim_time/10);
      //fprintf(fp_lamda_bd,"lamda_bd_1: %Le    lamda_bd_2:  %Le   simtime:  %f\n",lamda_bd[0],lamda_bd[1],sim_time*1e5);
      //printf("rel_bd[0] %Le rel_bd[1] %Le rel_pre[0] %Le rel_pre[1] %Le\n",rel_bd[0],rel_bd[1],rel_bd_pre[0],rel_bd_pre[1]);
      lamda_bd_count=0;
      for(i=0;i<1;i++){
	// temp_lamda_pre[i]=temp_lamda[i];//save previous interval temperature for obd cal
	//     temp_lamda2_pre[i]=temp_lamda2[i];
	     temp_lamda[i]=0;
	     temp_lamda2[i]=0;
      }
      
      //tmpr_acc+=(tmax[0]-tmax[1])*0.01;
      tmpr_acc+=(rel_bd_pre[0]-rel_bd[0])-(rel_bd_pre[1]-rel_bd[1]);
      if(schedule_switch==0){
	if(fabs(tmpr_acc)>swap_threshold){
	  fprintf(fp_lamda,"%f T1: %f  T2: %f  U1: %f  U2:  %f  %f ",(double)sim_time/10.0,tmax[0],tmax[1],core_U[0],core_U[1],tmpr_acc);
	  if(tmpr_acc>0){
	    for(i=0;i<sizeof(tasklist)/sizeof(struct task);i++){
	      if(swap_f[tasklist[i].number]==1&&core_U[1]+tasklist[i].utilization<=1){
		swap_f[tasklist[i].number]=2;
		core_U[1]+=tasklist[i].utilization;
		core_U[0]-=tasklist[i].utilization;
		break;
	      }
	    }
	  }
	  else{
	    for(i=0;i<sizeof(tasklist)/sizeof(struct task);i++){
	      if(swap_f[tasklist[i].number]==2&&core_U[0]+tasklist[i].utilization<=1){
		swap_f[tasklist[i].number]=1;
		core_U[0]+=tasklist[i].utilization;
		core_U[1]-=tasklist[i].utilization;
		break;
	      }
	    }
	  }
	tmpr_acc=0;
	fprintf(fp_lamda,"U1_n:  %f  U2_n:  %f\n",core_U[0],core_U[1]);
	//return 0;
	}
      }
      //return 0;
    }
  if(lamda_bd_count_out<3600){
    lamda_bd_count_out++;
  } 
  else{
    fprintf(fp_lamda_bd,"%f %f simtime:  %f  rel_bd[0]: %Le  pre: %Le  rel_bd[1]: %Le  pre:%Le total:  %Le\n",eff_time[0],eff_time2[0],(double)sim_time/10,rel_bd[0],rel_bd_pre[0],rel_bd[1],rel_bd_pre[1],rel_bd[0]*rel_bd[1]);
    lamda_bd_count_out=0;
  }
    //printf("ok lamda_bd2   %Le   %Le\n",rel_bd[0],rel_bd_pre[0]);   

    //printf("removing finished tasks\n");
    //remove finished task
    if(exe_queue!=NULL){
      //printf("sim_exe %d\n",exe_queue->sim_exe);
      if(exe_queue->sim_exe>speed[0]){
	//	printf("sim_exe %d 01\n",exe_queue->sim_exe);
	exe_queue->sim_exe=exe_queue->sim_exe-speed[0];
	//debug_count+=speed[0];
	exe_queue->execution+=speed[0];
      }
      //printf("sim_exe %f\n",exe_queue->sim_exe);
      else{
	//printf("sim_exe %d 02\n",exe_queue->sim_exe);
	//real_c[0].utilization-=exe_queue->utilization;
	//	printf("sim_exe %d 03\n",exe_queue->sim_exe);
	if(exe_queue->next==NULL){
	  fclose(exe_queue->powerfile);
	  q_free=exe_queue;
	  exe_queue=NULL;
	  free(q_free);
	  speed_limit[0]=0;
	    
	}else{
	  fclose(exe_queue->powerfile);
	  exe_queue->next->sim_exe-=(speed[0]-exe_queue->sim_exe);
	  q_free=exe_queue;
	  exe_queue=exe_queue->next;
	  free(q_free);
	  
	}
      }
      //printf("sim_exe %d 2\n",exe_queue->sim_exe);
    }
    if(exe_queue2!=NULL){
      //printf("sim_exe %d\n",exe_queue2->sim_exe);
      if(exe_queue2->sim_exe>speed[1]){
	exe_queue2->sim_exe=exe_queue2->sim_exe-speed[1];
	debug_count+=speed[1];
	exe_queue2->execution+=speed[1];
      }
      //printf("sim_exe %f\n",exe_queue->sim_exe);
      else{
	//real_c[1].utilization-=exe_queue2->utilization;
	if(exe_queue2->next==NULL){
	  fclose(exe_queue2->powerfile);
	  q_free=exe_queue2;
	  exe_queue2=NULL;
	  free(q_free);
	  //speed_limit[1]=0;
	  // printf("exe_queue NULL\n");
	}else{
	  fclose(exe_queue2->powerfile);
	  exe_queue2->next->sim_exe-=(speed[1]-exe_queue2->sim_exe);
	  q_free=exe_queue2;
	  exe_queue2=exe_queue2->next;
	  free(q_free);
	  
	}
      }
    }
    
    //printf("no problem5\n");
    
    
    //printf("ok 8\n");
    
    sim_time+=1;
    //printf("sim_time: %d\n",sim_time);
    speed[0]=real_c[0].speed;speed[1]=real_c[1].speed;
     
    if(sorted_queue==NULL){all_task_in=1;}
    else{all_task_in=0;}
    if(exe_queue==NULL&&exe_queue2==NULL){available_finished=1;}
    else{available_finished=0;}
    //printf("no problem final\n");
    // fp[0]=NULL;
    //fp[1]=NULL;
    //printf("%d  %d\n",all_task_in,available_finished);
   
    if(output_count==86400){
      //fprintf(fp_out2,"%d   ",sim_time);
      for(i=0;i<2;i++){
	fprintf(fp_out2,"%lf  ",tmpr[i]-273.15);
      }
      
      //fprintf(fp_out2,"\n");
      for(i=0;i<2;i++){
	fprintf(fp_out,"%lf  ",pwr[i]);
      }
      fprintf(fp_out,"\n");
      output_count=0;
      //printf("output ok1\n");
      if(exe_queue!=NULL){
	//printf("output ok in\n");
	//fprintf(fp_out4,"deadline:%f sim_t:%f speed %d,%d,  %d  ",exe_queue->deadline,sim_time,real_c[0].speed,speed_limit[0],exe_queue->number);
	fprintf(fp_out2,"  %d",exe_queue->number);
	//printf("output ok in2\n");
      }else{
	//printf("output ok in3\n");
	fprintf(fp_out2," 0  ");
	//fprintf(fp_out4,"core_0 idle tmpr ");
      }
      //printf("output ok2\n");
      if(exe_queue2!=NULL){
	//fprintf(fp_out4,"%d speed %d ",exe_queue2->number,speed[1]);
	fprintf(fp_out2,"  %d",exe_queue2->number);
      }else{
	//fprintf(fp_out4,"core_1 idle tmpr %f",tmax[1]);
	fprintf(fp_out2,"  0 ");
      }
      //fprintf(fp_out4,"\n");
      fprintf(fp_out2,"\n");
      
    }
    else{output_count++;}
    //printf("output ok\n");
    //post check
    /*
    q_temp=exe_queue;
    //debug_count2=0;

    while(q_temp!=NULL){
      printf("task %d  %d  ",q_temp->number,q_temp->sim_exe);
      //debug_count2+=q_temp->sim_exe;
      q_temp=q_temp->next;
    }
    //printf("\n");
    //fprintf(fp_out4,"left total: %d \n",debug_count2);
    q_temp=exe_queue2;
    while(q_temp!=NULL){
      printf("task[2] %d %d  ",q_temp->number,q_temp->sim_exe);
      q_temp=q_temp->next;
    }
    
    printf("\n");*/
    //printf("utiliaztion: %f   %f\n",real_c[0].utilization,real_c[1].utilization);
    /*
    if(schedule_switch==0){
      if(fabs(rel_bd[0]-rel_bd[1])*2/(2-rel_bd[0]-rel_bd[1])>=FINE){
	schedule_switch=1;
      }
    }
    else{
      if(fabs(rel_bd[0]-rel_bd[1])*2/(2-rel_bd[0]-rel_bd[1])<FINE/5){
	schedule_switch=0;
      }
    }*/
    /*
    if(schedule_switch==0){
      if(fabs(rel_bd[0]-rel_bd[1])>=FINE){
  schedule_switch=1;
      }
    }
    else{
      if(fabs(rel_bd[0]-rel_bd[1])<FINE/5){
  schedule_switch=0;
      }
    }*/
    /*
    printf("swap  ");
    for(i=0;i<5;i++){
      printf("%d  ",swap_counter[i]);
    }*/
    //printf("\n");
    //printf("cycle finished  %d\n",schedule_switch);
    //fprintf(fp_lamda,"core_U[0]  %f   core_U[1]  %f  U1  %f  U2   %f  %d  %Le  %Le\n",core_U[0],core_U[1],real_c[0].utilization,real_c[1].utilization,schedule_switch,rel_bd[0],rel_bd[1]);
    if(rel_bd[0]*rel_bd[1]<=0.9){
      return 0;
    }
    if(real_c[0].utilization>1 || real_c[1].utilization>1){
    	printf("utilization larger than 1\n");
    	return 0;
    }
    
    for(i=0;i<2;i++){
      pwr[i]=0;
    }
    no_of_task=0;
    
  }
  
  
   return 0;
}




/*****************************************************************************************************/
void init(){
   int i;
   for(i=0;i<2;i++){ // block tmpr init
     tmpr[i]=T_INIT+273.15;
     pwr_vec[i]=0;
     pwr[i]=0;
   }
   for(i=2;i<2+EXTRA;i++){ //extra block tmpr init
     tmpr[i]=T_INIT+273.15;
     pwr_vec[i]=0;
     pwr[i]=0;
   }
   for(i=0;i<20;i++){
     temp_lamda[i]=0;
     temp_lamda2[i]=0;
     temp_lamda_pre[i]=T_INIT+273.15;//save previous interval temperature for obd cal
     temp_lamda2_pre[i]=T_INIT+273.15;
     eff_time[i]=0;
     eff_time2[i]=0;
   }
   rel_bd[0]=pow(1,20);rel_bd[1]=pow(1,20); lamda_bd[0]=0;lamda_bd[1]=0;
   for(i=0;i<20;i++){ //block rel init for bd model
     rel_bd_blk[i]=1;
     rel_bd_blk2[i]=1;
     rel_bd_blk_pre[i]=1;
     rel_bd_blk2_pre[i]=1;
   }
   for(i=0;i<2;i++){  // core init
     real_c[i].tmpr=T_INIT;
     real_c[i].speed=0;
     real_c[i].utilization=0;
     real_c[i].number=i+1;
   }
   global_out_count=0;
   for(i=0;i<sizeof(tasklist)/sizeof(struct task);i++){
     tasklist[i].utilization=(double) tasklist[i].sim_exe/(tasklist[i].period*10);
   } 
   
   
     
     
}

int get_lcm(int sub_lcm, int input){
  int f1,f2;
  f1=sub_lcm;f2=input;
  if(sub_lcm==0){return input;}
  else{
    while(f1!=f2){
      if(f1<f2){
	while(f1<f2){
	  f1+=sub_lcm;
	}
      }
      else{
	while(f2<f1){
	  f2+=input;
	}
      }
    }
  }
  return f1;
}

struct task *add_task(int period_lcm){

  int p=0; 
  unsigned int task_num,i,max_cycle,cycle_number,cycle;
  struct task *p1,*p2,*head;
  task_num=sizeof(tasklist)/sizeof(struct task);
  for(i=0;i<task_num;i++){
    cycle=0;
    p=0;
    while(period_lcm-p>0){
      p+=tasklist[i].period;
      cycle++;
    }
    tasklist[i].cycle=cycle;
    //printf("task %d %d cycles\n",i,tasklist[i].cycle);
  }
  max_cycle=tasklist[0].cycle;
  for(i=0;i<task_num;i++){
    if(tasklist[i].cycle>max_cycle)max_cycle=tasklist[i].cycle;
    else continue;
  }
  cycle_number=0;
  while(max_cycle>0){
    for(i=0;i<task_num;i++){
      if(tasklist[i].cycle>0){
	
	p2=(struct task*)malloc(LEN);
	strcpy(p2->name,tasklist[i].name);
	p2->arrival=tasklist[i].arrival+cycle_number*tasklist[i].period;
	p2->deadline=tasklist[i].arrival+tasklist[i].deadline+cycle_number*tasklist[i].period;
	p2->sim_exe=tasklist[i].sim_exe;
	p2->execution=tasklist[i].execution;
	p2->Tss=tasklist[i].Tss;
	p2->cycle=tasklist[i].cycle;
	p2->assign=tasklist[i].assign;
	p2->number=tasklist[i].number;
	p2->period=tasklist[i].period;
	p2->powerfile=tasklist[i].powerfile;
	p2->Tss=tasklist[i].Tss;
	p2->utilization=tasklist[i].utilization;
	p2->next=tasklist[i].next;
	p2->temporary=tasklist[i].temporary;
  p2->avg_pwr=tasklist[i].avg_pwr;
	tasklist[i].cycle-=1;
	if(n==0)head=p2;
	else p1->next=p2;
	p1=p2;
	n+=1;
	      
           
      }
      else continue;
    }
    max_cycle-=1;
    cycle_number+=1;
   }
  return head;
}

struct task *form_list(struct task *sorted,double lcm){
  struct task *head,*p1,*p2,*p3;
  int random,rand_temp;
  double x,v,u;
  head=NULL;p1=NULL;p2=NULL;p3=NULL;
  p1=sorted;
  while(p1!=NULL){
    do{
      u=(double)rand()/RAND_MAX;
      v=(double)rand()/RAND_MAX;
      x=sqrt(-2*log(u))*cos(M_PI*2*v);
      rand_temp=(int)sigma*x+exe_mean[p1->number]+0.5;
    }while(rand_temp<=0 || rand_temp>10);
    random=rand_temp;
    p2=(struct task*)malloc(LEN);
    //printf("start add ok 2\n");
    strcpy(p2->name,p1->name);
    p2->arrival=p1->arrival+lcm;
    p2->deadline=p1->deadline+lcm;
    p2->sim_exe=p1->sim_exe/10*random;
    p2->execution=p1->execution;
    p2->Tss=p1->Tss;
    p2->cycle=p1->cycle;
    p2->assign=p1->assign;
    p2->number=p1->number;
    p2->period=p1->period;
    p2->powerfile=p1->powerfile;
    p2->Tss=p1->Tss;
    p2->utilization=p1->utilization;
    p2->temporary=p1->temporary;
    p2->avg_pwr=p1->avg_pwr;
    p2->next=NULL;
    if(head==NULL){
      head=p2;p3=head;p2=NULL;}
    else{
      p3->next=p2;p3=p3->next;p2=NULL;
    }
    p1=p1->next;
    fprintf(fp_schedule,"%d\n",random);
  }
  return head;
}


struct queue_util *update_util(unsigned int sim_time,struct queue_util *queue,int number){
  struct queue_util *pre,*temp;
 
  if(queue==NULL){
    return queue;
  }
  else{
    temp=queue;
    //printf("ok 1\n");
    while(temp!=NULL){
      // if(temp->next!=NULL){
      if(temp->deadline>sim_time){
	pre=temp;
	temp=temp->next;
      }
      else{
	if(temp->temporary==1){
	  core_U[number]-=temp->utilization;
	}
	if(temp==queue){
	  queue=queue->next;
	  free(temp);
	  temp=queue;
	}
	else{
	  pre->next=temp->next;
	  free(temp);
	  temp=pre->next;
	}
      }
    }
  }
  return queue;
}  

struct queue_util *add_u(unsigned long deadline,double utilization,struct queue_util *queue,int temporary){
  struct queue_util *temp,*new;
  new=(struct queue_util*)malloc(sizeof(struct queue_util));
  new->deadline=deadline;
  new->utilization=utilization;
  new->temporary=temporary;
  new->next=NULL;
if(queue==NULL){
     queue=new;
  }
  else{
    temp=queue;
    while(temp!=NULL){
      if(temp->next==NULL){
	temp->next=new;
	temp=NULL;
      }else{
	temp=temp->next;
      }
    }
  }
  return queue;
}
double cal_u(struct queue_util *queue){
  struct queue_util *temp;
  double u_sum=0;
  temp=queue;
  //printf("ok add_u\n");
  if(temp==NULL){
     //printf("ok add_u1\n");
    u_sum=0;
  }
  else{
    while(temp!=NULL){
      u_sum+=temp->utilization;
      temp=temp->next;
    }
  }
 //printf("ok add_u2\n");
  return u_sum;
}

void print_matrix(double *a,FILE *fp,int m,int n){
  int i,j;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      fprintf(fp,"%f   ",a[i*n+j]);
    }
    fprintf(fp,"\n");
    printf("%d  ",i);
  }
}
