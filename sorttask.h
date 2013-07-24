#ifndef __SORTTASK_H_
#define __SORTTASK_H_
#include <stdlib.h>
#include <stdio.h>
#include "RC.h" 
#include "flp.h"
struct task{
  char name[20];
  int arrival;
  unsigned int deadline;
  int sim_exe;
  int execution;// with Tss below belong to the exp estimation
  double Tss;
  double T_begin;
  int period;
  int cycle;
  int assign;
  int number;
  double utilization;
  //double mean_u;
  int temporary;
  double avg_pwr;
  FILE *powerfile;
  struct task *next;
}task;
struct core{  //core data structure, for future convenience
  struct task *current_task;
  int speed;
  double tmpr;
  double utilization;
  //double mean_u;
  int number;
}core;
double lamda_cal(double temp[40+EXTRA],double speed,flp_t *flp,int core_n,long double lamda);
struct task *sort_task(struct task *initial_queue);
struct task *sort_deadline(struct task *sorted_queue,struct task *exe_queue);
void change_assign(struct task *queue);
int lowest_speed(struct task *queue,int sim_time);
void amda_cal_bd(double temp_lamda[20],double temp_lamda2[20],double temp_lamda_pre[20],double temp_lamda2_pre[20],double interval,double sim_time);
void q_sort(struct task array[],int size);
void quick_sort(struct task array[],int left,int right);
int partition(struct task array[],int left,int right);
void swap(struct task array[],int i,int j);
//qsort for core
void q_sort_c(struct core array[],int size);
void quick_sort_c(struct core array[],int left,int right);
int partition_c(struct core array[],int left,int right);
void swap_c(struct core array[],int i,int j);
#endif
