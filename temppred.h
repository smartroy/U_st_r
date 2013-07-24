#ifndef __TEMPPRED_H_
#define __TEMPPRED_H_
#include "sorttask.h"
float exp_predict(float T_ss, float time_past, float T_init);
void tmpr_table();
void set_ratio(double t0,double t1,double base,int index);
//void readPower1(struct task *exe_queue);
//void readPower2(struct task *exe_queue);
void readPower_core1(struct task *exe_queue);
void readPower_core2(struct task *exe_queue);
void init_standby_power();
struct task *send_to_exe_queue(struct task *exe_queue);
void init_rel_ratio();
void power_input1(struct task *queue ,int speed);
void power_input2(struct task *queue ,int speed);
void power_input1_simple(struct task *queue ,int speed);
void power_input2_simple(struct task *queue ,int speed);
//extern double tmpr_array[4][6][6];
//extern double ratio[2];


#endif
