cc=gcc
test:task_test.o temperature.o flp.o util.o RCutil.o temppred.o sort_task.o tilts.o
	$(cc) -o test -D_FILE_OFFSET_BITS=64 -lm task_test.o temperature.o RCutil.o flp.o util.o temppred.o sort_task.o tilts.o -lm
task_test.o:task_test.c   
	$(cc) -c task_test.c -lm
temperature.o:temperature.c RC.h flp.h util.h 
	$(cc)  -DDEBUG -O0 -g -Wall -c temperature.c -lm
RCutil.o:RCutil.c RC.h flp.h util.h 
	$(cc) -DDEBUG -O0 -g -Wall -c RCutil.c -lm
flp.o:flp.c flp.h util.h 
	$(cc)  -DDEBUG -O0 -g -Wall -c flp.c -lm
util.o:util.c util.h
	$(cc)  -DDEBUG -O0 -g -Wall -c util.c -lm
temppred.o:temppred.c temppred.h 
	$(cc)  -DDEBUG -O0 -g -Wall -c temppred.c -lm
sort_task.o:sort_task.c sorttask.h
	$(cc)  -DDEBUG -O0 -g -Wall -c sort_task.c -lm
tilts.o:tilts.c tilts.h
	$(cc)  -DDEBUG -O0 -g -Wall -c tilts.c -lm



