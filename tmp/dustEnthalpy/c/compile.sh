#/bin/bash

CFLAGS="-O3 -lm -pg"

gcc -c string_utils.c -o string_utils.o $CFLAGS
gcc -c run_utils.c -o run_utils.o $CFLAGS
gcc -c Qabs_data.c -o Qabs_data.o $CFLAGS
gcc -c dustRadiation.c -o dustRadiation.o $CFLAGS 
gcc -c TemperatureDistribution.c -o TemperatureDistribution.o $CFLAGS
gcc -c test_main.c -o test_main.o $CFLAGS 
gcc -o tester test_main.o Qabs_data.o string_utils.o run_utils.o dustRadiation.o TemperatureDistribution.o $CFLAGS 
