#!/bin/bash

# With a total population of 2,000,000 agents:
# max-bound = 2000 -> ~0 neighbor per agent
# max-bound = 900 -> ~1 neighbor per agent
# max-bound = 600 -> ~3.5 neighbors per agent
# max-bound = 500 -> ~6 neighbors per agent
# max-bound = 400 -> ~11.5 neighbors per agent
# max-bound = 350 -> ~17 neighbors per agent
# max-bound = 300 -> ~27 neighbors per agent
# max-bound = 275 -> ~35 neighbors per agent
# max-bound = 250 -> ~47 neighbors per agent

POPULATION=2000000
ITERATIONS=20

echo "Warming up GPU..."
for it in {1..5}; do
  ./build/neighborhood_density $POPULATION $ITERATIONS 400;
done

# Run on one NUMA domain (0-15, 48-63 on olgpu-01)
for MB in 900 600 500 400 350 300 275 250; do
  FILENAME=gpu.csv
  echo "Running benchmark with max-bound = $MB"
  for it in {1..5}; do
    OMP_NUM_THREADS=32 taskset --cpu-list 0-15,48-63 ./build/neighborhood_density $POPULATION $ITERATIONS $MB >> $FILENAME
    echo -n "," >> $FILENAME
  done
  echo "" >> $FILENAME
done
