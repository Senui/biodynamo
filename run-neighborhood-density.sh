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
ITERATIONS=5

echo "Warming up GPU..."
for it in {1..5}; do
  ./build/neighborhood_density $POPULATION $ITERATIONS 400;
done

for MB in 900 600 500 400 350 300 275 250; do
  echo "Running benchmark with max-bound = $MB"
  for it in {1..5}; do
    ./build/neighborhood_density $POPULATION $ITERATIONS $MB >> gpu.csv
    echo -n "," >> gpu.csv
  done
  echo "" >> gpu.csv
done

