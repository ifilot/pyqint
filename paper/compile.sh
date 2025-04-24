#!/bin/bash

# Perform Docker compilation under Windows using Git BASH

# provide the path to the location where the paper is stored; note that under
# Windows Git Bash, two slashes need to be used instead of one
PAPERPATH="$PWD"

docker run --rm \
--volume $PAPERPATH:/data \
--user $(id -u):$(id -g) \
<<<<<<< HEAD
--env JOURNAL=joss \
=======
--env JOURNAL=jose \
>>>>>>> 7c9267227a66a7ad95712d3a19d40134857f93dd
openjournals/inara
