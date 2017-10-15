#!/bin/bash

rsync -avz --exclude 'results' * hz8228@raijin.nci.org.au:~/advection/
rsync -avz hz8228@raijin.nci.org.au:~/advection/results/* ./results/
