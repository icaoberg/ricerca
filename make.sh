#!/bin/bash

rm MANIFEST
python setup.py sdist
#scp -r dist/coconut*tar.gz bhcho@omepslid.compbio.cs.cmu.edu:~/
#scp -r dist/coconut*tar.gz icaoberg@lanec1.compbio.cs.cmu.edu:~/
