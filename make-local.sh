#!/bin/bash

rm MANIFEST
python setup.py sdist

cd dist
tar -xvzf pslid-0.0.tar.gz
cd pslid-0.0
sudo python setup.py install
cd ../../
