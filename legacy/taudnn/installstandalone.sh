#! /bin/bash
#needs python 3.6.5
python3 -m pip install --user pip
python3 -m pip install virtualenv
python3 -m venv myenv
. myenv/bin/activate
python -m pip install pip --upgrade
pip install coffea==0.6.47
pip install parsl
pip install h5py
pip install tensorflow==1.14
pip install sklearn
python scripts/analysis_tau.py --chunk 1000 --maxchunk 2 /pnfs/psi.ch/cms/trivcat/store/user/ytakahas/forVini/DNN_v4/2017/0000/flatTuple_124.root --save_h5 --h5folder=h5/ --samples=test
cd scripts/
python train.py --year 17
