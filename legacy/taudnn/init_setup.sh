#!/usr/bin/env bash

<<<<<<< Updated upstream
python3.6 -m pip install --user pip
python3.6 -m pip install virtualenv
python3.6 -m virtualenv --python=/bin/python3.6 myenv
. myenv/bin/activate
python -m pip install pip --upgrade
pip install coffea==0.6.47
pip install parsl
pip install h5py
pip install tensorflow==1.14
pip install sklearn
=======
pyenv virtualenv COFFEA
source activate COFFEA
python -m pip install --user pip --upgrade
>>>>>>> Stashed changes

if [ ! -d "$HOME/bin" ]
then
    echo "Creating bin folder under the home area."
    mkdir $HOME/bin
fi

export PARSLINSTALL=$(python -c "import parsl; print(parsl.__path__[0])")
ln -s $PARSLINSTALL/executors/high_throughput/process_worker_pool.py $HOME/bin
chmod +x $HOME/bin/process_worker_pool.py
