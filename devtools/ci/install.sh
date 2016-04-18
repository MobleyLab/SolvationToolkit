MINICONDA=Miniconda2-latest-Linux-x86_64.sh
MINICONDA_MD5=$(curl -s http://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget http://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b

export PATH=$HOME/miniconda/bin:$PATH

sudo apt-get update
sudo apt-get install -qq -y g++ gfortran csh g++-multilib gcc-multilib openbabel git

conda update --yes conda
conda config --add channels http://conda.binstar.org/omnia
conda config --add channels https://conda.binstar.org/rdkit

echo "HERE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
conda install --yes mdtraj
echo "HERE2>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
"
conda install --yes jinja2 binstar pip openmoltools packmol pytables 


git clone https://github.com/ParmEd/ParmEd.git
cd ParmEd && python setup.py install
cd ..

conda create -y -n myenv python=$PYTHON_VERSION openmoltools packmol numpy scipy netcdf4 pandas nose openmm pytables mdtraj

source activate myenv
