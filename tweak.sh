# tweak the matlab python interface installer into accepting python 3.8
MATLABROOT="/home/anton/pfiles/MATLAB"
# run matlab -nodesktop -nosplash -nojvm -batch matlabroot to obtain matlab root directory.
INSTALLER_PATH="$MATLABROOT/extern/engines/python"

cd $INSTALLER_PATH
sed -i.bak s/"'2.7', '3.6', '3.7'"/"'2.7', '3.6', '3.7', '3.8'"/g setup.py
python setup.py build

cd build/lib/matlab/engine
sed -i.bak s/"'2_7', '3_6', '3_7'"/"'2_7', '3_6', '3_7', '3_8'"/g __init__.py
sed -i s/"_PYTHONVERSION = _version"/"_PYTHONVERSION = '3_7'"/g __init__.py

cd $INSTALLER_PATH
cd dist/matlab/engine
sed -i.bak s/"'2_7', '3_6', '3_7'"/"'2_7', '3_6', '3_7', '3_8'"/g __init__.py
sed -i s/"_PYTHONVERSION = _version"/"_PYTHONVERSION = '3_7'"/g __init__.py
