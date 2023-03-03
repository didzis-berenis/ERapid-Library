# ERapid-Library
Libraries for coupling [Elmer FEM](https://www.csc.fi/web/elmer) and [RapidCFD](https://github.com/Atizar/RapidCFD-dev) + test cases. This library is based on [EOF-Library](https://github.com/jvencels/EOF-Library) but adapted for use with [RapidCFD](https://github.com/Atizar/RapidCFD-dev).

## Requirements ##
* [Ubuntu 20.04](https://releases.ubuntu.com/focal/) or [Ubuntu 18.04](https://releases.ubuntu.com/bionic/)
* [EOF-Library](https://github.com/jvencels/EOF-Library) with [OpenFOAM-7](https://openfoam.org/release/7/) and [Elmer FEM](https://www.csc.fi/web/elmer) (Recommended)

## Installation steps ##
Make sure that each step is executed correctly before proceding to the next one.

#### 1. Install required packages #### 
1.1 
```
sudo apt-get update
sudo apt-get install git gmsh cmake libblas-dev liblapack-dev build-essential curl flex gfortran libboost-system-dev libboost-thread-dev -y
sudo apt-get install bison automake zlib1g-dev libopenmpi-dev openmpi-bin libreadline-dev libncurses-dev libxt-dev freeglut3-dev -y
```

1.2a **Ubuntu 18.04**
```
sudo apt-get install qt4-dev-tools libqt4-dev libqt4-opengl-dev libqtwebkit-dev -y
```

1.2b **Ubuntu 20.04**
```
sudo apt-get install qt5-default gcc-7 g++-7 -y
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 800 --slave /usr/bin/g++ g++ /usr/bin/g++-9
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 900 --slave /usr/bin/g++ g++ /usr/bin/g++-7
```

#### 2. Install Cuda Toolkit ####
**CUDA Toolkit 10.1** tested on **Ubuntu 18.04**
**CUDA Toolkit 11.1** tested on **Ubuntu 20.04**
See [CUDA Toolkit Archive](https://developer.nvidia.com/cuda-toolkit-archive) for documentation

* **CUDA Toolkit 10.1 on Ubuntu 18.04**
```
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-ubuntu1804.pin
sudo mv cuda-ubuntu1804.pin /etc/apt/preferences.d/cuda-repository-pin-600
sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
sudo add-apt-repository "deb http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/ /"
sudo apt-get update
sudo apt-get install -y cuda-toolkit-10-1
export PATH=$PATH:/usr/local/cuda-10.1/bin
export CUDADIR=/usr/local/cuda-10.1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-10.1/lib64
```

* **CUDA Toolkit 11.1 on Ubuntu 20.04**
```
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/7fa2af80.pub
sudo add-apt-repository "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/ /"
sudo apt-get update
sudo apt-get -y install cuda-toolkit-11-1
sudo apt-get update
export PATH=$PATH:/usr/local/cuda-11.1/bin
export CUDADIR=/usr/local/cuda-11.1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-11.1/lib64
```

#### 3. Install ThirdParty and RapidCFD ####
3.1
```
cd $HOME
mkdir RapidCFD
cd RapidCFD
git clone https://github.com/Atizar/RapidCFD-dev
wget -O - http://dl.openfoam.org/third-party/7 | tar xvz
mv ThirdParty-7-version-7 ThirdParty-dev
export WM_PROJECT=RapidCFD
export WM_PROJECT_VERSION=dev
export FOAM_INST_DIR=$HOME/$WM_PROJECT
```

3.2a **With CUDA 10.0**
Openmpi version should match with the OS (Ubuntu 18.04) by default

3.2b **With CUDA 11.1**
Make sure to patch file
```
$HOME/RapidCFD/RapidCFD-dev/etc/config/settings.sh
```
On line #403 change
``    export FOAM_MPI=openmpi-4.0.2``
to
``    export FOAM_MPI=openmpi-4.0.5``

3.3 Update compilation rules according to your [GPU Compute Capability](https://developer.nvidia.com/cuda-gpus)
Patch file
```
$HOME/RapidCFD/RapidCFD-dev/wmake/rules/linux64Nvcc/c
```
according to your compute capability. For example if your compute capability is 7.5, then on line #5
```    cc          = nvcc -m64 -arch=sm_30 
```
change ``-arch`` value
```    cc          = nvcc -m64 -arch=sm_75 
```

Additionally patch 
```
$HOME/RapidCFD/RapidCFD-dev/wmake/rules/linux64Nvcc/c++
```
and change line #10
```    CC          = nvcc -Xptxas -dlcm=cg -std=c++11 -m64 -arch=sm_30
```

3.4a **With CUDA 10.0**
change ``-arch`` value
```    CC          = nvcc -Xptxas -dlcm=cg -std=c++11 -m64 -arch=sm_75
```

3.4b **With CUDA 11.1**
change ``-arch`` and ``-std`` value
```    CC          = nvcc -Xptxas -dlcm=cg -std=c++14 -m64 -arch=sm_75
```
3.5 Source project bashrc
```
source $FOAM_INST_DIR/$WM_PROJECT-$WM_PROJECT_VERSION/etc/bashrc
```

3.6 Update ThirdParty compilation configuration by applying patch to file
```
$HOME/RapidCFD/ThirdParty-dev/Allwmake
```
On line #85 change
``    configOpt="--with-sge"``
to
``    configOpt="--with-cuda"``
and on line #100 change
``    --enable-mpi-fortran=none \``
to
``    --enable-mpi-fortran \``
3.7 Navigate to ThirdParty installation directory
```
cd ThirdParty-dev
```

3.8a **With CUDA 10.1**
```
wget -O - https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.2.tar.gz | tar xvz
```
3.8b **With CUDA 11.1**
```
wget -O - https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.5.tar.gz | tar xvz
```
3.9 Installation of ThirdParty
``./Allwmake``

3.10 If you plan to use ERapid-Library with only one GPU, then make sure to apply patch to file
```
$HOME/RapidCFD/RapidCFD-dev/src/Pstream/mpi/UPstream.C
```
On line #87 change 
``    if (numprocs <= 1)``
to
``    if (numprocs < 1)``

3.11 Installation of RapidCFD
```
cd ../RapidCFD-dev
export WM_NCOMPPROCS=16
./Allwmake
```

#### 4. Install Elmer ####
4.1
```
cd ../..
mkdir elmerRapid
cd elmerRapid
git clone https://github.com/ElmerCSC/elmerfem
mkdir build
cd build
```
4.2a **With CUDA 10.1**
```
cmake -DWITH_MPI=TRUE --with-mpi-inc-dir=$HOME/RapidCFD/ThirdParty-dev/openmpi-4.0.2/ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../elmerfem
```
4.2b **With CUDA 11.1**

```
cmake -DWITH_MPI=TRUE --with-mpi-inc-dir=$HOME/RapidCFD/ThirdParty-dev/openmpi-4.0.5/ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../elmerfem
```
4.3 Installation of Elmer
``make -j install``

#### 5. Download and install ERapid-Library ####

```
cd ../..
git clone https://github.com/didzis-berenis/ERapid-Library
. ERapid-Library/etc/bashrc
export ELMER_HOME=$HOME/elmerRapid/install
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ELMER_HOME/lib
export PATH=$PATH:$ELMER_HOME/bin
eofCompile
```

#### 6. Configure environment variables ####
* Configure tools for switching between OpenFOAM-7 and RapidCFD
```
cp $HOME/ERapid-Library/setupTools/PROJECTDIR $HOME/
cp $HOME/ERapid-Library/setupTools/sourceOF $HOME/
cp $HOME/ERapid-Library/setupTools/sourceRapid $HOME/
'export OLD_PATH=$PATH' >> $HOME/.bashrc
'export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH' >> $HOME/.bashrc
'export WM_PROJECT=OpenFOAM' >> $HOME/.bashrc
'export WM_PROJECT_VERSION=7' >> $HOME/.bashrc
'export FOAM_INST_DIR=$HOME/$WM_PROJECT' >> $HOME/.bashrc
'export ELMER_HOME=$HOME/elmer/install' >> $HOME/.bashrc
'export EOF_HOME=$HOME/EOF-Library' >> $HOME/.bashrc
'source $HOME/PROJECTDIR' >> $HOME/.bashrc
```

* (Optional) Decrease solver log amount.
In file ``$HOME/RapidCFD/RapidCFD-dev/etc/controlDict``
change value of ``GAMGAgglomeration`` , ``SolverPerformance`` and ``lduMatrix``  to 0

## Work with ERapid-Library ##

* Before first use make sure to compile solvers
```
source $HOME/sourceRapid
cd ERapid-Library/solvers
./Allwmake
```

* Create an empty folder
```
cd ../
mkdir runs
cd runs
```

* Copy test simulation
``` 
cp -r tests/mhdFurnace runs
```
* Prepare case
```
cd runs/mhdFurnace
. sourceOF
gmsh -3 furnace.geo -o furnace.msh -format msh2
gmshToFoam furnace.msh
changeDictionary
decomposePar
```

* Extract Elmer mesh (or make your own with [Salome](https://www.salome-platform.org/) and [salomeToElmer](https://github.com/jvencels/salomeToElmer) script)
``tar -xf meshElmer.tar.xz``

* (Recommended) Create O2E pair files with EOF-Library before running the case with ERapid-Library

* Run rapidCFD simulation on 2 GPUs and Elmer simulation on 6 CPU cores/threads:
```
source $HOME/sourceRapid
ElmerGrid 2 2 meshElmer -metis 6
mpirun -np 2 MHDPimpleAverage -parallel -devices "(0 1)" : -np 6 ElmerSolver_mpi case.sif > log.txt &
```
* Alternatively run rapidCFD simulation on 1 GPU and Elmer simulation on 4 CPU cores/threads:
```
source $HOME/sourceRapid
ElmerGrid 2 2 meshElmer -metis 4
mpirun -oversubscribe -np 1 MHDPimpleAverage -parallel : -np 4 ElmerSolver_mpi case.sif > log.txt &
```
* Simulation results will appear in host system *runs* folder

