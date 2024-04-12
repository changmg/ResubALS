# Efficient Resubstitution-Based Approximate Logic Synthesis

This is an approximate logic synthesis flow that automatically generates approximate circuits under a user-specified error constraint.

<!-- Its overall flow is shown below: -->

<!-- <img src="fig/flow.jpg" alt="flow" style="zoom: 100%;" /> -->

<!-- For more details, you can refer to the following paper: -->

<!-- [Chang Meng, Hanyu Wang, Yuqi Mai, Weikang Qian, and Giovanni De Micheli, "*VACSEM: Verifying average errors in approximate circuits using simulation-enhanced model counting*," to appear in Design, Automation & Test in Europe Conference & Exhibition (DATE), Valencia, Spain, 2024.](https://github.com/changmg/changmg.github.io/blob/master/preprint/Meng_Wang_Mai_Qian_DeMicheli_VACSEM-Verifying_Average_Errors_in_Approximate_Circuits_Using_Simulation-Enhanced_Model_Counting.pdf) -->

<!-- ## Project Structure

There are two important folders in the root directory of VACSEM project:

1. Circuit2Cnf: This folder contains the first sub-project of VACSEM, which corresponds to Phase 1 in the above figure, i.e., circuit-aware construction of #SAT problems

2. SimSharpSat: This folder contains the second sub-project of VACSEM, which corresponds to Phase 2 in the above figure, i.e., simulation-enhanced model counting.

   ***The sub-project "SimSharpSat" is developed based on the open-source #SAT tool GANAK. We are deeply grateful to authors of GANAK for their clear and easy-to-use code!***
   
   ***Note that in our paper, we are comparing to the [GANAK](https://github.com/meelgroup/ganak/tree/9b20008ea07c6be9f7be6ca40384ff01cc026413) with the git commit ID "9b20008ea07c6be9f7be6ca40384ff01cc026413", released on July 9, 2023, which works on a single CPU thread.*** -->

## Dependencies 

- Reference environment, **Ubuntu 20.04 LTS** with the following tools and libraries:

  - Tools: gcc 10.3.0 & g++ 10.3.0 & [cmake](https://cmake.org/) 3.16.3

    You can install these tools with the following command:

    ```shell
    sudo apt install gcc-10
    sudo apt install g++-10
    sudo apt install cmake
    ```

    You also need to check whether the default versions of gcc and g++ are 10.3.0:

    ```shell
    gcc --version
    g++ --version
    ```

    If the default versions of gcc and g++ are not 10.3.0, please change them to 10.3.0.

  - Libraries: [libboost](https://www.boost.org/) 1.74.0, libreadline 8.0-4

    You can install these libraries with the following command:

    ```shell
    sudo apt install libboost1.74-all-dev
    sudo apt install libreadline-dev
    ```

  - Other tool: yosys

    ```shell
    sudo apt install yosys
    ```

- **Alternatively, we prepare a docker image containing the dependencies:**

  https://hub.docker.com/r/changmeng/als_min

## Download

This project contains a submodule: open-source logic synthesis and verification tool *abc*

To download the project:

```shell
git clone --recursive https://github.com/changmg/VACSEM.git
```

**Please do add the parameter "--recursive", or the submodule will not be downloaded**

## Build

- To build, go to the root directory of the project, and then execute:

```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..
```

If you compile successfully, you will obtain the following executable program:

resubals.out

- To clean up, go to the root directory of the project, and then execute:

```
rm -rf build
```

## Run

To seek help, you can run:

```shell
./resubals.out -h
```

You will get the following illustration:
```
usage: ./resubals.out --accCirc=string [options] ... 
options:
      --accCirc              path to accurate circuit (string)
      --standCell            path to standard cell library (string [=./input/standard-cell/nangate_45nm_typ.lib])
      --outpPath             path to approximate circuits (string [=tmp])
      --metrType             error metric type: ER, MED, NMED, MSE, MHD, NMHD (string [=NMED])
      --distrType            error distribution type: UNIF, ENUM (string [=UNIF])
      --seed                 seed for randomness (unsigned int [=0])
      --errUppBound          error upper bound (double [=0.15])
      --nFrame               #Monte Carlo samples, nFrame should be an integer multiple of 64 (int [=102400])
      --nFrame4ResubGen      #patterns for AppResub Generation (int [=64])
      --maxCandResub         max #candidate AppResubs (int [=100000])
      --nThread              number of threads (int [=4])
      --enableFastErrEst     when this option is enabled, the program performs faster approximate error estimation;
                                otherwise, the program performs slower accurate error estimation
      --enableMeasureMode    when this option is enabled, the program measures the quality of the circuit specified by `appCirc' option;
                                otherwise, the program performs approximate logic synthesis
      --appCirc              path to approximate circuit,
                                this option is only used in the case when the `enableMeasureMode' option is active  (string [=])
  -h, --help                 print this message
```

**Example command 1**

```shell
./resubals.out --accCirc ./input/benchmark/iscas/c880.blif --standCell ./input/standard-cell/nangate_45nm_typ.lib --outpPath ./tmp/ --metrType ER --errUppBound 0.05 --nThread 4
```

In this example, 

- The original circuit is "./input/benchmark/iscas/c880.blif"
- The stand cell library is "./input/standard-cell/nangate_45nm_typ.lib"
- The approximate circuits will be outputed to "./tmp/"
- The error metric is error rate (ER)
- The error upper bound is 0.05
- The applied CPU threads are 4