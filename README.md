# Adiabatic Quantum Computer Simulation

This program simulates the performance of an Adiabatic Quantum Computer (AQC), going from an initial state to the target state, which is the ground state of a system of particles. This simulation does not try to solve the problem of finding this ground state, but to simulate how an AQC should operate in the finding of such state. Thus, the simulation is more didactic than functional.

This should be the first file you read before trying to run the code or understand it. It will give you an idea of how the repository is structured, how to proceed with the installation of the code and understanding what the program does.

## Structure

The repository is divided into the next parts:

- The file Adiabatic_Quantum_Computer_Simulation.pdf contains the explanation about the problem that is trying to be solved by the simulation. It also contains examples of the results that are obtained when executing the code. It is recommended to read it before running the code.

- The file Simulation_parameters.txt contains the input parameters that are used by the code.

- The file Simulation_parameters_explanation.txt contains the explanation of how to interpret Simulation_parameters.txt and how to modify the input parameters.

- The file Simulation.py contains the code, with all the functions necessary to perform the simulation. It reads the input parameters present on Simulation_parameters.txt and generates the initial and target Hamiltonian. Then, it performs the simulation and yields the results.

- The folder Graphs contains the two figures obtained with the default parameters in Simulation_parameters.txt, and will store the graphs generated when executing the program.

- The file test_Simulation.py contains the testing of all the functions in Simulation.py to ensure they work properly.

## Installing and Running

The procedure to understand, install and run the program is the following:

- Download the repository saving all it contains in a certain folder of your computer.

- First, it is suggested to open "Adiabatic_Quantum_Computer_Simulation.pdf" file, to get an understanding of the problem and to know which results you can expect when running the code.

- Once you have read this file, you can proceed with the installations. The code was written on Python 3.7.6, but you can download the latest version of Python.

- Download it and select the installation that includes IDLE, pip and documentation.
**WARNING: It is very important when installing to check the box that says "Add Python 3.x to PATH" to ensure that the interpreter will be placed in your execution path.

- Open the Simulation_parameters.txt. This file contains the initial parameters that the program must read when it is executed.

- To understand how to interpret and modify Simulation_parameters.txt open the "Simulation_parameters_explanation.txt" file.

- Once modified (or not, it depends on what you want) Simulation_parameters.txt file, save the changes.

- Open IDLE (which you have already installed when downloading the last version of python).

- In IDLE click on File and then click on Open.

- Open Simulation.py (which will be on the path where you downloaded the repository).

- The program imports some packages, so, before running it, you must do the following:
  + Open the command line or console. To open it quickly on Windows you can press Windows+R.
  + Write "pip install numpy" and press Enter. A downloading bar appears. Wait until it has been installed.
  + Write "pip install scipy" and press Enter. A downloading bar appears. Wait until it has been installed.
  + Write "pip install matplotlib" and press Enter. A downloading bar appears. Wait until it has been installed.

- If those commands do not initialize any installation it means that during the installation of python the option "Add Python 3.x
to PATH" has not been checked. An option is to reinstall and check it this time. Then try to install the packages again, following
the previous step.

- Once the packages numpy, scipy and matplotlib have been installed, go Simulation.py, which has been opened with IDLE, and select
Run, and then Run Module (Or just press F5).

- The code will be run, and the program will return the time required for the Adiabatic Quantum Computing. Two plots will also be  created, which will be stored in the folder "Graphs", which is found in the same directory as the code. If you don't move the plots to another folder or do not change their name before executing the program again the plots will be replaced by the plots  created in the last execution.

- The code can also be run from the command line:
  + To access the command line press Windows+R (if you have a Windows operating system).
  + In the command line go to the directory in which you have downloaded the archives of the repository. To do it use the command cd.
  Example: If the archives are on a directory with the path C:\Users\User\Documents\Physics\Adiabatic-Quantum-Computer-Simulation    and in the command line you are located in a directory with path C:\Users\User, you must write: cd Documents\Physics\Adiabatic-  Quantum-Computer-Simulation.
  + Once you are located in the directory, type "python Simulation.py" and press Enter.
  + The code will be executed and program will return the things mentioned previously.

## Testing

This repository also contains a testing suite, which ensures that the program is actually doing what it is intended to do. Those tests are written in the file test_Simulation.py. In case you want to run the tests, the procedure is the following:

- Open the command line (Press Windows+R if you have a Windows operating system).

- Write "pip install pytest" and press Enter. A downloading bar appears. Wait until it has been installed.

- If this commands does not initialize any installation it means that during the installation of python the option "Add Python 3.x
to PATH" has not been checked. An option is to reinstall, check it this time and then repeat the previous step (installing   pytest).

- Once pytest is successfully installed go to the directory in which you have downloaded the repository. To do it use the command cd, as explained on the section Installing and Running.

- Once you are located in the directory of the repository on the command line, write on the command line "pytest -v" and press Enter. This will initialize the testing, and if everything works as it should, all test will be shown as "PASSED". For the testing to work it is also important that the file Simulation_parameters.txt is correctly written, or errors could arise. To make sure it is correctly written you can check it and read Simulation_parameters_explanation.txt to have an idea of how the parameters must be modified.
