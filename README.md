# Simple Variational Monte Carlo solve for FYS4411

Example class structure for the first VMC project of FYS4411 (spring 2016). You may, if you wish, fork this repository and make it the basis of your project. If you choose to do this, you will have to implement a lot of the crucial functions yourself. The relevant functions you need to implement are spread throughout the project, but they are all commented with a note saying what each function should do.

Please note that this is only a start, and when you have implemented all of these functions you will only have completed the first exercise. However, once this is done, you will have a very good basis for further work, and adding functionality will be easier with a good class structure.

If you want to write your own code from scratch, you are of course welcome to do so, and feel free to use this code as inspiration for your own class structure.

- If you choose to use this code as a basis for your work, the first thing you should do is fork it, pull it down to your computer, and make sure it compiles and runs. See the next section on how to compile and run the project. After this you should spend at least 10 minutes looking at the structure and familiarizing yourself with how the classes interact with eachother. 
- A good way to do this may be to simply start at the top of the main.cpp file, and go through all the calls to the System class functions. Consider also the base classes WaveFunction, Hamiltonian, and InitialState and see which functions are virtual (which functions NEED to be implemented by any sub-class).
- You can skip over the output function in the Sampler class and the entire Random class.


## Compilling and running the project
There are now several options you can use for compiling the project. If you use QT Creator, you can import this project into the IDE and point it to the `.pro`-file. If not, you can use CMake to create a Makefile for you which you can then run. You can install CMake through one of the Linux package managers, e.g., `apt install cmake`, `pacman -S cmake`, etc. For Mac you can install using `brew install cmake`. Other ways of installing are shown here: [https://cmake.org/install/](https://cmake.org/install/).

### Compiling the project using CMake
In a Linux/Mac terminal this can be done by the following commands
```bash
# Create build-directory
mkdir build

# Move into the build-directory
cd build

# Run CMake to create a Makefile
cmake ../

# Make the Makefile using two threads
make -j2

# Move the executable to the top-directory
mv vmc ..
```
Or, simply run the script `compile_project` via
```bash
./compile_project
```
and the same set of commands are done for you. Now the project can be run by executing
```bash
./vmc
```
in the top-directory.

#### Cleaning the directory
Run `make clean` in the top-directory to remove the executable `vmc` and the `build`-directory.

#### Windows
Compilation of the project using Windows is still an open question to me, but please include a pull-request if you've got an example. CMake should be OS-independent, but `make` does not work on Windows.

## Completing the missing parts ##
Here follows a suggestion for how you can work to complete the missing parts of the code:
- Start by implementing the Gaussian wave function: Write the evaluate function. Assume for now that the number of particles is always one, and the number of dimensions is always 1. Next, compute the Laplacian analytically, and implement the double derivative function.
- Secondly, use the Random class (or your own favorite random number generator, should you have one) to implement the missing part of the setupInitialState part of the RandomUniform class. Note that this should be pretty straight forward and simple.
- Next, implement the metropolisStep function in the System class. Implement also the small missing part of the runMetropolisSteps function.
- Now, the last big thing needed is to implement the energy calculation. This is done by the Hamiltonian sub-class HarmonicOscillator. Here you will have to use the Laplacian you calculated for the wave function earlier.
- Now the code should be functioning and you should see (somewhat) reasonable results. Try to set the oscillator frequency to 1 and calculate analytically the energy of the oscillator. Recall the form of the ground state wave function of the harmonic oscillator, and set the parameter alpha accordingly. What is the resulting energy?
- If this energy is NOT correct, the last bit missing is to take a look at the computeAverages function in the Sampler class. What is missing here?

