# VolCE

## What is it?
VolCE is a tool designed for computing or estimating the size of the solution space of an SMT(LA) formula. If you are not familiar with VolCE, you can read the [manual](manual.pdf). It is licensed under the [GNU General Public License](COPYING).

## Directories
| Name           | Description   |
|  ------------- | ------------- |
| [release_64bit/](release_64bit/)	| Binary files compiled on Ubuntu 14.04 64-bit, including Z3, Vinci and LattE |
| [src/](src/) | Sorce codes |
| [benchmarks.zip](examples.zip) | The complete set of benchmarks |
| [build.sh](build.sh) | Shell for building VolCE |
| [COPYING](COPYING) | GNU GPL lincese |
| [makefile](makefile) | makefile |
| [manual.pdf](manual.pdf) | Manual for VolCE |
| [vinci-1.0.5.zip](vinci-1.0.5.zip) | Source codes for a modified version of Vinci |
| [z3-master.zip](z3-master.zip) | Source codes of Z3 |

## Build status
This release of VolCE has been successfully built on the following operating systems:
* Ubuntu 18.04 on 64-bit with g++ 7.3.0
* Ubuntu 14.04 on 64-bit with g++ 4.8.4
* Ubuntu 12.04 on 32-bit with g++ 4.8.1

## Building VolCE
* Step 1: Make sure that g++ (version 4.8 or higher version) is installed on your machine (you can type "g++ -v" to check this).
* Step 2: The functionality of VolCE is dependent on some other libraries: [boost](http://www.boost.org/), [glpk](http://www.gnu.org/software/glpk/), and [Armadillo](http://arma.sourceforge.net/).
* Step 3: Execute:
```bash
sh build.sh
```
* Step 4: Build and install [LattE](https://www.math.ucdavis.edu/~latte/). Then move the executable files (*count* and *scdd\_gmp*) into directory *bin/*. For 64-bit user, one could copy *count* and *scdd\_gmp* directly from [compiled binary files](release_64bit/volce3_release_64bit.zip).

Quick test, simply execute:
```bash
./volce3 -h
```
VolCE should pop up the help menu by this command.

**Note:** Move or copy *volce3* with directory *bin/* together since VolCE3 requires the tools in *bin/*.

### Quick guide for building on Ubuntu

Execute:

```bash
sudo apt-get install g++
sudo apt-get install libglpk-dev
sudo apt-get install libboost-dev
sudo apt-get install libarmadillo-dev
sh build.sh
```

Build and install [LattE](https://www.math.ucdavis.edu/~latte/), then move the executable files (*count* and
*scdd\_gmp*) into directory *bin/*. For 64-bit user, one could copy *count* and *scdd\_gmp* directly from [compiled binary files](release_64bit/volce3_release_64bit.zip).

**Note:** On older versions of Ubuntu, you may need install g++-4.8 (or higher version) by hand.

### Questions/Feedback/Comments ###
Please contact:

  1. Cunjing Ge ([gecj@ios.ac.cn](mailto:gecj@ios.ac.cn))


Enjoy!



