# Project Dependencies

This project depends on FastJet and its extensions fjcontrib and IFNPlugin, as well as boost Library. Below are the steps to install and configure these dependencies:

## FastJet

1. Download and install FastJet:
   ```bash
   wget http://fastjet.fr/repo/fastjet-3.4.2.tar.gz
   tar xzf fastjet-3.4.2.tar.gz
   cd fastjet-3.4.2
   ./configure --prefix=$HOME/.local/
   make
   make check
   make install
   cd ..
   ```
   This will install FastJet in the `.local` folder under your home directory.

## fjcontrib

1. Download and install fjcontrib (make sure to use a version compatible with your FastJet):
   ```bash
   wget https://fastjet.hepforge.org/contrib/downloads/fjcontrib-1.054.tar.gz
   tar xzf fjcontrib-1.054.tar.gz
   cd fjcontrib-1.054
   ./configure --fastjet-config=$HOME/.local/bin/fastjet-config
   make
   make check
   make install
   cd ..
   ```

## IFNPlugin

IFNPlugin is a FastJet plugin for implementing interleaved flavour neutralisation. Clone it from the following location:
```bash
git clone https://github.com/jetflav/IFNPlugin.git
cd IFNPlugin
make
make install
cd ..
```

## GHSAlgo

GHSAlgo is another FastJet plugin for implementing interleaved flavour neutralisation. Clone it from the following location:
```bash
git clone https://github.com/jetflav/GHSAlgo.git
cd GHSAlgo
make
make install
cd ..
```

## Boost Library

The project also requires support from the Boost library:

```bash
wget https://boostorg.jfrog.io/artifactory/main/release/1.85.0/source/boost_1_85_0_rc3.tar.bz2
tar xjf boost_1_85_0_rc3.tar.bz2
#cp /home/xiaoxue/code/2.SpinCorrelation/history/spin/ways/Machine_Learning/boost_1_85_0_rc3.tar.bz2 .
cd boost_1_85_0
./bootstrap.sh --prefix=$HOME/.local/
./b2 install
cd ..
```

## Adding the `bin` Directory to Your PATH

To add the `bin` directory of the current project to your `PATH` environment variable, run the following command in your terminal:

```sh
echo "export PATH=\"\$PATH:$(pwd)/bin\"" >> ~/.zshrc
source ~/.zshrc
```
Edit the `bin/compile` script and change line 15 to your own path:
Adding the `IFN PATH` to Your bin/compile
```sh
IFNPLUGIN_PATH=/path/to/IFNPlugin
For example:
IFNPLUGIN_PATH=/home/zlin/app/IFNPlugin
```
Adding the `GHS PATH` to Your bin/compile
```sh
IFNPLUGIN_PATH=/path/to/GHSAlgo
For example:
IFNPLUGIN_PATH=/home/zlin/app/GHSAlgo
```

## Creating the Directory
Make sure you have a large folder in your computer, then typing:
```sh
ln -s [path] data_storage
```