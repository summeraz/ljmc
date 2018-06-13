### Installation
To install, clone this repo via
```git clone https://github.com/summeraz/ljmc.git```

### Running MC simulations
Users should only need to interact with the `run.py` file, which
acts as the driver for MC simulations. Here one can specify the
force field, system, and simulation parameters.

To run, simply execute:
```python run.py```

As this is pure Python code, performance is quite poor. Individual
modules have been Cythonized, which has led to modest improvements.
However, this code is far to inefficient to be used for research-grade
simulations, and should be used for demonstration purposes only.
