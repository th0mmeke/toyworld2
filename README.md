# ToyWorld2

## Prerequisites

1. Python 2.7.x
1. PyMunk
  - `pip install pymunk`
1. RDKit
1. NetworkX
  - `pip install networkx`
1. If using weighting_functions.replicant_weighting(), python package [python-levenshtein](https://pypi.python.org/pypi/python-Levenshtein/0.12.0) (C-extensions supported with Python 2.2+)
  - `pip install python-levenshtein`
  - If import fails, check permissions on installed package directory (on Linux, `/usr/local/lib/python2.7/dist-packages/Levenshtein`)

## Usage
### Local

`python main.py`

### Using Docker

```
docker run -v <local directory>:/toyworld2/data th0mmeke/toyworld2
```
    
 