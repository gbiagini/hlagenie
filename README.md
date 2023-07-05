# HLAgenie
Python package for dealing with HLA sequence data

## Quick-Start Guide

## Installation

### Install from PyPI

``` pip install HLAgenie ```

### Install from source

``` 
python3 -m venv venv
source venv/bin/activate
python setup.py install
```

## Available loci:

1. A
2. B
3. C
4. DPA1
5. DPB1
6. DQA1
7. DQB1
8. DRB1

## Supported functions

### Initialize the HLAGenie object

Instantiate the object with the IPD-IMGT/HLA database version. This will create a SQLite database of sequences in the /tmp directory.

```python
import hlagenie

genie = hlagenie.init("3510")
```

### Retrieve amino acid position from mature protein sequence

Get the amino acid residue at a given position in a given allele's mature protein sequence. This is 1-indexed, not 0-indexed. 

```python
genie.getAAposition("A*01:01",1)
```

## TODO

- [ ] Add documentation for other functions
- [ ] Add example output