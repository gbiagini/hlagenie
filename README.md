# HLAGenie
Python package for dealing with HLA sequence data

![hlagenie_logo.png](images/hlagenie_logo.png)

Heavily inspired by and dependent on the [pyARD](https://www.github.com/nmdp-bioinformatics/py-ard) package. The `HLAGenie` package strives to streamline and standardize the bulk handling of HLA sequence data and provide functions to simplify matching analysis.

## Table of Contents
- [Installation](#installation)
  - [Install from PyPI](#install-from-pypi)
  - [Install from source](#install-from-source)
- [Using HLAGenie](#using-hlagenie)
  - [Using `hlagenie` from Python](#using-hlagenie-from-python)
    - [Initialize `hlagenie`](#initialize-hlagenie)
    - [Accessing amino acid sequence dictionaries for HLA alleles](#accessing-amino-acid-sequence-dictionaries-for-hla-alleles)
    - [Retrieve amino acid position from mature protein sequence](#retrieve-amino-acid-position-from-mature-protein-sequence)
    - [Retrieve amino acid substring from mature protein sequence](#retrieve-amino-acid-substring-from-mature-protein-sequence)
    - [Retrieve epitope from mature protein sequence](#retrieve-epitope-from-mature-protein-sequence)
    - [Check if two alleles have a mismatch at a given position](#check-if-two-alleles-have-a-mismatch-at-a-given-position)
    - [Count the number of amino acid mismatches at a position between donor and recipient](#count-the-number-of-amino-acid-mismatches-at-a-position-between-donor-and-recipient)
    - [Count the number of amino acid mismatches between donor and recipient at a given position given the alleles](#count-the-number-of-amino-acid-mismatches-between-donor-and-recipient-at-a-given-position-given-the-alleles)

## Installation

### Install from PyPI

``` pip install hlagenie ```

### Install from source

``` 
python3 -m venv venv
source venv/bin/activate
python setup.py install
```

## Using HLAGenie

### Using `hlagenie` from Python

`hlagenie` is intended to simplify the handling of HLA sequence data.

#### Initialize `hlagenie`

Import `hlagenie` package

```python
import hlagenie
```

Initialize `GENIE` object with a version of the IMGT/HLA database.

```python
import hlagenie

genie = hlagenie.init("3510")
```

The default behavior is to use ungapped sequences for HLA alleles. If you would prefer the alleles to maintain the gaps that exist in the sequence alignment, pass `ungap = False` to the `init` function.

```python
import hlagenie

genie = hlagenie.init("3510", ungap = False)
```

The first time an object is instantiated with a given IMGT/HLA database version, the package will download the appropriate MSF files from the IMGT/HLA GitHub repository and create a SQLite database in the `/tmp` folder.

#### Accessing amino acid sequence dictionaries for HLA alleles

The `GENIE` object contains dictionaries of amino acid sequences for each HLA allele. The keys for the dictionaries are the HLA allele names. The values are the amino acid sequences for the full or mature protein sequence dictionary.

All of the keys are two-field alleles given the shared protein sequence of these alleles.

```python
genie.seqs # mature protein sequences
genie.full_seqs # full protein sequences
```

#### Retrieve amino acid position from mature protein sequence

To get a given amino acid position from a given HLA allele, you can use the `getAAposition` function. This function is 1-indexed to match standard IMGT/HLA database nomenclature.

For this and following functions, if a three- or four-field allele name is passed, a call is made to `py-ard` to reduce to the two-field level.

```python
genie.getAAposition("A*01:01",1) # returns "G"
```

#### Retrieve amino acid substring from mature protein sequence

To get a given amino acid substring from a given HLA allele, you can use the `getAAsubstring` function. This function is 1-indexed as well and is inclusive of the start and end positions.

```python
genie.getAAsubstring("A*01:01",1,10) # returns "GSHSMRYFFT"
```

#### Retrieve epitope from mature protein sequence

If you pass an allele name and a list of positions to the `getEpitope` function, `hlagenie` will return a formatted epitope string.

```python
genie.getEpitope("A*01:01",[1,2,3,4,5]) # returns "1G_2S_3H_4S_5M"
```

#### Check if two alleles have a mismatch at a given position

If you pass two allele names and a position to the `isPositionMismatched` function, `hlagenie` will return a boolean indicating whether or not the two alleles have a mismatch at that position.

```python
genie.isPositionMismatched("A*01:01","A*01:02",1) # returns False (the amino acid is the same at position 1 for both alleles)
```

#### Count the number of amino acid mismatches at a position between donor and recipient

The `countAAMismatches` function takes as input four amino acids, two from the donor and two from the recipient (in order). The function returns the number of mismatches between the donor and recipient at that position.

This uses the logic that if a recipient has the donor's amino acid as either of their two amino acids, it is not a mismatch.

```python
genie.countAAMismatches("A","G","G","G") # returns 1
```

#### Count the number of amino acid mismatches between donor and recipient at a given position given the alleles

The `countAAMismatchesAllele` function takes as input four alleles and an amino acid position. The input order is: Donor Allele 1, Donor Allele 2, Recipient Allele 1, Recipient Allele 2, amino acid position. The function returns the number of mismatches between the donor and recipient at that position, adjusting for donor homozygosity.

```python
genie.countAAMismatchesAllele("A*02:01","A*02:01","A*01:01","A*01:01", 44) # returns 1
```