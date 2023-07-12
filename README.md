# HLAGenie
Python package for dealing with HLA sequence data

![hlagenie_logo.png](images/hlagenie_logo.png)

Heavily inspired by and dependent on the fantastic [pyARD](https://www.github.com/nmdp-bioinformatics/py-ard) package. The `HLAGenie` package strives to streamline and standardize the bulk handling of HLA sequence data and provide functions to simplify matching analysis.

## Table of Contents
- [Installation](#installation)
  - [Install from PyPI](#install-from-pypi)
  - [Install from source](#install-from-source)
- [Using HLAGenie](#using-hlagenie)
  - [Using `hlagenie` from Python](#using-hlagenie-from-python)
    - [Initialize `hlagenie`](#initialize-hlagenie)
    - [Accessing sequence dictionaries for HLA alleles](#accessing-sequence-dictionaries-for-hla-alleles)
    - [Retrieve amino acid or nucleotide position from mature protein sequence](#retrieve-amino-acid-or-nucleotide-position-from-mature-protein-sequence)
    - [Retrieve amino acid substring from mature protein sequence](#retrieve-amino-acid-substring-from-mature-protein-sequence)
    - [Retrieve epitope from mature protein sequence](#retrieve-epitope-from-mature-protein-sequence)
    - [Check if two alleles have a mismatch at a given position](#check-if-two-alleles-have-a-mismatch-at-a-given-position)
    - [Count the number of amino acid mismatches at a position between donor and recipient](#count-the-number-of-amino-acid-mismatches-at-a-position-between-donor-and-recipient)
    - [Count the number of amino acid mismatches between donor and recipient at a given position given the alleles](#count-the-number-of-amino-acid-mismatches-between-donor-and-recipient-at-a-given-position-given-the-alleles)
    - [Get the antigen recognition domain sequence of an allele](#get-the-antigen-recognition-domain-sequence-of-an-allele)
    - [Get the extracellular domain sequence of an allele](#get-the-extracellular-domain-sequence-of-an-allele)
  - [Using `hlagenie` from the command line](#using-hlagenie-from-the-command-line)
    - [Retrieval of specific amino acid positions](#retrieval-of-specific-amino-acid-positions)
    - [Retrieval of ARD sequence](#retrieval-of-ard-sequence)
    - [Retrieval of XRD sequence](#retrieval-of-xrd-sequence)
    - [Retrieval of mature protein sequence](#retrieval-of-mature-protein-sequence)
    - [Checking if positions are mismatched between two alleles](#checking-if-positions-are-mismatched-between-two-alleles)
    - [Counting mismatches between donor and recipient given two sets of alleles](#counting-mismatches-between-donor-and-recipient-given-two-sets-of-alleles)
- [Feature requests](#feature-requests)
- [Contributing](#contributing)

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

#### Accessing sequence dictionaries for HLA alleles

The `GENIE` object contains dictionaries of amino acid and nucleotide sequences for each HLA allele. The keys for the dictionaries are the HLA allele names. The values are the genetic sequences.

All of the keys are two-field alleles given the shared protein sequence of these alleles.

```python
genie.seqs # mature protein sequences
genie.full_seqs # full protein sequences
genie.nuc_seqs # full nucleotide sequences
```

#### Retrieve amino acid or nucleotide position from mature protein sequence

To get a given amino acid (or nucleotide) position from a given HLA allele, you can use the `getAA` or `getNuc` functions. These functions are 1-indexed to match standard IMGT/HLA database nomenclature.

For this and following functions, if a three- or four-field allele name is passed, a call is made to `py-ard` to reduce to the two-field level.

```python
genie.getAA("A*01:01",1) # returns "G"
genie.getNuc("A*01:01",1) # returns "A"
```

#### Retrieve amino acid substring from mature protein sequence

To get a given amino acid substring from a given HLA allele, you can use the `getPeptide` function. This function is 1-indexed as well and is inclusive of the start and end positions.

```python
genie.getPeptide("A*01:01",1,10) # returns "GSHSMRYFFT"
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

#### Get the antigen recognition domain sequence of an allele

The `getARD` function takes as input an allele name and returns the antigen recognition domain sequence of that allele.

```python
genie.getARD("A*01:01")
```

#### Get the extracellular domain sequence of an allele

The `getXRD` function takes as input an allele name and returns the extracellular domain sequence of that allele.

```python
genie.getXRD("A*01:01")
```

### Using `hlagenie` from the command line

Some command-line functions are now available.

- Retrieval of specific amino acid positions
- Retrieval of ARD sequence
- Retrieval of XRD sequence
- Retrieval of mature protein sequence

Note that the gapped sequences can be retrieved for any of the following by passing the `--gapped` flag.

#### Retrieval of specific amino acid positions

Calling `hlagenie` from the command line with an allele name and space-delimited list of positions will allow you to retrieve the amino acid at those positions.

```bash
hlagenie -a "A*01:01" -p 1 2 3 4 5 # returns 1G_2S_3H_4S_5M
```

#### Retrieval of ARD sequence

Calling `hlagenie` from the command line with an allele name and the `--ard` flag will allow you to retrieve the ARD sequence of that allele.

```bash
hlagenie -a "A*01:01" --ard
```

#### Retrieval of XRD sequence

Calling `hlagenie` from the command line with an allele name and the `--xrd` flag will allow you to retrieve the XRD sequence of that allele.

```bash
hlagene -a "A*01:01" --xrd
```

#### Retrieval of mature protein sequence

Calling `hlagenie` from the command line with only an allele name will allow you to retrieve the mature protein sequence of that allele.

```bash

hlagenie -a "A*01:01"
```

#### Checking if positions are mismatched between two alleles

Calling `hlagenie-match` from the command line with two allele names and a position will allow you to check if the two alleles have a mismatch at that position. This is based on the mature protein sequence of the alleles. This uses the gapped sequences by default to best assess matching.

```bash
hlagenie-match --allele1 "A*01:01" --allele2 "A*01:02" --positions 1 # returns Matched
```

Supplying a space-delimited list of positions will provide a count of mismatches between the alleles.

```bash
hlagenie-match --allele1 "A*01:01" --allele2 "A*01:02" --positions 1 2 3 4 5 # returns 0
```

Calling `hlagenie-match` from the command line with two allele names and no specified positions will allow you to count the total number of mismatches between the two alleles.

```bash
hlagenie-match --allele1 "A*01:01" --allele2 "A*01:02" # returns 2
```

Alternatively, you can use the `--ard` or `--xrd` flags to check for mismatches in the ARD or XRD sequences, respectively.

```bash
hlagenie-match --allele1 "A*01:01" --allele2 "A*68:02" --ard # returns 23
hlagenie-match --allele1 "A*01:01" --allele2 "A*68:02" --ard # returns 29
```

#### Counting mismatches between donor and recipient given two sets of alleles

Calling `hlagenie-match` from the command line with a recipient genotype and a donor genotype and a set of positions will allow you to retrieve the number of mismatches between the donor and recipient at those positions. This is based on the mature protein sequence of the alleles. This uses the gapped sequences by default to best assess matching. This value is adjusted for donor homozygosity.

```bash
hlagenie-match --recip-haplo "A*01:01+A*01:02" --donor-haplo "A*02:01+A*01:02" --positions 1 2 3 4 5 # returns 0
```

This can also be done with a single position.

```bash
hlagenie-match --recip-geno "A*01:01+A*01:02" --donor-geno "A*02:01+A*01:02" --positions 1 # returns 0
```

Supplying no positions gets the total number of mismatches between the donor and recipient, adjusted for donor homozygosity.

```bash
hlagenie-match --recip-geno "A*01:01+A*01:02" --donor-geno "A*02:01+A*01:02" # returns 32
```

Alternatively, you can use the `--ard` or `--xrd` flags to get the mismatches between the donor and recipient at the ARD or XRD level, respectively.

```bash
hlagenie-match --recip-geno "A*01:01+A*01:02" --donor-geno "A*02:01+A*01:02" --ard # returns 24
hlagenie-match --recip-geno "A*01:01+A*01:02" --donor-geno "A*02:01+A*01:02" --xrd # returns 29
```

## Feature requests

If you have a feature request, please feel free to open a new discussion in the [Ideas](https://github.com/gbiagini/hlagenie/discussions/categories/ideas) page of the Discussions tab. Doing so allows for a more open discussion of the feature and allows others to chime in with their thoughts.

## Contributing

Contributions are welcome. Please feel free to open a pull request with your changes. If you are unsure of how to do this, please feel free to open a discussion in the [Q&A](https://github.com/gbiagini/hlagenie/discussions/categories/q-a). If you are interested in contributing but are unsure of where to start, please check out the [Issues](https://github.com/gbiagini/hlagenie/issues) tab.

Bug reporting is also welcome in the [Issues](https://github.com/gbiagini/hlagenie/issues) tab. Please include as much information as possible, including the version of `hlagenie` you are using, the version of Python you are using, and the operating system you are using. If you are able to provide a minimal reproducible example, that would be very helpful.