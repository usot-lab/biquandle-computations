# Code for virtual knotoid–biquandle computations

This repository contains the code used to compute the invariants
appearing in the paper:

**"Biquandle Virtual Brackets and Virtual Knotoids"**

by Neslihan Gügümcü and Hamdi Kayaslan.

The computations are combinatorially intensive and are not feasible
to carry out by hand.

---

## Requirements

- Python 3
- SymPy
- NumPy

(Developed using Anaconda.)

---

## Files

- `main.py` – command-line interface for running computations
- `functions.py` – core computational routines
- `dataset.py` – virtual knotoid codes, biquandles, and biquandle virtual brackets used in the paper

---

## Usage

List available virtual knotoid codes (labels) and biquandles (indices) in the dataset:
```bash
python main.py count --list
```

Compute the biquandle counting invariant and the biquandle counting matrix invariant of a virtual knotoid with the given biquandle:
```bash
python main.py count --code <LABEL> --biquandle <INDEX>
```
Example :
```bash
python main.py count --code 2.1.1 --biquandle 4
```

List available virtual knotoid codes (labels) and biquandle virtual brackets (indices) in the dataset:
```bash
python main.py bvm_matrix --list
```

Compute the biquandle virtual bracket matrix invariant of a virtual knotoid with respect to the given biquandle virtual bracket:
```bash
python main.py bvb_matrix --code <LABEL> --bvb <INDEX>
```
Example :
```bash
python main.py bvb_matrix --code 2.1.1 --bvb 4
```

Validate a biquandle in the dataset or check if a new data you add to the biquandle storage actually is a valid biquandle:
```bash
python main.py validate --biquandle <INDEX>
```

Validate a biquandle virtual bracket in the dataset or check if a new data you add to the biquandle virtual bracket storage actually is a valid biquandle virtual bracket:
```bash
python main.py validate --bvb <INDEX>
```

Compute biquandle counting invariants of all virtual knotoids in the dataset with a given biquandle:
```bash
python main.py all --biquandle <INDEX>
```

Compute biquandle virtual bracket matrix invariants of all virtual knotoids in the dataset with a given biquandle virtual bracket:
```bash
python main.py all --bvb <INDEX>
```

Examples in the paper :

Example 3.10
```bash
python main.py count --code 2.1.1 --biquandle 7
```

Example 3.13 
```bash
python main.py count --code 2.1.1 --biquandle 6
```
```bash
python main.py count --code 3.1.2 --biquandle 6
```

Example 3.21
```bash
python main.py count --code 2.1.2 --biquandle 6
```
```bash
python main.py count --code 4.1.2 --biquandle 6
```

Example 4.9
```bash
python main.py bvb_matrix --code 2.1.1 --bvb 2
```
```bash
python main.py bvb_matrix --code 4.1.1 --bvb 2
```

Example 4.15
```bash
python main.py count --code 3.1.1 --biquandle 11
```
```bash
python main.py count --code 3.1.3 --biquandle 11
```
```bash
python main.py bvb_matrix --code 3.1.1 --bvb 2
```
```bash
python main.py bvb_matrix --code 3.1.3 --bvb 2
```

Example 4.16 / Table 1
```bash
python main.py count --code 3.1.7 --biquandle 3
```
```bash
python main.py count --code 3.1.9 --biquandle 3
```
```bash
python main.py bvb_matrix --code 3.1.7 --bvb 17
```
```bash
python main.py bvb_matrix --code 3.1.9 --bvb 17
```

Table 2
```bash
python main.py all --bvb 17
```

Notes :

Biquandles and biquandle virtual brackets are indexed by integers
corresponding to their position in the dataset.

## How to run the code

1. Clone the repository:
   ```bash
   git clone https://github.com/usot-lab/biquandle-computations.git
   cd biquandle-computations
   ```
2. Install dependencies:
   ```bash
   pip install sympy numpy
   ```
3. Run a computation:
   ```bash
   python main.py count --list
   ```


