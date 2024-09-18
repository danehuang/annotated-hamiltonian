# The Annotated Hamiltonian 

This series of notebooks introduces Analog Quantum Computing (QC), including related concepts of Adiabatic QC and Neutral Atom QC.

```

    Gate-Based QC  <---------------|                ^
(Unitary Transformation)           |                |
            ^                      |             Abstract
 ===========|======================|=======================
            |                      |             Physical
            v                      v                |
        Analog QC <----------- Adiabatic QC         v
 (Schrödinger's Equation)  (Adiabatic Theorem)
            ^                      ^
            |-----|          |-----|
                  |          |
                  Neutral Atom
              (Rydberg Hamiltonian)

```


## Installation

```
python -m venv annotated-hamiltonian
source annotated-hamiltonian/bin/activate
pip install -r requirements.txt
```


## Citation

If you found our work useful, please consider citing it:

```
@misc{huang2024annotatedhamiltonian,
  author       = {Daniel Huang},
  title        = {The Annotated Hamiltonian},
  year         = {2024},
  url          = {https://github.com/danehuang/annotated-hamiltonian},
  note         = {GitHub repository},
}
```
