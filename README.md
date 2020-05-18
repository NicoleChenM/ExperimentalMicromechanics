# Experimental Micromechanics

Python tools for experimental micromechanics. Tools for:
- SEM images
- Nanoindentation
- etc.

Documentation, Tutorials are at [github-pages](https://steffenbrinckmann.github.io/ExperimentalMicromechanics/)

# Installation (just files)
1. Requirements
numpy, matplotlib, scipy, Pillow, scikit-image, opencv-python
```
pip3 install numpy matplotlib scipy Pillow scikit-image opencv-python
```
2. Download files from 'src' folder and use


# Installation (everything)
1. Requirements
numpy, matplotlib, scipy, Pillow, scikit-image, opencv-python
```
pip3 install numpy matplotlib scipy Pillow scikit-image opencv-python
```
2. clone this gitlab repository
```
git clone git@jugit.fz-juelich.de:s.brinckmann/experimetal-micromechanics.git
```
3. run doctest to verify that everything works.
E.g. to test the Tif class
```
python3 -m doctest -o 'NORMALIZE_WHITESPACE' Tif.doctest
```
4. create documentation (optional)
  - install doxygen
  - install doxypypy "pip3 install doxypypy"
  - run "python3 buildDoxy.py"
  - run "doxygen doxygen/Doxyfile"

# Contact information
```
Dr. Steffen Brinckmann
Micromechanics and Tribology
Microstructure and Properties of Materials (IEK-2)
Forschungszentrum Jülich
Tel.: +49 2461 61-4425
```