# Goal
Learn the mathematics of single particle electron cryomicroscopy (cryoEM)

# Resources
## cryoEM
* Singer, A., & Sigworth, F. J. (2020). Computational Methods for Single-Particle Cryo-EM, 1–40.
* Nelson, P. C. (2019). Chapter 12 : Single Particle Reconstruction in Cryo-electron Microscopy. In Physical Models of Living Systems (pp. 305–325).
* Singer, A. (2018). Mathematics for cryo-electron microscopy.
* Grant, J. (2018). Getting Started in Cryo-EM. http://cryo-em-course.caltech.edu/videos
* Glaeser, R. (2007). Electron Crystallography of Biological Macromolecules. Oxford University Press.
* Frank, J. (1996). Electron Microscopy of Macromolecular Assemblies. In Three-Dimensional Electron Microscopy of Macromolecular Assemblies (pp. 12–53). Elsevier. http://doi.org/10.1016/B978-012265040-6/50002-3

## Math
* Hobbie, R. K., & Roth, B. J. (2007). Intermediate Physics for Medicine and Biology (4th ed.). New York, NY: Springer New York. http://doi.org/10.1007/978-0-387-49885-0
  * Chapter 11: Method of Least Squares and Signal Analysis
  * Chapter 12: Images

# Why?
We can see living atoms with electrons. Electron microscopes use hundreds of thousands of volts to speed up single electrons to three quarters the speed of light. At such high speeds electrons have picometer wavelengths and can resolve the distances between atoms in biomolecules. Electron cryomicroscopy (cryoEM) won the 2017 Nobel prize in Chemistry, and pharmaceutical companies have invested in this technology for applications like structure based drug design. 

After biochemical sample preparation of a purified biomolecule, two dimensional images are collected on electron microscopes. Images of single biomolecules are very noisy and computer algorithms average tens of thousands to millions of 2D images to reconstruct a 3D discrete scalar map that represents the Coulomb density (what the electron feels). 

Historically cryoEM computational workflows have drawn from digital signal and image processing theory. With the recent popularization of cryoEM, many researchers have entered the field who do not have a strong physics/engineering background and treat algorithms as black boxes. Unfortunately, this can limit their intuition of what data processing strategies might work, or how to troubleshoot when confusing results are generated computationally. Furthermore, there is a growing trend of applied mathematicians, statisticians and computer scientists to approach cryoEM using the familiar methods of their discipline, providing novel breakthroughs. Therefore, a more transparent and pedagogical treatment of the mathematics of single particle electron cryomicroscopy would be welcomed by researchers who analyze cryoEM data and/or develop computational approaches to cryoEM data.

The mathematics of cryoEM spans several disciplines, from the physics of electron microscopes, through digital Fourier transformations, to Bayesian inference and beyond. This broad range of specialties makes the mathematics of cryoEM a challenge to master, and open to contributions from many fields.

In this repo I would like to open up the black box of cryoEM computation, and exhibit the mathematical objects contained inside. Those sufficiently curious and motivated can teach themselves more about cryoEM data processing by playing with simple models in interactive programming notebooks, made publicly available here.
