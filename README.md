# Goal
Learn the mathematics of single particle electron cryomicroscopy (cryoEM)

# Resources
Singer, A., & Sigworth, F. J. (2020). Computational Methods for Single-Particle Cryo-EM, 1–40.
Technical University of Denmark. (2020). Introduction to advanced tomography.
Nelson, P. C. (2019). Chapter 12 : Single Particle Reconstruction in Cryo-electron Microscopy. In Physical Models of Living Systems (pp. 305–325).
Singer, A. (2018). Mathematics for cryo-electron microscopy.
Sigworth, F. J. (2016). Principles of cryo-EM single-particle image processing. Microscopy (Oxford, England), 65(1), 57–67. http://doi.org/10.1093/jmicro/dfv370
Jensen, G. J. 2014. Getting Started in Cryo-EM. Retrieved from http://cryo-em-course.caltech.edu/videos
Jensen, G. J. (2010). Preface. In Methods in enzymology (Vol. 482, pp. xv–xvi). http://doi.org/10.1016/S0076-6879(10)82018-X
Leong, P. A., Yu, X., Zhou, Z. H., & Jensen, G. J. (2010). Correcting for the Ewald Sphere in High-Resolution Single-Particle Reconstructions. In Methods in Enzymology (1st ed., Vol. 482, pp. 369–380). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82015-4
Scheres, S. H. W. (2010). Classification of Structural Heterogeneity by Maximum-Likelihood Methods. In Methods in Enzymology (1st ed., Vol. 482, pp. 295–320). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82012-9
Cong, Y., & Ludtke, S. J. (2010). Single Particle Analysis at High Resolution. In Methods in Enzymology (1st ed., Vol. 482, pp. 211–235). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82009-9
Penczek, P. A. (2010). Resolution Measures in Molecular Electron Microscopy. In Methods in Enzymology (1st ed., Vol. 482, pp. 73–100). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82003-8
Penczek, P. A. (2010). Fundamentals of Three-Dimensional Reconstruction from Projections. In Methods in Enzymology (1st ed., Vol. 482, pp. 1–33). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82001-4
Sigworth, F. J., Doerschuk, P. C., Carazo, J.-M., & Scheres, S. H. W. (2010). An Introduction to Maximum-Likelihood Methods in Cryo-EM. In Methods in Enzymology (1st ed., Vol. 482, pp. 263–294). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82011-7
Penczek, P. A. (2010). Image restoration in cryo-electron microscopy. Methods in Enzymology (1st ed., Vol. 482). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82002-6
Glaeser, R. (2007). Electron Crystallography of Biological Macromolecules. Oxford University Press.
Frank, J. (1996). Electron Microscopy of Macromolecular Assemblies. In Three-Dimensional Electron Microscopy of Macromolecular Assemblies (pp. 12–53). Elsevier. http://doi.org/10.1016/B978-012265040-6/50002-3

# Why?
We can see living atoms with electrons. Electron microscopes use hundreds of thousands of volts to speed up single electrons to three quarters the speed of light. At such high speeds electrons have picometer wavelengths and can resolve the distances between atoms in biomolecules. Electron cryomicroscopy (cryoEM) won the 2017 Nobel prize in Chemistry, and pharmaceutical companies have invested in this technology for applications like structure based drug design. 

After biochemical sample preparation of a purified biomolecule, two dimensional images are collected on electron microscopes. Images of single biomolecules are very noisy and computer algorithms average tens of thousands to millions of 2D images to reconstruct a 3D discrete scalar map that represents the Coulomb density (what the electron feels). 

Historically cryoEM computational workflows have drawn from digital signal and image processing theory. With the recent popularization of cryoEM, many researchers have entered the field who do not have a strong physics/engineering background and treat algorithms as black boxes. Unfortunately, this can limit their intuition of what data processing strategies might work, or how to troubleshoot when confusing results are generated computationally. Furthermore, there is a growing trend of applied mathematicians, statisticians and computer scientists to approach cryoEM using the familiar methods of their discipline, providing novel breakthroughs. Therefore, a more transparent and pedagogical treatment of the mathematics of single particle electron cryomicroscopy would be welcomed by researchers who analyze cryoEM data and/or develop computational approaches to cryoEM data.

The mathematics of cryoEM spans several disciplines, from the physics of electron microscopes, through digital Fourier transformations, to Bayesian inference and beyond. This broad range of specialties makes the mathematics of cryoEM a challenge to master, and open to contributions from many fields.

In this repo I would like to open up the black box of cryoEM computation, and exhibit the mathematical objects contained inside. Those sufficiently curious and motivated can teach themselves more about cryoEM data processing by playing with simple models in interactive programming notebooks, made publicly available here.
