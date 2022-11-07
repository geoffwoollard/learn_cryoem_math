# Goal
Learn the mathematics of single particle electron cryomicroscopy (cryoEM)

# Why?
We can see living atoms with electrons. Electron microscopes use hundreds of thousands of volts to speed up single electrons to three quarters the speed of light. At such high speeds electrons have picometer wavelengths and can resolve the distances between atoms in biomolecules. Electron cryomicroscopy (cryoEM) won the 2017 Nobel prize in Chemistry, and pharmaceutical companies have invested in this technology for applications like structure based drug design. 

After biochemical sample preparation of a purified biomolecule, two dimensional images are collected on electron microscopes. Images of single biomolecules are very noisy and computer algorithms average tens of thousands to millions of 2D images to reconstruct a 3D discrete scalar map that represents the Coulomb density (what the electron feels). 

Historically cryoEM computational workflows have drawn from digital signal and image processing theory. With the recent popularization of cryoEM, many researchers have entered the field who do not have a strong physics/engineering background and treat algorithms as black boxes. Unfortunately, this can limit their intuition of what data processing strategies might work, or how to troubleshoot when confusing results are generated computationally. Furthermore, there is a growing trend of applied mathematicians, statisticians and computer scientists to approach cryoEM using the familiar methods of their discipline, providing novel breakthroughs. Therefore, a more transparent and pedagogical treatment of the mathematics of single particle electron cryomicroscopy would be welcomed by researchers who analyze cryoEM data and/or develop computational approaches to cryoEM data.

The mathematics of cryoEM spans several disciplines, from the physics of electron microscopes, through digital Fourier transformations, to Bayesian inference and beyond. This broad range of specialties makes the mathematics of cryoEM a challenge to master, and open to contributions from many fields.

In this repo I would like to open up the black box of cryoEM computation, and exhibit the mathematical objects contained inside. Those sufficiently curious and motivated can teach themselves more about cryoEM data processing by playing with simple models in interactive programming notebooks, made publicly available here.

**NB: If there is a particular computational aspect of cryoEM analysis that you are curious about, feel free to request (via [Issues](https://github.com/geoffwoollard/learn_cryoem_math/issues) in this repo) some learning materials. I will do my best to develop some code and data to play with for interactive learning.**

# Resources
## Online Resources
* Jensen, G. J. (2014). Getting Started in Cryo-EM. Retrieved from http://cryo-em-course.caltech.edu/videos
   * Recommended place to start for total beginners. Coursera style.
* MRC-LMB Electron Cryo-microscopy Course 2017 (recorded lectures). Retrieved from https://www2.mrc-lmb.cam.ac.uk/research/scientific-training/electron-microscopy/
* Jensen, G. J., & Vos, M. (2020). Electron Microscopy University. Retrieved from https://em-learning.com/
  * 70 hours of videos. From intro to very advanced. Theoretical white board lectures and videos of someone operating the microscope. Includes important concepts hard to find elsewhere (vacuum theory, heuristics of data collection / experimental design).
* Sigworth, F., & Tagare, H. (2020). Cryo-EM Principles. Retrieved from https://cryoemprinciples.yale.edu/
  * Pedagogical step-by-step introduction to math
* Shen, P., Iwasa, J., Thuesen, A., Thoms, J., & Mayson, P. (2020). CryoEM 101. Retrieved from https://cryoem101.org/
  * Visual animations help build intuition

  
## CryoEM Textbooks
* Frank, J. (1996). Electron Microscopy of Macromolecular Assemblies. In Three-Dimensional Electron Microscopy of Macromolecular Assemblies (pp. 12–53). Elsevier. http://doi.org/10.1016/B978-012265040-6/50002-3
* Glaeser, R. (2007). Electron Crystallography of Biological Macromolecules. Oxford University Press.
   * Good historical overview and mathematical excursions. Narrative style. Works out analytical calculations. Before the major gains in resolution, and emphasis on 2D crystals.
* Glaeser, R. M., Nogales, E., & Chiu, W. (Eds.). (2021). Single-particle Cryo-EM of Biological Macromolecules. IOP Publishing. http://doi.org/10.1088/978-0-7503-3039-8
  * 1 Introduction and overview
    * 1.1 Visualizing biological molecules to understand life’s principles. Gives an honest overview of the current state of the art by citing recent studies and putting them in historical perspective.
    * 1.2 Recovery of 3D structures from images of weak-phase objects. Commentary on weak phase and projection approximations and explanation of what the mathematical equations are physically modeling the image. Comparisons to X-ray crystalography.
  * 2 Sample preparation
    * 2.2 Initial screening of samples in negative stain. Includes detailed protocols for the "black-art" of negative staining, and advice on what can be interpreted about the sample from data processing.
    * 2.3 Standard method of making grids for cryo-EM. Grids, glow-discharging, freezing.
    * 2.4 Requirement to make very thin specimens for cryo-EM. Inelastic scattering, too thick/thin. Useful advice for screening and data collection. Helpful to build intuition for what is physically happening to the sample and why this prevents good quality data from being collected.
    * 2.5 Current strategies for optimizing preparation of cryo-grids. Air-water interface and physics of thin films. How to troubleshoot. Review of different freezing hardware devices.
  * 3 Data collection
    * 3.2 Radiation damage in cryo-EM. Sources of radiation damage, role of electron energy, exposure weighting, beam-induced motion, charging. This chapter helps untangle these related ideas and cites the empirical studies that support the current consensus about what is physically causing damage and how it quantitatively behaves.
    * 3.3 Low-dose protocols for recording images. Automatic data collection, grid atlas, focusing, etc.
    * 3.4 Practical considerations: defocus, stigmation, coma-free illumination, and phase plates. Focus on phase plates. Loss of information from defocus explained with a simple figure. Literature cited on volta and laser phase plate development.
    * 3.5 Practical considerations: movie-mode data acquisition. Magnification, dose rate, total dose. Brief overview.
  * 4 Data processing
    * 4.2 Automated extraction of particles. Going from micrographs to particles. Challenges of particle picking: manual, template, neural network. Extraction. Cleanup with 2D class and 3D ab initio.
    * 4.3 CTF estimation and image correction (restoration). Mathematical details of CTF equations, including astigmatism. Estimating the CTF from experimental images. Correcting images with Wiener-filtering. Brief remarks on higher order corrections.
    * 4.4 Merging data from structurally homogeneous subsets. How many particles are needed: in ideal theoretical cases and in practice. Review of 3D reconstruction: algebraic, Fourier methods based on central-slice theorem, tilting, common lines, random initialization, validation.
    * 4.5 3D classification of structurally heterogeneous particles. Brief summary of masking, symmetry, continuous motion.
    * 4.6 Preferred orientation: how to recognize and deal with adverse effects. Causes, effects, and mitigation of preferred orientation. Ample literature cited. Links to software tools developed by the author to quantitate preferred orientations, their effect on resolution, etc.
    * 4.7 B factors and map sharpening. Commentary on what causes the B-factor and what it is modeling. How to estimate and advice on sharpening in practice. Table of diverse reporting of B factors in literature.
    * 4.8 Optical aberrations and Ewald sphere curvature. Symmetrical and antisymmetrical components of the phase shift, γ. Order of aberrations, Zernike polynomials: astigmatism and defocus; axial coma; trefoil or three-fold astigmatism; spherical aberration and; four-fold astigmatism or tetrafoil. Practical considerations when estimating, including how estimates of different orders influence each other. At what point a single molecule is thick enough to feel the curvature of the Ewald sphere.
  * 5 Map validation
    * 5.2 Measures of resolution: FSC and local resolution. FSC and SSNR, pathologies to FSC from masking, limited defoci diversity. Local resolution. Instructive chapter helping interpret statistical meaning of metrics and what sort of conclusions they can and can't justify.
    * 5.3 Recognizing the effect of bias and over-fitting. Build-up of extraneous map features due to aligned noise during iterative refinement. Effect on FSC. Local overfitting and what map artefacts look like. New algorithmic data processing solutions and why/how they work at a high level.
    * 5.4 Estimates of alignment accuracy. SNR dependency on masking and box size, investigating alignment errors with synthetic data, effect of alignment accuracy on resolution. Some basic quantitative treatment of alignment errors, a topic that is not discussed much since there is no ground truth alignment of experimental data.
  * 6 Model building and validation
    * 6.2 Using known components or homologs: model building. Review of software for medium to low resolution maps: rigid body, flexible for resolutions worse than ~4 A.
    * 6.3 Building atomistic models in cryo-EM density maps. Flow chart of practical advice. Overview of ~22 software tools, with a few sentences to a paragraph written about each. Fit to density, map-model FSC.
    * 6.4 Quality evaluation of cryo-EM map-derived models. Overview of different map-model metrics: map-model FSC, atom inclusion, average density value, cross correlation, Z-score, EMRinger, Q-score. Q-score resolution dependency and per amino acid / nucleotide base, water molecules and ions. Model only metrics.
    * 6.5 How algorithms from crystallography are helping electron cryo- microscopy. Review of phenix software and how they can be modularly used with each other in map improvement, map interpretation, model building, model optimization, validation. Tools include phenix.auto_sharpen, phenix.resolve_cryo_em, phenix.combine_focused_maps, phenix.auto_sharpen, phenix.real_space_refine, phenix.superpose_maps, phenix.dock_in_map, phenix.map_to_model, phenix.trace_and_build, phenix.combine_models, phenix.sequence_from_map, phenix.map_box, phenix.mtriage, phenix.comprehensive_validation, phenix.reduce, phenix.probe, phenix.kinemage, phenix.mr_rosetta. Lot's of literature cited and useful heuristics given.
    * 6.6 Archiving structures and data. EMPIAR, EMDB, PDB. Data validation. Sample sequence and ligands. OneDep. Metadata.

## Math

### Rotations
* Quine, J. R. (n.d.). Mathematical techniques in structural biology, 1–86.
  * Good gentle introduction from structural biologists on orthonormal rotations. Derivation of familiar properties of rotations. Exercises.
* Hartley, R., Trumpf, J., Dai, Y., & Li, H. (2013). Rotation Averaging. International Journal of Computer Vision, 103(3), 267–305. http://doi.org/10.1007/s11263-012-0601-0
  * Long review on rotational encodings. Encodings are mainly $\mathcal{R}^{3x3}$, (axis-angle, quaternions). SO(3) is discussed from concepts from Lie groups, and differential geometry. Distance measures are derived for different encodings and focus on the Frobenius norm of the coordinates in that encoding, and the angle/gendesic. Three basic computer vision problem types are considered: 1. averaring multiple rotations into a single "average rotation", 2. estimating how to convert between frames from multiple rotations pairs (e.g. from calibrated cameras), and 3. estimating mutliple rotations from difference frames from sparse relative rotations between paris of frames.

### Fourier
* Lighthill, M. J. (1958). An Introduction to Fourier Analysis and Generalised Functions. An Introduction to Fourier Analysis and Generalised Functions. Cambridge University Press. http://doi.org/10.1017/CBO9781139171427.001
  * subtelties about delta functions. Helps for not getting confused between discrete and continous FT 
* Vembu, S. (1961). Fourier Transformation of the n-Dimensional Radial Delta Function. The Quarterly Journal of Mathematics, 12(1), 165–168. http://doi.org/10.1093/qmath/12.1.165
* Briggs, W. L., & Hensen, V. E. (1995). The DFT: an owner’s manual for the discrete Fourier transform (1st ed.). SIAM.
  * very pedagogical textbook on DFTs. Starts with gentle introduction and goes to advanced level.
* Goodman, J. W. (1996). Introduction to Fourier Optics (2nd ed.).
  * Classic optics text. Incudes end of chapter problems. See math chapter, Chapter 2: Analysis of Two-Dimensional Signals and Systems.
* Bracewell, R. (2000). The Fourier Transform And Its Applications. The Fourier Transform and It’s Applications (3rd ed.). http://doi.org/10.1017/cbo9780511623813.014
* Khuri, A. I. (2004). Applications of Dirac’s delta function in statistics. International Journal of Mathematical Education in Science and Technology, 35(2), 185–195. http://doi.org/10.1080/00207390310001638313
* Hobbie, R. K., & Roth, B. J. (2007). Intermediate Physics for Medicine and Biology (4th ed.). New York, NY: Springer New York. http://doi.org/10.1007/978-0-387-49885-0
  * Discussion of convolution theorems and slice theorems. Includes very thoughtful end of chapter problems. Author maintains blog where he works out some problems in detail.
* Osgood, B. (2007). Lecture Notes for EE 261: The Fourier Transform and its Applications. Lecture Notes for EE 261 - The Fourier Transform and its Applications.
* Filedman, J. (2007). Discrete-Time Fourier Series and Fourier Transform. Retrieved from https://www.math.ubc.ca/~feldman/m267/dft.pdf
* Gonzalez, R. C., & Woods, R. E. (2008). Digital Image Processing (2nd ed.).
  * Classic image analysis chapter. Chapter of Fourier filtering helps build intuition of how to use DFT/FFT in practice.
* Lam, E. (2008). Some Properties of Fourier Transform 1.
* Wang, Q., Ronneberger, O., & Burkhardt, H. (2008). Fourier Analysis in Polar and Spherical Coordinates. Technical Report, University of Freiburg, (Internal Report 1/08).
  * Advances, but details (not terse) mathamatical analysis of 2D and 3D transforms (Fourier, Zernike, Spherical Harmonics, etc).
* Jeong, D. (2010). Appendix A: Fourier transforms. In Cosmology with high (z > 1) redshift galaxy surveys (pp. 238–253). CRC Press.
* Bright, A., & Wang, R. (2010). Delta Functions Generated by Complex Exponentials. Retrieved from http://fourier.eng.hmc.edu/e102/lectures/ExponentialDelta.pdf
* Baddour, N. (2011). Two-Dimensional Fourier Transforms in Polar Coordinates. In Advances in Imaging and Electron Physics (Vol. 165, pp. 1–45). Elsevier Inc. http://doi.org/10.1016/B978-0-12-385861-0.00001-4
* Berendsen, H. J. C. (2012). Fourier transforms. In Simulating the Physical World (pp. 315–334). Cambridge: Cambridge University Press. http://doi.org/10.1017/CBO9780511815348.014
* Henning, A. J., Huntley, J. M., & Giusca, C. L. (2015). Obtaining the Transfer Function of optical instruments using large calibrated reference objects. Optics Express, 23(13), 16617. http://doi.org/10.1364/OE.23.016617
* Haber, H. E. (2018). A Gaussian integral with a purely imaginary argument. Retrieved from http://scipp.ucsc.edu/~haber/ph215/Gaussian.pdf
* Baddour. (2019). Discrete Two-Dimensional Fourier Transform in Polar Coordinates Part I: Theory and Operational Rules. Mathematics, 7(8), 698. http://doi.org/10.3390/math7080698
* Yao, X., & Baddour, N. (2020). Discrete two dimensional Fourier transform in polar coordinates part II: numerical computation and approximation of the continuous transform. PeerJ Computer Science, 6, e257. http://doi.org/10.7717/peerj-cs.257

### Poisson (stats, pmf, random variable)
* Khuri, A. I. (2004). Applications of Dirac’s delta function in statistics. International Journal of Mathematical Education in Science and Technology, 35(2), 185–195. http://doi.org/10.1080/00207390310001638313
* Ward, M. D. (2005). Sums of independent Poisson random variables are Poisson random variables. Retrieved from https://llc.stat.purdue.edu/2014/41600/notes/prob1805.pdf
* Foi, A., Trimeche, M., Katkovnik, V., & Egiazarian, K. (2008). Practical Poissonian-Gaussian noise modeling and fitting for single-image raw-data. IEEE Transactions on Image Processing, 17(10), 1737–1754. http://doi.org/10.1109/TIP.2008.2001399
* Peacock, J. (2012). Astronomical Statistics. Retrieved from https://www.roe.ac.uk/japwww/teaching/astrostats/astrostats2012_part2.pdf
* Loh, N. D. (2014). A minimal view of single-particle imaging with X-ray lasers. Philosophical Transactions of the Royal Society B: Biological Sciences, 369(1647), 20130328. http://doi.org/10.1098/rstb.2013.0328
  * Poisson stats used to model XFEL. Includes alignment algorithm and background noise vs signal.
* Hasinoff, S. W. (2014). Photon, Poisson Noise. Computer Vision, 608–610. http://doi.org/10.1007/978-0-387-31439-6_482
* Foi, A. (2014). ICIP 2014 Tutorial T7 : Signal-Dependent Noise and Stabilization of Variance. Retrieved from http://www.cs.tut.fi/~foi/papers/icip2014_tutorial_t7_sigdepvst-foi-single.pdf
* Liu, L. T., Dobriban, E., & Singer, A. (2018). ePCA: High dimensional exponential family PCA. Annals of Applied Statistics, 12(4), 2121–2150. http://doi.org/10.1214/18-AOAS1146
* Kardar, M. (2019). II. Probability. Retrieved from http://web.mit.edu/8.333/www/lectures/lec5.pdf
* Grob, P., Bean, D., Typke, D., Li, X., Nogales, E., & Glaeser, R. M. (2013). Ranking TEM cameras by their response to electron shot noise. Ultramicroscopy, 133, 1–7. http://doi.org/10.1016/j.ultramic.2013.01.003
  * Details of stats of DFT/FFT of noise image. Poisson counting stats assumed, in high electron count limit, where the detector detects the energy from the counts, and each depositing event (measured contribution to the pixel value) is a general random variable with only the mean and variance known/fixed.
* Sauter, N. K., Kern, J., Yano, J., & Holton, J. M. (2020). Towards the spatial resolution of metalloprotein charge states by detailed modeling of XFEL crystallographic diffraction. Acta Crystallographica Section D: Structural Biology, 76, 176–192. http://doi.org/10.1107/S2059798320000418
* Jensen, G. J., & Vos, M. (2020). Electron Microscopy University. Retrieved from https://em-learning.com/
  * includes a whole section of videos on detectors. A good introduction to detector physics.
*  Brunton, Steve. (2021). Fourier Analysis. https://youtube.com/playlist?list=PLMrJAkhIeNNT_Xh3Oy0Y4LTj0Oxo8GqsC
  * Youtube play list. Fourier basics. Code in matlab and python. Fourier series. Complex Fourier series. Gibbs. Convolution. Parseval's. DFT. FFT. Denoising with FFT. Derivatives wtih FFT. PDEs wtih FFT. Spectrogram. Uncertainty. Image compression. Wavelets. Laplace.

### Linear Algebra
* Abdi, H. (2007). The Eigen-Decomposition: Eigenvalues and Eigenvectors. Encyclopedia of Measurements and Statistics, 1–10.
  * Helpful for math in Grob et al (2013)

### Geometry
* eigenchris (2017). Tensors For Beginners. Youtube. https://youtu.be/8ptMTLzV4-I
  * 19 part video seris (each ~ 10 min) on tensors, forward and backward transformations, vectors, covectors, linear maps, metric tensors, bilinear forms, tensor products, kronecker product, raising/lowering indexes

### Complex Random Variables
* Halliwell, L. J. (2015). Complex random variables, 1–66.
* Miller, S. L., & Childers, D. (2012). Pairs of Random Variables. In Probability and Random Processes (pp. 177–243). Elsevier. http://doi.org/10.1016/B978-0-12-386981-4.50008-4
* Silvestrov, D. (2014). Lecture 10 : Characteristic Functions.
* Eriksson, J., Ollila, E., & Koivunen, V. (2009). Statistics for complex random variables revisited. In 2009 IEEE International Conference on Acoustics, Speech and Signal Processing (Vol. 2, pp. 3565–3568). IEEE. http://doi.org/10.1109/ICASSP.2009.4960396
* Wikipedia. (2020). Complex Random Variable. Retrieved from https://en.wikipedia.org/wiki/Complex_random_variable

### Machine Learning (incl. general statistics / probability theory & related topics)
* Duda RO, Hart PE, Stork DG. (2000) Pattern Classification. John Wiley & Sons.
  * a somewhat dated view on machine learning & pattern recognition, but nevertheless a good introduction to the basic algorithms. 
* MacKay DJC. (2003) Information Theory, Inference, and Learning Algorithms. Cambridge University Press. http://www.inference.org.uk/mackay/itila/book.html (free PDF download)
* Bishop CM. (2006) Pattern Recognition and Machine Learning. Springer Berlin / Heidelberg
* Hastie T, Tibshirani R, Friedman J. (2008) The Elements of Statistical Learning: Data Mining, Inference, and Prediction. https://web.stanford.edu/~hastie/ElemStatLearn/ (free PDF download)

## Physics (electon optics, lenses, detector)
* Shaw, R. (1978). Evaluating the efficient of imaging processes. Reports on Progress in Physics, 41(7), 1103–1155. http://doi.org/10.1088/0034-4885/41/7/003
* Rabbani, M., Van Metter, R., & Shaw, R. (1987). Detective quantum efficiency of imaging systems with amplifying and scattering mechanisms. Journal of the Optical Society of America A, 4(5), 895. http://doi.org/10.1364/JOSAA.4.000895
* Snyder, D. L., White, R. L., & Hammoud, A. M. (1993). Image recovery from data acquired with a charge-coupled-device camera. Journal of the Optical Society of America A, 10(5), 1014. http://doi.org/10.1364/JOSAA.10.001014
* Vliet, L. J. Van, Boddeke, F. R., Sudar, D., & Young, I. T. (1998). Image Detectors for Digital Image Microscopy. In Digital Image Analysis of Microbes; Imaging, Morphometry, Fluorometry and Motility Techniques and Applications (pp. 37–64).
  * readout noise as Gaussian noise from electronics. semi techincal explanation of how CCDs work. pedagogical outline of sources of noise: 1. Readout noise: on-chip electronics, pre-amplifier and KTC noise (Gaussian additive) 2. Quantization noise (uniform, additive) 3. Photon (ie shot) (Poisson) 4. Thermal (ie dark current and hot pixels) (Poisson). See Table 3 for their distribution, dependencies, SNR, and other remakrs. Detailed discussion of gain, but likely specific to CCD cameras.
* Meyer, R. R., & Kirkland, A. (1998). The effects of electron and photon scattering on signal and noise transfer properties of scintillators in CCD cameras used for electron detection. Ultramicroscopy, 75(1), 23–33. http://doi.org/10.1016/S0304-3991(98)00051-5
* Faruqi, A. R., Henderson, R., & Tlustos, L. (2005). Noiseless direct detection of electrons in Medipix2 for electron microscopy. Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, 546(1–2), 160–163. http://doi.org/10.1016/j.nima.2005.03.022
  * Experimental characterization of Medipix2. Dose as low as < 0.01 e/pixel. Counting stats "essentially Poissonian". Discussions of MTF.
* Faruqi, A. R., Henderson, R., Pryddetch, M., Allport, P., & Evans, A. (2005). Direct single electron detection with a CMOS detector for electron microscopy. Nuclear Instruments and Methods in Physics Research, Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, 546(1–2), 170–175. http://doi.org/10.1016/j.nima.2005.03.023
  * Experimental characterization of MAPS detector. Single electon events detected and pixel values of nearby pixels studied: profile of sorrounding pixels, SNR.
* Evans, D. A., Allport, P. P., Casse, G., Faruqi, A. R., Gallop, B., Henderson, R., … Waltham, N. (2005). CMOS active pixel sensors for ionising radiation. Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, 546(1–2), 281–285. http://doi.org/10.1016/j.nima.2005.03.108
  * Plots of "Electron Transfer Curve", ie deposited energy per electron event.
* Kirkland, E. J. (2006). Image Simulation in Transmission Electron Microscopy, 1–14.
* Glaeser, R. (2007). Electron Crystallography of Biological Macromolecules. Oxford University Press.
* McMullan, G., Cattermole, D. M., Chen, S., Henderson, R., Llopart, X., Summerfield, C., … Faruqi, A. R. (2007). Electron imaging with Medipix2 hybrid pixel detector. Ultramicroscopy, 107(4–5), 401–413. http://doi.org/10.1016/j.ultramic.2006.10.005
  * Experimental characterizatoin of the Medipix2, originally a high energy particle physics detector, and repositioned to detection electrons in TEM. They derive a simple and useful analytical result of the DQE (at zero spatial frequency) for counting detectors. It perhaps could be modified for other intermadiate spatial frequencies. They investigate the effect of threshold energy (to call a count) on the DQE with Monte Carlo simulations, and suggest a threshold of half the incident energy (e.g. 100 keV for a 200 keV scope). Dong MC simulations allows them to investigate how many pixels are being miscounted through backscattering, or would not be observed experimentally (the number of zero count events).
* Turchetta, R. (2007). CMOS monolithic active pixel sensors (MAPS) for scientific applications: Some notes about radiation hardness. Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, 583(1), 131–133. http://doi.org/10.1016/j.nima.2007.08.226
  * review the state-of-art of MAPS as imaging sensors. historical perspective. consumer and scientific applications. references to electronics literature.
* Mooney, P. (2009). A Noise-Tolerant Method for Measuring MTF from Found-Object Edges in a TEM. Microscopy and Microanalysis, 15(S2), 234–235. http://doi.org/10.1017/S1431927609098407
* Baxter, W. T., Grassucci, R. A., Gao, H., & Frank, J. (2009). Determination of signal-to-noise ratios and spectral SNRs in cryo-EM low-dose imaging of molecules. Journal of Structural Biology, 166(2), 126–132. http://doi.org/10.1016/j.jsb.2009.02.012
* McMullan, G., Chen, S., Henderson, R., & Faruqi, A. R. (2009). Detective quantum efficiency of electron area detectors in electron microscopy. Ultramicroscopy, 109(9), 1126–1143. http://doi.org/10.1016/j.ultramic.2009.04.002
  * Rigerous and pedagogical (and not too long) overview of how to measure modulation transfer function (MTF), noise power spectrum (NPS), detective quantum efficiency (DQE) from images (knife method). Analytical results for edge spread function, line spread function given which are derived/presented in more detail in references (textbooks, etc); see equations 10-13. special care is taken for DQE(w=0), because there is much noise in the NPS at low spatial frequency and the MTF can drop off quickly at low spatial frequency. Methods to measure the NPS from the difference of digitized frames is presented. All this theory is used to compare film, CCD, Medipix2 (hybrid pixel detector) and a MAPS (monolithic active pixel sensor) detector. For detectors (such as MAPS) that can distinguish individual electron events, one can calculate the DQE from the probability density of the energy (Landau plot), although this density is hard to estimate. The authors do Monte Carlo simulations at two levels of theory (“continuous slowing down approximation”, and another they call “Full Monte Carlo” that treats both elastic and inelastic collisions as stochastic), and give references to a previous study and a textbook. Using this theory they simulate electron trajectories, deposited energy, and the resulting image in film, and MAPS. They are then able connect detector performance (DQE, etc) to physical/causal explanations, and suggest how to make the ideal detector. History has vindicated their calls for how to make a detector of choice.
* Rullgård, H., Öfverstedt, L.-G., Masich, S., Daneholt, B., & Öktem, O. (2011). Simulation of transmission electron microscope images of biological specimens. Journal of Microscopy, 243(3), 234–256. http://doi.org/10.1111/j.1365-2818.2011.03497.x
  * Early forward model of image formation. Includes: Poisson shot noise, inverse Gaussian (or Wald) for the Landau distribution for the amount of energy an electron deposits.
* Ghadimi, R., Daberkow, I., Kofler, C., Sparlinek, P., & Tietz, H. (2011). Characterization of 16 MegaPixel CMOS Detector for TEM by Evaluating Single Events of Primary Electrons. Microscopy and Microanalysis, 17(S2), 1208–1209. http://doi.org/10.1017/S143192761100691X
* Ruskin, R. S., Yu, Z., & Grigorieff, N. (2013). Quantitative characterization of electron detectors for transmission electron microscopy. Journal of Structural Biology, 184(3), 385–393. http://doi.org/10.1016/j.jsb.2013.10.016
  * Details on how to experimentally measure the DQE for different detectors (CMOS, CCD, film). Mathematical details about modulation transfer function, noise poser spectrum, signal to noise ratio
* Vulović, M., Ravelli, R. B. G., van Vliet, L. J., Koster, A. J., Lazić, I., Lücken, U., Rullgård, H., Öktem, O., Rieger, B. (2013). Image formation modeling in cryo-electron microscopy. Journal of Structural Biology, 183(1), 19–32. http://doi.org/10.1016/j.jsb.2013.05.008
  * detailed motivation of an image formation (forward model). Authors include employees of the FEI company, which makes electron microscopes.
* Grob, P., Bean, D., Typke, D., Li, X., Nogales, E., & Glaeser, R. M. (2013). Ranking TEM cameras by their response to electron shot noise. Ultramicroscopy, 133, 1–7. http://doi.org/10.1016/j.ultramic.2013.01.003
  * Details mathematical analysis of detector noise, based on minimal assumptions (poisson stats of the number of general random variable events for the count registered, with mean and std given). See supplementary info for detailed proofs.
* Jensen, G. J. (2014). Getting Started in Cryo-EM. Retrieved from http://cryo-em-course.caltech.edu/videos
* Vulović, M., Voortman, L. M., Van Vliet, L. J., & Rieger, B. (2014). When to use the projection assumption and the weak-phase object approximation in phase contrast cryo-EM. Ultramicroscopy, 136, 61–66. http://doi.org/10.1016/j.ultramic.2013.08.002
* Öktem, O. (2015). Mathematics of Electron Tomography. In Handbook of Mathematical Methods in Imaging (pp. 937–1031). New York, NY: Springer New York. http://doi.org/10.1007/978-1-4939-0790-8_43
* Afanasyev, P., Ravelli, R. B. G., Matadeen, R., De Carlo, S., Van Duinen, G., Alewijnse, B., … Van Heel, M. (2015). A posteriori correction of camera characteristics from large image data sets. Scientific Reports, 5, 1–9. http://doi.org/10.1038/srep10317
  * Correction of detector defects (dust, electronics) from large batch of images. Post-factum gain correction. `camera_norm` code in IMAGIC. Discussion of what causes these artefacts (dust, etc), equations of how to estimate the parameters to correct, and examples from read datasets of 1000s of micrographs. 
* Kirkland, E. J. (2016). Computation in electron microscopy. Acta Crystallographica Section A Foundations and Advances, 72(1), 1–27. http://doi.org/10.1107/S205327331501757X
* McMullan, G., Faruqi, A. R., & Henderson, R. (2016). Direct Electron Detectors. In Methods in Enzymology (1st ed., Vol. 579, pp. 1–17). Elsevier Inc. http://doi.org/10.1016/bs.mie.2016.05.056
  * Excellent historical overview of improvements in electron detectors, especially the 1990s and 2000s. Good place to start for further delving into detector literature, especially in the journals *Ultramicroscopy* and *Nuclear Instruments & Methods in Physics Research Section A—Accelerators, Spectrometers, Detectors and Associated Equipment*.
* Maigné, A., & Wolf, M. (2018). Low-dose electron energy-loss spectroscopy using electron counting direct detectors. Microscopy, 67(suppl_1), i86–i97. http://doi.org/10.1093/jmicro/dfx088
  * Characterization of a Gatan K2 for EELS on a biological protein embedded in amorphous ice. Promising applications for time resolved spectra measurments on cryoEM samples.
* Booth, C. (2019). Detection Technologies for Cryo-Electron Microscopy.
* Datta, A., Ban, Y., Ding, M., Chee, S. W., Shi, J., & Loh, N. D. (2019). ReCoDe: A Data Reduction and Compression Description for High Throughput Time-Resolved Electron Microscopy.
* Kirkland, E. J. (2020). Advanced Computing in Electron Microscopy. Cham: Springer International Publishing. http://doi.org/10.1007/978-3-030-33260-0
  * Superb textbook! Very pedagical treatment of physics and computation. Emphasis on TEM and SEM as well as materials science and biological cryoEM. Computer programs are available online accompanying the text. Derivation and analysis of: wave equation for fast electrons, multislice (thick specimens). Important details about appropriate discretizations (rules of thumb and the reason why) and other matters surrounding numerical computations.
* Nakane, T., Kotecha, A., Sente, A., McMullan, G., Masiulis, S., Brown, P. M. G. E., … Scheres, S. H. W. (2020). Single-particle cryo-EM at atomic resolution. BioRxiv, 2020.05.22.110189. http://doi.org/10.1101/2020.05.22.110189
  * electron-event representation data format for electron counting
* Saxton, W. O. (2020). Advances in Imaging and Electron Physics: Computer Techniques for Image Processing in Electron Microscopy. (M. Hÿtch & P. W. Hawkes, Eds.).
  * Multichapter reprints from Advances in Electronics and Electron Physics, Supplement 10, 1978. Although much of the numerical implementation details are dated, the modelling choices and their motivations can be discerned.
* Sigworth, F., & Tagare, H. (2020). Cryo-EM Principles. Retrieved from https://cryoemprinciples.yale.edu/
* Himes, B. A., & Grigorieff, N. (2021). Cryo-TEM simulations of amorphous radiation-sensitive samples using multislice wave propagation. BioRxiv, 6. http://doi.org/10.1101/2021.02.19.431636
  * advanced forward model. Theory section with abundant references to published literature, especially the physics literature. Focuses on an the "accurate representation of molecular density", "compensating for the isolated atom superposition approximation", "modeling radiation damage", modeling the solvent envelope", the "accurate representation of solvent noise" and "amplitude contrast". Results show simulated micrographs and show how they follow physical trends that are more close to experiments, and also on the same quantitiative scale. Currently in C++ in cisTEM. Future work includes speeding up wtih GPU impelementation.
  
## CryoEM algorithms / data processing

### Overview of algorithms
* Cong, Y., & Ludtke, S. J. (2010). Single Particle Analysis at High Resolution. In Methods in Enzymology (1st ed., Vol. 482, pp. 211–235). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82009-9
* Scheres, S. H. W. (2010). Classification of Structural Heterogeneity by Maximum-Likelihood Methods. In Methods in Enzymology (1st ed., Vol. 482, pp. 295–320). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82012-9
* Jensen, G. J. (2014). Getting Started in Cryo-EM. Retrieved from http://cryo-em-course.caltech.edu/videos
* Sigworth, F. J. (2016). Principles of cryo-EM single-particle image processing. Microscopy (Oxford, England), 65(1), 57–67. http://doi.org/10.1093/jmicro/dfv370
* Singer, A. (2018). Mathematics for cryo-electron microscopy.
* Sigworth, F., & Tagare, H. (2020). Cryo-EM Principles. Retrieved from https://cryoemprinciples.yale.edu/
* Singer, A., & Sigworth, F. J. (2020). Computational Methods for Single-Particle Cryo-EM, 1–40.
* Bendory, T., Bartesaghi, A., & Singer, A. (2020). Single-Particle Cryo-Electron Microscopy: Mathematical Theory, Computational Challenges, and Opportunities. IEEE Signal Processing Magazine, 37(2), 58–76. http://doi.org/10.1109/MSP.2019.2957822

### Expectation-Maximization
* Sigworth, F. J. (1998). A Maximum-Likelihood Approach to Single-Particle Image Refinement. Journal of Structural Biology, 122(3), 328–339. http://doi.org/10.1006/jsbi.1998.4014
  * First application of expecation-maximization in cryoEM, although different language (maximum-liklihood) used by the author. Synthetic data analyzed in successfull attempt to overcome template bias in 2D class averages.
* Do, C. B., & Batzoglou, S. (2008). What is the expectation maximization algorithm? Nature Biotechnology, 26(8), 897–899. http://doi.org/10.1038/nbt1406
* Sigworth, F. J., Doerschuk, P. C., Carazo, J.-M., & Scheres, S. H. W. (2010). An Introduction to Maximum-Likelihood Methods in Cryo-EM. In Methods in Enzymology (1st ed., Vol. 482, pp. 263–294). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82011-7
* Tagare, H. D., Barthel, A., & Sigworth, F. J. (2010). An adaptive Expectation–Maximization algorithm with GPU implementation for electron cryomicroscopy. Journal of Structural Biology, 171(3), 256–265. http://doi.org/10.1016/j.jsb.2010.06.004
* Scheres, S. H. W. (2010). Classification of Structural Heterogeneity by Maximum-Likelihood Methods. In Methods in Enzymology (1st ed., Vol. 482, pp. 295–320). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82012-9
* Nelson, P. C. (2019). Chapter 12 : Single Particle Reconstruction in Cryo-electron Microscopy. In Physical Models of Living Systems (pp. 305–325).
  * Superb! Very good place to start. Very pedagogical treatment and motivation of the EM algorithm, in real space under Gaussian stats. 1D case (shifts) discussed first to motivate 2D case (shits and rotation). Care taken to explain what the motation means (subscripts, what is a vector, what is a scalar, etc.). Notebooks used to make the textbook are available on this repo by permission from the authors. 
  

### 2D Classification
* Rao, R., Moscovich, A., & Singer, A. (2020). Wasserstein K-Means for Clustering Tomographic Projections, (2016), 1–11.
  * The authors use a displacement/transport based loss, in place of the typical pixel based loss (e.g. L2). They propose an rotationally invariant k-means alignment algorithm. A fast approximation to the Earthmover's distance is employed, which is based on a weighted L1 distance between wavelet transforms. Compared with the Gaussian white noise case, where the L2 loss can be converted to a normalized probability, there is no such analytically convenient probabilitic transformation/interpretation with the Wasserstein loss. The appendix proves a useful result that provides a bound on how the L2 loss can be well behaved for smooth signals.

### Resolution (FSC, SSNR, etc)
* Harauz G, van Heel M. (1986) Exact filters for general geometry three dimensional reconstruction. Optik (Stuttg) 4: 146-156.
* Penczek, P. A. (2010). Resolution Measures in Molecular Electron Microscopy. In Methods in Enzymology (1st ed., Vol. 482, pp. 73–100). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82003-8
* Kucukelbir, A., Sigworth, F. J., & Tagare, H. D. (2014). Quantifying the local resolution of cryo-EM density maps. Nature Methods, 11(1), 63–65. http://doi.org/10.1038/nmeth.2727
* Penczek, P. A. (2020). Reliable cryo-EM resolution estimation with modified Fourier shell correlation. IUCrJ, 7(6), 1–14. http://doi.org/10.1107/s2052252520011574
  * modified FSC with real space masking done at end. Fourier shell still taken of each half map, but then goes back to real space and applies mask before correlation. Have to do an inverse FT for each shell, but can be parallelized.
* van Heel, M., & Schatz, M. (2020). Information: to Harvest, to Have and to Hold, 1–43.

### B-factor
* Rosenthal, P. B., & Henderson, R. (2003). Optimal Determination of Particle Orientation, Absolute Hand, and Contrast Loss in Single-particle Electron Cryomicroscopy. Journal of Molecular Biology, 333(4), 721–745. http://doi.org/10.1016/j.jmb.2003.07.013
  * Classic paper that studies "amplitude spectrum" B-factor, which can be computed from a map: "the loss of contrast in the map at high resolution [depends on ...] causes such as radiation damage, imaging imperfections, and errors in the reconstruction procedure." Includes good "Theoretical Background" section on "Guinier analysis", "Loss of contrast", and "Contrast restoration".
 * Singer, A. (2021). Wilson Statistics : Derivation , Generalization , and Applications to Electron Cryomicroscopy. BioRxiv, 1–16.
   * "the first rigorous mathematical derivation of Wilson statistics", which is the high frequency regime of the power spectra (Guinier is the low frequency)

### Ewald Sphere
* Wolf M, DeRosier DJ, Grigorieff N. (2006) Ewald sphere correction for single-particle electron microscopy. Ultramicroscopy 4-5: 376-82.
* Leong, P. A., Yu, X., Zhou, Z. H., & Jensen, G. J. (2010). Correcting for the Ewald Sphere in High-Resolution Single-Particle Reconstructions. In Methods in Enzymology (1st ed., Vol. 482, pp. 369–380). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82015-4
* Russo CJ, Henderson R. (2018) Ewald sphere correction using a single side-band image processing algorithm. Ultramicroscopy: 26-33.

### PCA
* Sorzano, C. O. S., & Carazo, J. M. (2021). Principal component analysis is limited to low-resolution analysis in cryoEM. Acta Crystallographica Section D Structural Biology, 77(6), 835–839. http://doi.org/10.1107/s2059798321002291
  * 4 page paper that goes over some of the math and gives advice on interpretability (reconstructured volume vs components along eigenvolumes, or eivenvolumes alone). Connections mentioned in passing with other literature on 3D deformations: Taylor expansions of the volume, Laplacian analysis of a graph, generalized Fourier transforms to arbitrary geometrical shapes (spherical harmonics for surface of sphere, Bessel functions for cylindar, prolate spheroidals for sphere). Simple pedagogical 2D example of moving line (metranome like).
* Punjani, A., & Fleet, D. J. (2021). 3D variability analysis: Resolving continuous flexibility and discrete heterogeneity from single particle cryo-EM. Journal of Structural Biology, 213(2), 107702. http://doi.org/10.1016/j.jsb.2021.107702
  * examples include: solving multiple dimensions of heterogeneity, high resolution flexible motion of small proteins, multiple high resolution modes of bending in a sodium ion channel, symmetric and asymmetric flexible motion at high resolution, large flexible motions of large complexes, directly resolving discrete heterogeneity. In practice the software implementation is fast and thus practically allows for a high amount of components (8+ vs 2-3).
* Tagare, H. D., Kucukelbir, A., Sigworth, F. J., Wang, H., & Rao, M. (2015). Directly reconstructing principal components of heterogeneous particles from cryo-EM images. Journal of Structural Biology, 191(2), 245–262. http://doi.org/10.1016/j.jsb.2015.05.007
  * detailed theoretical/analytical analysis of the Fourier slice theorem for covaniances: how do we estimate how pixels in 3D are correlating with each other when all we haev is information on their 2D projections? Rich in mathematical details.

### Deep generative 3D reconstruction
* Punjani, A., & Fleet, D. (2021). 3D Flexible Refinement : Structure and Motion of Flexible Proteins from Cryo-EM. BioRxiv, 1–21. http://doi.org/10.1101/2021.04.22.440893
  * Auto-decoder approach using one 3D reference map (learned from data), modified by vector field (convection, also learned from data). Each particle has its own convection (latent variable learned from data). Nice discussion of design choices that inspired confidence and is useful for methods developers in this area.

#### (Variational) Autoencoders
* Doersch, C. (2016). Tutorial on Variational Autoencoders, 1–23.
* Bepler, T., Zhong, E. D., Kelley, K., Brignole, E., & Berger, B. (2019). Explicitly disentangling image content from translation and rotation with spatial-VAE, (NeurIPS 2019).
* Zhong, E. D., Bepler, T., Davis, J. H., & Berger, B. (2019). Reconstructing continuous distributions of 3D protein structure from cryo-EM images, 1–20.
* Miolane, N., Poitevin, F., Holmes, S., & Li, Y. T. (2019). Estimation of orientation and camera parameters from cryo-electron microscopy images with variational autoencoders and generative adversarial networks. ArXiv.
  * Ribosome (one simulated dataset in same 2D pose, three empirical datasets from different 2D classes) analyzed with a VAE-GAN architecture. In latent space, defocus and rotations are disentagled (in unsupervised way).
* Nashed, Y., Peck, A., Martel, J., Levy, A., Koo, B., Wetzstein, G., … Poitevin, F. (2022). Heterogeneous reconstruction of deformable atomic models in Cryo-EM.
  * The authors demonstrate that the anisotropic network model modes could capture the heterogeneity of adenylate kinase transitioning between open and closed states. A continuous trajectory was discreatly sampled in 50 states with a tool in a molecular viewer program to generate synthetic data. An autoencoder was used to estimate the up to 14 normal mode components, where the estimated values were used in a physics decoder. Other latent variables (rotation, CTF defocus, open conformation of atomic model) were provided and not inferred. The physics decoder represented the full atomic model with multiple Gaussians with atom-type specific parameters from established tabulated values. The elastic network model is computed on a subset of atoms, and then interpolated for the remaining atoms, avoiding the expensive diagonalization of the 3N × 3N Hessian, where N is the number of atoms.
  
### Reconstruction
* Ludike, S. J., & Wah Chiu. (2002). Image restoration in sets of noisy electron micrographs. In Proceedings IEEE International Symposium on Biomedical Imaging (pp. 745–748). IEEE. http://doi.org/10.1109/ISBI.2002.1029365
* Penczek, P. A. (2010). Fundamentals of Three-Dimensional Reconstruction from Projections. In Methods in Enzymology (1st ed., Vol. 482, pp. 1–33). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82001-4
* Penczek, P. A. (2010). Image restoration in cryo-electron microscopy. Methods in Enzymology (1st ed., Vol. 482). Elsevier Inc. http://doi.org/10.1016/S0076-6879(10)82002-6



## Tomography
* Natterer, F. (1986). VII: Mathematical Tools. In The Mathematics of Computerized Tomography (pp. 180–212).
  * Terse notation requires familiarity with domain and mathematical symbols.
* Block, W. (2004). Computed Tomography Notes, Part 1 Challenges with Projection X-ray Systems. Retrieved from https://www.medphysics.wisc.edu/~block/bme530lectures/ct1.pdf
* Technical University of Denmark. (2020). Introduction to advanced tomography. https://www.coursera.org/learn/cinemaxe
