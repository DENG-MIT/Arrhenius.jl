# About the package

Arrhenius.jl is a Julia package for reacting-flow and combustion modeling, including differentiable chemistry workflows.

## Publications

- W. Ji, X. Su, B. Pang, S. J. Cassady, A. M. Ferris, Y. Li, Z. Ren, R. Hanson, S. Deng, [“Arrhenius.jl: A Differentiable Combustion Simulation Package”](https://arxiv.org/abs/2107.06172), arXiv:2107.06172 [physics.chem-ph], 2021.
- W. Ji, S. Deng, [“Autonomous Discovery of Unknown Reaction Pathways from Data by Chemical Reaction Neural Network”](https://doi.org/10.1021/acs.jpca.0c09316), *The Journal of Physical Chemistry A* 125(4):1082–1092, 2021.
- W. Ji, J. Zanders, J.-W. Park, S. Deng, [“Data-Driven Approaches to Learn HyChem Models”](https://doi.org/10.1115/icef2021-67925), *ASME 2021 Internal Combustion Engine Division Fall Technical Conference* (ICEF2021), 2021.
- X. Su, W. Ji, L. Zhang, W. Wu, Z. Ren, S. Deng, [“Neural Differential Equations for Inverse Modeling in Model Combustors”](https://doi.org/10.1115/imece2021-69657), *ASME 2021 International Mechanical Engineering Congress and Exposition* (IMECE2021), Volume 11: Heat Transfer and Thermal Engineering, 2021.
- W. Ji, X. Su, B. Pang, Y. Li, Z. Ren, S. Deng, [“SGD-based optimization in modeling combustion kinetics: Case studies in tuning mechanistic and hybrid kinetic models”](https://doi.org/10.1016/j.fuel.2022.124560), *Fuel* 324:124560, 2022.
- W. Ji, F. Richter, M. J. Gollner, S. Deng, [“Autonomous kinetic modeling of biomass pyrolysis using chemical reaction neural networks”](https://doi.org/10.1016/j.combustflame.2022.111992), *Combustion and Flame* 240:111992, 2022.

```bibtex
@misc{ji2021arrhenius,
  title = {Arrhenius.jl: A Differentiable Combustion Simulation Package},
  author = {Ji, Weiqi and Su, Xingyu and Pang, Bin and Cassady, Sean Joseph and Ferris, Alison M. and Li, Yujuan and Ren, Zhuyin and Hanson, Ronald and Deng, Sili},
  year = {2021}, eprint = {2107.06172}, archivePrefix = {arXiv}, primaryClass = {physics.chem-ph},
  url = {https://arxiv.org/abs/2107.06172}
}
@article{ji2021crnn,
  title = {Autonomous Discovery of Unknown Reaction Pathways from Data by Chemical Reaction Neural Network},
  author = {Ji, Weiqi and Deng, Sili},
  journal = {The Journal of Physical Chemistry A}, volume = {125}, number = {4}, pages = {1082--1092}, year = {2021},
  publisher = {American Chemical Society}, doi = {10.1021/acs.jpca.0c09316}
}
@inproceedings{ji2021hychem,
  title = {Data-Driven Approaches to Learn HyChem Models},
  author = {Ji, Weiqi and Zanders, Julian and Park, Ji-Woong and Deng, Sili},
  booktitle = {ASME 2021 Internal Combustion Engine Division Fall Technical Conference}, series = {ICEF2021}, year = {2021},
  publisher = {American Society of Mechanical Engineers}, doi = {10.1115/icef2021-67925}
}
@inproceedings{su2021neuralde,
  title = {Neural Differential Equations for Inverse Modeling in Model Combustors},
  author = {Su, Xingyu and Ji, Weiqi and Zhang, Long and Wu, Wantong and Ren, Zhuyin and Deng, Sili},
  booktitle = {Volume 11: Heat Transfer and Thermal Engineering}, series = {IMECE2021}, year = {2021},
  publisher = {American Society of Mechanical Engineers}, doi = {10.1115/imece2021-69657}
}
@article{ji2022sgd,
  title = {SGD-based optimization in modeling combustion kinetics: Case studies in tuning mechanistic and hybrid kinetic models},
  author = {Ji, Weiqi and Su, Xingyu and Pang, Bin and Li, Yujuan and Ren, Zhuyin and Deng, Sili},
  journal = {Fuel}, volume = {324}, pages = {124560}, year = {2022}, publisher = {Elsevier}, doi = {10.1016/j.fuel.2022.124560}
}
@article{ji2022biomass,
  title = {Autonomous kinetic modeling of biomass pyrolysis using chemical reaction neural networks},
  author = {Ji, Weiqi and Richter, Franz and Gollner, Michael J. and Deng, Sili},
  journal = {Combustion and Flame}, volume = {240}, pages = {111992}, year = {2022}, publisher = {Elsevier}, doi = {10.1016/j.combustflame.2022.111992}
}
```

## Related projects

- [ArrheniusActiveSubspace](https://github.com/DENG-MIT/ArrheniusActiveSubspace): sensitivity analysis and uncertainty quantification.
- [Arrhenius_Flame_1D](https://github.com/DENG-MIT/Arrhenius_Flame_1D): one-dimensional flame sensitivity analysis.
- [CRNN_HyChem](https://github.com/DENG-MIT/CRNN_HyChem): neural-network approaches to learning kinetic mechanisms.
- [DeepReduction](https://github.com/DENG-MIT/DeepReduction): mechanism reduction workflows.

Related software includes [Cantera](https://cantera.org/) for mechanism processing and reference calculations and [ReactionMechanismSimulator.jl](https://github.com/ReactionMechanismGenerator/ReactionMechanismSimulator.jl) for chemical mechanism simulation.

For the package license, see the [MIT License](https://github.com/DENG-MIT/Arrhenius.jl/blob/main/LICENSE).
