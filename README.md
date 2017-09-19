# Structural-Kernel-Learning-for-HDBBO
This is the MATLAB code repository associated with the paper [_Batched High-dimensional Bayesian Optimization via Structural Kernel Learning_](https://arxiv.org/pdf/1703.01973.pdf).

## System Requirement
We tested our code with MATLAB R2015b on Ubuntu 14.04 LTS (64-bit). Our code should work out of box without installing additional packages. However, please make sure you installed the GNU Scientific Library ([GSL](http://www.gnu.org/software/gsl/)). On Ubuntu, you can install GSL by 

```
sudo apt-get install libgsl0-dev
```

In MATLAB command line, you can mex the c files, for example, by

```
mex chol2invchol.c -lgsl -lblas
```

To run Box2D related code in test-functions/robot-pushing/, please install [Pybox2d](https://github.com/pybox2d/pybox2d).

## Example
example_addgp.m is a simple example using Bayesian Optimization (BO) to maximize a high-dimensional black-box function. Please see the comments in the code for more details about the usage.

add_gpopt.m is the function for BO with additive Gaussian processes for high dimensional problems.

example_batch_addgp.m is an example using *Batch* Bayesian Optimization (BBO) to maximize a high-dimensional black-box function. Please see the comments in the code for more details about the usage. Different methods could be used by choosing the value of `bo_method` from 

* `batch_rand`: Each point for exploration is chosen uniformly randomly
* `batch_ucb_dpp_ori`: All acquisition functions are UCB. Exploration is done via DPP with posterior covariance kernels for each group. Combination is done uniformly randomly.
* `batch_ucb_pe_ori`: All acquisition functions are UCB. Exploration is done via PE with posterior covariance kernels for each group. Combination is done uniformly randomly.
* `batch_ucb_dpp_fnc`: All acquisition functions are UCB. Exploration is done via DPP with posterior covariance kernels for each group. Combination is done by maximizing the quality function, which is UCB in our case.
* `batch_ucb_pe_fnc`: All acquisition functions are UCB. Exploration is done via PE with posterior covariance kernels for each group. Combination is done by maximizing the quality function, which is UCB in our case.

batch_add_gpopt.m is the function for BBO with additive Gaussian processes for high dimensional problems.

## Example functionals for optimization
In test_functions/, we provide some functionals one can use to test Bayesian optimization algorithms. 

1. Approximated functions sampled from Gaussian processes: sample_GPprior.m, sample_addGP.m.
2. Optimization test functions: egg.m, michNf.m, shekelMf.m.
3. Tuning hyper-parameters for neural networks: nnet_boston.m, nnet_cancer.m.
4. Active learning for robot pushing: robot_pushing_3.m, robot_pushing_4.m, two_robot_pushing.m.
5. Tuning the walking speed of a planar bipedal robot: walker_speed.m (Westervelt el al., 2007).

## Citation
Please cite our work if you would like to use the code.
```
@inproceedings{wang2017batched,
  title={Batched High-dimensional Bayesian Optimization via Structural Kernel Learning},
  author={Wang, Zi and Li, Chengtao and Jegelka, Stefanie and Kohli, Pushmeet},
  booktitle={International Conference on Machine Learning (ICML)},
  year={2017}
}
```
