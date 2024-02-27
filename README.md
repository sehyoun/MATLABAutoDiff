# Automatic Differentiation for MATLAB

This toolbox implements automatic/algorithmic differentiation for matlab using sparse representation for jacobians.

For example usage/syntax, check
* [Syntax Example](http://sehyoun.com/EXAMPLE_AutoDiff_Syntax.html)
* [Speed test against symbolic differentiation](http://sehyoun.com/EXAMPLE_test_symbolic.html)

It is now (Feb, 2024) possible to automatically differentiation of Schur/QZ/SVD matrix decomposition. They are only unique up to rotations. Refer to the appendix G of this [paper](https://www.norges-bank.no/en/news-events/news-publications/Papers/Working-Papers/2023/wp-11-2023/) for details for the Schur decomposition. As the expressions are only unique up to rotations, interface with automatic differentiation requires a different syntax than the "just use the same code and call myAD." Refer to the example files for the syntax.

For more detailed explanation and list of supported files, read the [documentation](https://github.com/sehyoun/MATLABAutoDiff/blob/master/README.pdf).
