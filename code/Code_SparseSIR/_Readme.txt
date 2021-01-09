Description for the codes:

1) The m-file a_plot_run.m can  produce Fig 1 – Fig 4 (in the paper), and Fig 1 – Fig 3 (in the supplementary material). 
2) The m-file a_tab_run.m can produce Tab 1 – Tab 2 (in the paper), and Tab 1 – Tab 3 (in the supplementary material). 
3) The m-file ssir_natural.m produces the natural estimator by solving Eqn (8) through ADMM algorithm.
4) The m-file ssir_refine.m produces the refined estimator by solving Eqn (20) through group lasso algorithm.


Steps for running the codes:

Step 1. Install the Matlab package TFOCS.

The attached file includes a file named TFOCS-master, which is downloaded from http://cvxr.com/tfocs/download/. This is a popular tool in convex optimization and is critical when solving Eqn (8) in the submitted paper. It can be easily installed into Matlab by the process “Set Path -> Add with subfolders -> choose the package -> Save”.

Step 2. Install the Matlab package glmnet.

The attached file has another package called glmnet, which is from https://web.stanford.edu/~hastie/glmnet_matlab/download.html. It’s a well-known package for implementing Lasso and regularized generalized linear regression. In order to solve Eqn (20) in the submitted paper, we changed the line 152 of glmnetSet.m (in the package of glmnet) to be “options.intr = false”, the reason for this change is that there is no intercept term when we perform group lasso for Eqn (20). It can also be easily installed into Matlab by the process “Set Path -> Add with subfolders -> choose the package -> Save”.

Step 3. Run a_plot_run.m file to obtain figures like Fig 1 – Fig 4.

Step 4. Run a_tab_run.m file to obtain Tables like Table 1 – Table 2.



Note: The results in the submitted papers are based on k = 100 repetitions, and the running time for producing Fig 1 and Table 1  are approximately 1 day and  3 days, respectively. For your convenience, the codes in the attached file are all set to run k = 2 times, and you can easily increase the number of repeating times to get more accurate results. The screenshots of running the code based on k = 2 are attached for ease of reference. 