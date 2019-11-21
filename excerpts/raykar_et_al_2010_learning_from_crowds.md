<!-- raykar_et_al_2010_learning_from_crowds.md -->
- Raykar et al. (2010) "Learning From Crowds"
    - *introduction* (pp. 1298-1300)
        - consider the context of repeated annotation of data by multiple, error-prone annotators, and the questions 
            (i) how to adapt conventional supervised learning algorithms to deal with noisy labels absent a gold standard?
            (ii) how to evaluate systems absent a gold standard?
            (iii) how to estimate how reliable different annotators are?
        - stress that the problem with majority voting (MV) is that it assumes that all annotators are equally reliable, but a majority of unreliable novices will always dominate experts (though novices may make unsystematic errors)
        - one faces a chivek-and-egg problem if the gold standard does not exist, hence the use of the EM algorithm, which treats item labels as missing data
        - make *novel contributions*:
            1. specifically address the problem of learning a classifier instead of just addressing the problem of estimating the ground truth
            2. the proposed algorithm estimates and learns the ground truth jointly
            3. method can be used in conjunction with any supervised learning algorithm
            4. model is intrinsically Bayesian (refer to Carpenter, 2008, for a full exposition of estimating Bayesian annotation models)
    - *binary classification* (pp. 1300-)
        - annotator performance measured as sensitivity and specificity wrt to the *unknown* gold standard
        - *a two coin analogy*
            - $j$ assigns label y_j \in \{0,1\} to the instance $\mathbf{x}$
            - $y$ is the (unobserved) true label of the instance
            - if $y = 1$, $j$ flips a coin with bias \alpha_j (sensitivity) and records $y_j = 0$ with probability $1 - \alpha_j$ (false negative)
            - if $y = 0$, $j$ flips a coin with bias $\beta_j$ (specificity) and records $y_j = 1$ with probability $1 - \beta_j$ (false positive)
        - it is assumed that $alpha_j$ and $\beta_j$ do not depend on the data (i.e., there are no per-instance effects on labeling performance)
        - *classification model*
            - consider the family of linear discriminating functions $\mathcal{F} = \{f_w\}$ where
             - $f_w(\mathbf{x}) = \mathbf{w}^\mbox{T}\mathbf{x})$, and 
             - $\Pr(y = 1 | \mathbf{w}, \mathbf{x}) = \sigma (\mathbf{w}^\mbox{T}\mathbf{x})$
            - in particular, logistic regression has \sigma (z) = 1/(1+e^{-z}) 
        - *learning problem*
            - training data:= $\mathcal{D} = \{x_i,y_{i,1},\ldots , y_{i,R}\}_ {i=1}^N$, with $R$, the number of annotators, and $N$, the number of training datasets  
            - the task is to estimate $\mathbf{w}$, $\boldsymbol{\alpha}$, and $\boldsymbol{\beta}$
        - MLE
            - assume that the training data is independently sampled 
            - the MLE of $\boldsymbol{\theta} = \{\mathbf{w}, \boldsymbol{\alpha}, \boldsymbol{\beta}\}$ maximizes $\Pr(\mathcal{D} | \boldsymbol{\theta}) = \prod_{i=1}^N \Pr(y_{i,1}, \ldots, y_{i,R} | \mathbf{x}_ i, \boldsymbol{\theta})$
            - $\boldsymbol{\theta}_ {ML} = \mbox{argmax}\{\ln \Pr(\mathcal{D} | \boldsymbol{\theta})\}$
        - the EM algorithm
            - refer to the EMA as "an efficient iterative procedure to compute the maximum-likelihood solution in presence of missing/hidden data"
            - consider true labels $\mathbf{y}$ as missing data 
            - if $\mathbf{y}$ were known, $\ln \Pr(\mathcal{D}, \mathbf{y} | \boldsymbol{\theta})$ could be specified (see p. 1302)
            - the EMA iteratively alternates between two steps:
                - expectation step
                    - it computes $\matchbb{E}\left[\ln \Pr(\mathcal{D}, \mathbf{y} | \boldsymbol{\theta})\right]$ with $\mu_i$ instead of $y_i$, where 
                        - $\mu_i = \frac{a_i p_i}{a_i p_i + b_i(1-p_i)}$,
                        - $p_i = \sigma (\mathbf{w}^\mbox{T}\mathbf{x})$
                        - $a_i = \prod_{j=1}^R \alpha_j^{y_{i,j}} (1-\alpha_j)^{1-y_{i,j}}$ (i.e., the product of annotators sensitivities wrt to their annotations)
                        - $b_i = \prod_{j=1}^R \beta^{1-y_{i,j}} (1-\beta_j)^{y_{i,j}}$ (i.e., the product of annotators specificities wrt to their annotations)
                - maximization step
                    - based on the current value of $\mu_i$ and the observations in $\mathcal{D}$, parameters are estimated by maximizing the conditional expectation:
                        - $\alpha_i = \frac{\sum_{i=1}^N \mu_i y_{i,j}}{\sum_{i=1}^N \mu_i}$
                        - $\beta_i = \frac{\sum_{i=1}^N (1 - \mu_i) (1 - y_{i,j})}{\sum_{i=1}^N (1 - \mu_i)}$
                        - $\mathbf{w}$ is obtained using Newtom-Raphson optimization, for no closed form exists due to the non-linearity of the sigmoid
            - propose to use majority voting as start values for $\mu_i$ in practice
            

