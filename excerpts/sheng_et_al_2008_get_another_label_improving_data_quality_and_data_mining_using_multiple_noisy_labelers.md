<!-- sheng_et_al_2008_get_another_label_improving_data_quality_and_data_mining_using_multiple_noisy_labelers.md -->
- Sheng et al. (2008) "Get Another Label? Improving Data Quality and Data Mining Using Multiple, Noisy Labelers"
    - address the problem of whether repeated acquisition of labels for data improves categorization of items when the labeling is imperfect (i.e., 'noisy')
    - show that repeated labeling
        (i) can improve label and model quality
        (ii) is preferable to single labeling when labels are noisy, even in the traditional setting where labels are not particularly cheap
        (iii) is advantageous whenever the cost of processing the unlabeled data is not free (i.e., the cost of acquisition of an additional instance vs. another label for an already labeled instance)
        (iv) can be deployed selectively to economize on resources by exploiting uncertainty to select data points for which quality should be improved by requesting additional labels by annotators that let expect the biggest possible gain in mutual information
    - [hli:] note that the authors have published an updated version of their article in 2014 as Ipeirotis (2014) "Repeated labeling using multiple noisy labelers"
    - *introduction* (pp. 614f.)
        - problem setting: 
            - obtain certain data values (“labels”) relatively cheaply, from multiple sources (“labelers”
            - labels may be noisy due to lack of expertise, dedication, interest, or other factors
        - *repeated labeling*:= obtaining multiple labels for some or all data points
        - stress that the quality of statistical learning models depends on *both* the quality of training labels and the number of labeled instance: hence a trade-off arises pitting repeated labeling against the labeling of a new data point (p. 615)
    - *related work* (p. 615)
        - refer to work by Dawid and Skene, Smyth an colleagues, on the problem of probabilistic learning 
    - *repeated labeling: basics* (pp. 615) 
        - *notation and assumptions* (p. 616) 
            - formally, they consider a special case of the problem of supervised induction of a (binary) classification model
                - training example $(y_i, x_i)$
                    - $x_i$ is the unlabeled feature of example $i$ 
                    - $y_i$ is the 'true' label of example $i$
                - labeling incurs cost $C_L$ 
                - labeling is error prone:
                    - each label $y_{ij}$ comes from an imperfect, noisy labeler $j$
                    - labelers *labeling quality* is measured as $p_j = \Pr(y_{ij} = y_i)$
            - assume that labeling quality is independent from the features of an example, i.e., $\Pr(y_{ij} = y_i) = \Pr(y_{ij} = y_i | x_i)$
        - *majority voting and label quality* (pp. 616)
            - *integrated label quality* (ILQ):= $q = \Pr(\hat{y} = y)$ is the quality of the integrated label of a labeled instance
            - the *integrated label* of an instance may be determined by different methods, e.g., voting
            - with common labeler quality (i.e., $p_j = p \forall j \in 1,...,J$) the ILQ is the sum of the probabilities of having more correct than incorrect labels  (equ. 1, p. 616)
                - generally, for $p > .5$, the ILQ increases as the number of labelers increases, but so does the marginal improvement (see Figure 2, p. 616)
            - with varying labeler quality, there exists a threshold in differences among labelers that determines whether it is better to aggregate labels or use the best labelers label (Figure 4, p. 617) 
        - *uncertainty-preserving labeling*
            - stress that majority voting (MJ) ignores information about labeler quality
            - propose  a *multiplied examples* (ME) approach to the cosntruction of training sets: for each multiset of labels on an item, represent the item with a certain label, and assign it a weight proportional to the number of times the instance has been annotated with this label
                - stress that many learning algorithms can incorporate case weights
    - *repeated-labeling and modeling* (pp. 617f.)
        - *experimental setup*
            - use 12 real-world datasets to examine how different noisy annotation models fare in recovering true labels
            - 70% of examples of each dataset are used to train a classifier, which is then evaluated on the hold-out set
            - experimentally manipulate labeler quality by inducing error in the labels of training examples at rate $1-p$
        - *generalized round-robin strategies* (pp. 617-9)
            - *generalized round-robin (GRR) strategy*:= acquisition of a fix number of $k$ labels for each example dataset
            - study the decision boundary between the options of (a) acquiring a new training example (total cost: $C_U + C_L$), and (b) getting another label for an existing example (cost: $C_L$)
            - case when $C_U << CL$:
                - demonstrate that repeated-labeling with MV on average outperforms single labeling (SL) if a minimum number of training datasets have been acquired, and if enough labels have been collected per training dataset (Figure 5(a), p. 618)
                - conclude that "repeated-labeling is a viable alternative to single-labeling, even when the cost of acquiring the 'feature' part of an example is negligible compared to the cost of label acquisition"
            - general costs when $C_U >> CL$:
                - total cost is $C_D = C_U T_   au + C_L N_L$
                - fix the cost ratio $\rho = C_U/C_L = 3$ (i.e., acquiring a new example dataset is three times as costly as repeated-labeling) 
                - compare the performance of SL, MV, and ME
                - demonstrate that for $p = .6$ 
                    - MV and ME dominate SL in virtually all datasets under examination (see, Figure 6, p. 619)
                    - averaged across $C_D$, ME performs always as least as good as MV, but tends to outperform MV at higher levels of $C_D$ (see, Figure 8, p. 619)
                - conclude that "when labeler quality is low, inductive modeling often can benefit from the explicit representation of the uncertainty incorporated in the multiset of labels for each example."
        - *selective repeated-labeling* (pp. 619)
            - explain why the intuitively straightforward procedure of selecting instances for repeated-labeling based on resolving impurity in the multiset of labels is misleading: "under noise [entropy and distance of the majority label to the decision threshold] do not really measure the uncertainty in the estimation of the class label"
                - e.g., the multiset $\{+, +, +\}$ is pure, but with p = .6 it is not 95% certain that the true label is $+$ (p. 619)
                - the round-robin strategy therefore outperforms entropy-based selective repeated-labeling strategies in the long-run, because: "with a high noise level, the long-run label mixture will be quite impure, even though the true class of the example may be quite certain ... [while] [m]ore-pure, but incorrect, label multisets are never revisited" (p. 619)
            - *estimating label uncertainty* (LU, pp. 619f.)
                - provide a Bayesian estimate of the uncertainty in the class of an example
                    - specify a uniform prior on the distribution of the true label quality $p(y)$
                    - the posterior is thus $\mbox{Beta}(1 + L_{pos}, 1 + L_{neg})$ 
                - equate the uncertainty to the *posterior error probability* (i.e., the probability mass below the labeling decision threshold)
                    - take the CDF of the regularized incomplete beta function $I_x(\alpha, \beta)$
                - the selection rule is thus defined as S_{LU} = \min\{I_{.5}(L_{pos}, L_{neg}), 1 - I_{.5}(L_{pos}, L_{neg})\}
                - compare the performance of $S_{LU}$ to the GRR strategy, and find that the former tends to outperform the latter (i.e., it systematically achieves higher labeling quality as the number of labels increases) on all datasets under examination
            - *using model uncertainty* (MU, p. 620)
                - apply active learning scoring to select instances for which additional labels will to be acquired:
                    $S_{MU} = .5 - \left|\frac{1}{m} \sum_{i=1}^m \Pr(+|x, H_i) - .5 \right|$, with $\Pr(+|x, H_i)$ defined as the probability that an instance with features $x$ is classified as $+$ by model $H_i$, and $m$ is the number of learned models
                - use $m = 10$ and $H$ = random forests 
            - *using label and model uncertainty* (LMU, p. 620)
                - define a combined uncertainty metric as $S_{LMU} = \sqrt{S_{LU} \times S_{MU}}
                - demonstrate that LMU strongly dominates all other selection strategies in terms of improvements in label quality as a function of the number of labels acquired (Figure 10, p. 620)
            - *model performance with selective ML* (pp. 620)
                - train models using different four different selection rules: GRR, LU, MU, and LMU
                - examine accuracy on hold-out sample as a function of the number of labels used for training
                - find taht LU and LMU are consistently better than GRR, and that combining label and model uncertainty is the best approach (Figure 11, p. 621)
    - *conclusions, limitations, and future work* (pp. 621f.)
        - "Repeated-labeling is a tool that should be considered whenever labeling might be noisy, but can be repeated" (p. 621)
        - "when quality is low, preserving the uncertainty in the label multisets
        for learning [25] can give considerable added value."
        - limitations:
            - do not discuss how to incorporate labelers varying quality in their selection procedure, e.g., for deciding which labeler should be tasked with repeated-labeling
            - do not address correlation in labelers errors (e.g., due to different item difficulty or systematic biases)
            - do not discuss how their LMU selection strategy could be adapted to incorporate information about item difficulty, e.g., in order to deliberately allocate repeated-labeling resources to most difficult items


