<!-- passonneau_and_carpenter_2014_the_benefits_of_a_model_of_annotation.md -->
- Passonneau and Carpenter (2014) "The Benefits of a Model of Annotation"
    - *introduction* (pp. 311f.)
        - standard validation methods of human-annotated datasets (textual and others) rely on measurements of inter-annotator consistency (e.g., the chance-adjusted inter-annotator agreement metric)
        - problems: 
            "Such measures, however, merely report how often annotators agree, with no direct measure of corpus quality, nor of the quality of individual items" (p. 311)
            - "high chance-adjusted interannotator agreement is neither necessary nor sufficient to ensure high quality gold-standard labels" (p. 311)
        - authors advocate the use of probabilistic annotation models which may be fit to the labels observed by multiple annotators to obtain estimates of 'true' class membership of items (p. 311) as well as the uncertainty associated with these estimates (p. 312), and hence to produce more reliable data about the quality of an annotated corpus
    - *agreement metrics versus a model* (pp. 312-317)
        - assumption: "high-confidence ground truth label for each annotated instance is the ultimate goal of annotation" (p. 312)
        - claim: "With many annotators to compare, the value of gathering a label can be quantified using information gain and mutual information" (p. 312)
        - *pairwise and chance-adjusted agreement measures* (pp. 312-314)
            - current best practice: an iterative process for the creation of annotation standards:
                1. (re)formulate the annotation task
                2. write/revise coding scheme and textual guidelines provided to annotators (incl. training examples) 
                3. 2+ annotators work independently to annotate a sample of datasets
                4. measure the inter-annotator agreement on the data sample
                    - if to low, start a new iteration at (2.)
                    - otherwise proceed
                5. create gold standard dataset where each item is annotated by a single annotator to obtain 'ground truth' labels
            - problems with this procedure:
                - not clear how to determine if a sample is representative of the corpus
                - agreement metrics 
                    - do not differentiate annotators based on their varying quality (i.e., accuracy)
                    - do not provide information about the quality of the 'ground truth' labels
            - notation:
                - *items*:= $i \in 1,..., I$ with the item-set being $\mathcal{I}$
                - *annotators*:= $j \in 1,..., J$  with the annotator-set being $\mathcal{J}$
                - *label classes in categorical coding scheme*:= $k \in {1,..., K}$  with the label-set being $\mathcal{K}$
                - *observed labels/annotations*:= $y_{i,j} \in {1,..., K}$ recorded for item $i$ by annotator $j$  with the annotation-set being $\mathcal{Y}$
            - agreement metrics:
                - *agreement* (p. 313):
                    - pairwise agreement between two annotator $m,n$ is the proportion of items both labeled and which they labeled identically: $A_{m,n} = \frac{1}{I} \sum_{i=1}^{I} \mathcal{I}(y_{i,m} = y_{i,n})$
                        - the MLE of agreement in a binomial model
                    - averaging over ${J \choose 2}$ pairwise agreements yields the agreement metric: $A_{m,n} = \frac{1}{{J \choose 2}} \sum_{m=1}^{J-1} \sum_{n=m+1}^{J-1} A_{m,n}$
                - *chance-adjusted agreement* (p. 313): "the proportion of observed agreements that are above the proportion expected by chance"
                    - *chance agreement*:= $IA_{m,n} = \frac{A_{m,n}-C_{m,n}}{1-C_{m,n}}$
                        - Cohen's $\kappa$:= $C_{m,n}^{\mbox{Cohen}} = \sum_{k=1}^K \psi_{m,k} \times \psi_{n,k}$, where $\psi_{j,k} = \frac{1}{I} \sum_{i=1}^I \mathcal{I}(y_{i,j}=k)$ gives the proportion of the label *k* in annotator $j$'s annotations
                            - assumes that "each annotator draws uniformly at random from her set of labels" (p. 313) 
                        - Krippendorff's $\alpha$:= $C_{m,n}^{\mbox{Krippendorff}} = \sum_{k=1}^K \phi_k^2$, where $phi_k$ is the proportion of label $k$ in the entire set of annotations
                            - "assumes each annotator draws uniformly at random from the pooled set of labels from all annotators"
                    - takes into account the (empirical) prevalence of labels $k$ (annotators may agree by either or both guessing in accordance with the majority class label)
            - shortcomings of agreement metrics:
                1. intrinsically pairwise (though comparison to a voted consensus or average over multiple pairwise agreements possible)
                2. two wrongs make an agreement
                3. when chance agreement is high, even high-accuracy annotators can have low chance-adjusted agreement
                    - see numerical example on p. 313
                    - therefore, "high chance-adjusted interannotator agreement is not a necessary condition for a high-quality corpus" (p. 313)
                4. agreement statistics implicitly assume annotators are unbiased
                    - bias of two annotators in the same direction leads to overestimation of reliability
                    - averaging over pairwise agreement can thus induces systematic measurement bias, too
                5. item-level effects (e.g., item difficulty) can inflate agreement-in-error
                6. no confidence intervals provided
                    - uncertainty in annotator accuracy propagates into the error bounds of agreement metrics, so that they may "span the rather arbitrary decision boundaries for acceptability employed for interannotator agreement statistics" (p. 314)
        - *a probabilistic annotation model* (pp. 314-317)
            - allows to estimate corpus attributes (e.g. 'true' prevalence and item difficulty), as well as annotator attributes (e.g., bias/accuracy, sensitivity and specificity, etc.)
            - referring to Dawid and Skene's seminal paper (1979)
                - model to determine a consensus among patient histories taken by multiple doctors
                - inference driven by accuracies estimated on a per-annotator, per-category basis
                - consider incomplete-panel case (i.e., annotators $K$ label items $I$ $n_{i,j} \in 0,...,N$ times)
                    - data structure: a $N$-by-3 matrix, with
                        - items $ii[n] \in \mathcal{I}$ in the first column,
                        - annotators $jj[n] \in \mathcal{J}$ in the second column, and
                        - labels $y[n] \in \{y_{i,j}\}_ {i \in \mathcal{I}, j \in \mathcal{J}}$ in the third column
                    - estimates:
                        - *'true' category of item $i$*:= $Z_i \in \mathcal{K}$
                        - *'true' prevalence of category $k$*:= $\pi_k \in [0,1]$ with $\boldsymbol{\pi}$ being a $K$-vector
                        - *per-annotator per-category accuracy/error*:= $\theta_{j,y,k} \in [0,1]$ with $y$ being the label recorded by $j$ for an item with true category $k$
                            - quantifies accuracy for $y = k$ and error for $y \new k$, respectively
                            - [hli:] specifies a confusion matrix for each annotator
            - authors implement Dawid and Skene's model in a Bayesian framework, using Dirichlet priors to smooth per-annotator and per-item parameters, respectively:
                - $z_i \sim \mbox{Categorical}(\boldsymbol{\pi})$ 
                - $y_n \sim \mbox{Categorical}(\boldsymbol{\theta}_ {jj[n], z[ii[n]]})$, where the argument to the distribution is the $K$-vector of the category-specific accuracy/error of annotator $j = jj[n]$ given 'true' category $z[ii[n]]$
                - $\boldsymbol{\theta}_ {j,k} \sim \mbox{Dirichlet}(\boldsymbol{\alpha}_ k)$
                - $\boldsymbol{\pi} \sim \mbox{Dirichlet}(\boldsymbol{\beta})$
            - specifying hyperparameters $\boldsymbol{\alpha}_ k$ and $\boldsymbol{\beta}$ to be unit-vectors yields the unsmoothed MLE, which is identical to the Maximum A-Posteriori (MAP) estimate
            - *estimating 'true' categories* (pp. 315f.)
                - $\Pr(z_i|\mathbf{y}, \boldsymbol{\theta}, \boldsymbol{\pi}) \propto \Pr(z_i | \boldsymbol{\pi})\Pr(\mathbf{y} | z_i,  \boldsymbol{\theta}) = \pi_{z[i]} \prod_{ii[n] = i} \theta_{jj[n], z[i], y[n]}$
                - accounting for annotators' biases/accuracies may induce to assign a different category than implied by majority vote, see MWE on page 315
                - advantages of the model
                    - adjusts for spam annotators (i.e., annotators, whose annotations are independent from items' true categories)
                    - exploit the (limited) information provided by biased annotators (i.e., annotators, who are more accurate on some categories compared to other categories)
                    - exploits information provided by adversarial annotators (i.e., annotators, whose annotations are negatively correlated with 'true' categories) 
                        - problem when too many adversarial annotators present in the data
            - *estimating the information value of annotations* (p. 316) 
                - comparing the uncertainty associated with $z_i$ before and after adding a new annotation $y_{i,j}$ allows to *compute the reduction in uncertainty* provided by an additional annotation by annotator $j$ for any given item $i$:
                    - *entropy*:= \mbox{H}(Z_i) = - \sum_{k=1}^K \Pr_{Z_i}(k) log(\Pr_{Z_i}(k)) with 
                        - $Z_i$ being the random variable corresponding to item $i$'s true category $k$, and 
                        - $\Pr_{Z_i}$ being the probability mass function of $Z_i$
                    - *conditional entropy*:= $\mbox{H}(Z_i | Y_n = k')$, where $Y_n$ is the random variable corresponding to $k'$, the annotation of annotator $j$ on $i$
                    - *mutual information*:= $\mbox{I}(Z_i, Y_n) = \mbox{H}(Z_i) - \mbox{H}(Z_i | Y_n)$ gives the expected reduction in entropy in $Z_i$ after observing an annotation $Y_n$
                - information quality: 
                    - spam annotators: \mbox{H}(Z_i | Y_n) = \mbox{H}(Z_i)
                    - perfect annotator: \mbox{H}(Z_i | Y_n) = 0 (i.e., reduces all uncertainty)
            - *implementation and priors* (pp. 316f.)
                - Dirichlet priors specified for stabilization purposes only (i.e., for cases where an annotator did not label any instances of a category)
                - modeling per-item random effects would be more realistic, and would generally cause confidence intervals on per-annotator estimates to widen (p. 316, referring to Carpenter 2008)
    - *two data collections*
        - fit the model to the MASC Word Sense Sentence Corpus and a crowd-sourced corpus for the same word-sense annotation task (R code distributed via Carpenter's Github account)
        - note that
            - estimating tight confidence intervals for the prevalence of low-count categories requires a sufficient number of annotations per item (as estimated error goes down as $\mathcal{O}(1/\sqrt{n})$ with $n$ independent annotations, p. 319)
        - distributed tasks only to crowd-workers with 98% lifetime approval and min. 20,000 successfully completed tasks to filter out spam annotators
        - include data from annotators who were blocked throughout the exercise due to low-accuracy performance to demonstrate that the model-based estimation manages to account for the bad quality of their annotations
    - *discussion* (pp. 321f.)
        - stress that only 10% in the expert-labeled data has been annotated by more than one annotator, preventing from computing agreement levels for each item
        - stress that experimental research has demonstrated that model-based estimation of categories fares better than majority vote (referring to Snow et al., 2008)
        - argue that "[there] is no evidence that a label from a trained annotator provides more information than a Turker's" (p. 321)
    - *conclusion* (pp. 323f.)
        - critique in a nutshell: "When two or more annotators have very high interannotator agreement on a task, unless they have perfect accuracy, there will be instances where they agreed incorrectly, and no way to predict which instances these are."
        - stress utility of model-based estimates of items' categories for machine learning applications: "Those who would use the corpus for training benefit because they can differentiate high from low confidence labels"


