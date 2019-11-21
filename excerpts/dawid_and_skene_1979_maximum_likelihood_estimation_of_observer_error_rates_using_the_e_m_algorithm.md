<!-- dawid_and_skene_1979_maximum_likelihood_estimation_of_observer_error_rates_using_the_e_m_algorithm.md -->
- Dawid and Skene (1979) "Maximum Likelihood Estimation of Observer Error-Rates Using the EM Algorithm"
    - consider the case where multiple patients $1,...,I$ are classified by multiple doctors/clinicians/observers $1,...,K$ into multiple classes $1,...,J$
        - [hli:] this is the (in)complete panel of multiple noisy annotators labeling multiple items  
    - *introduction*
        - assume that "many facets [of this process] are subject to errors of measurements" (p. 20) 
            - *individual error-rates* := $\pi_{j,l}^{(k)}$ is the probability that observer *k* records value *l* given that *j* is the true response/category 
                - includes case where $j = l$, where *k* records the correct value
                - assumed to be independent from the marginal distribution of the true response (p. 21)
        - enumerate different reasons why learning individual error-rates may be interesting (p. 21), especially that individual error-rates can be used to obtain a weighted majority vote on a patients class
        - when the true response $l$ is known, a sensible estimate of individual error-rates can be calculated as $\hat{\pi}^{(k)}$ can be computed for each combination of $j$ and $l$ by dividing 
            (i) the number of times observer $k$ records value l when j is correct by
            (ii) the number of patients seen by observer $k$ where j is correct
        - problem arises when estimating individual error-rates and true responses are not known (!!!)
    - *MLE*
        - *number of times observer $k$ gets response $j$ from patient $i$* := $\n_{i,l}^{(k)}$
        - assumptions: 
            1. $L^{(i)}$, the values recorded for a single patient by successive clinicians and/or to repeated questioning by a single clinician, are independent given the true response
            2. patients respond independently, i.e., $J = (J^{(1)}, ..., L^{(i)})$ are independent
            3 no patient-by-clinician interaction (equivalently, patients do not differ in 'difficulty', i.e., no clinician obtains a more helpful response from a patient than any other clinician)
        - **case 1: true responses are available**
            - *true responses for patient $i$*:= $\{T_{ij}: j = 1,...,J\}$
                - if $i$'s true response is $q$, then $T_{iq}=1$ and $T_{ij}=0\ \forall j \new q$
            - *true prevalence of $j$* (typically unknown): $p_j = \Pr(j=1)$
            - consider a patient-clinician tuple $(k,i)$:
                - if $i$'s true response were $q$ (i.e., $T_{iq}=1$), values recorded by $k$ would be distributed multinomial: the likelihood that $\Pr(\pi^{(k)}|n_{i,l},j,l=q) \propto \prod_{l=1}^J\left(\pi_{q,l}^{(k)}\right)^{n_{i,l}^{(k)}}$
                    - i.e., the probability that $k$ records $l=1$ is $\pi_{q,1}^{(k)}$ (her individual error-rate for detecting responses in category $q$), that she records $l=2$ is $\pi_{q,2}^{(k)}$, etc.
                        - [hli:] if $k$ is very sensitive to detecting response $q$ (i.e., $\pi_{q,l=q}^{(k)}$ is close to one), the probability mass of her multinomial will be centered on $q$
                    - [hli:] suppose $j \in \{1,2\}$
                        - if $j = 1$, $\Pr(l=j=1) = \pi_{1,1}^{(k)}$ and $\Pr(l=2) = 1 - \Pr(l=j=1)$
                        - it is as if $k$ were throwing a biased coin and recoding a value that matches the true response only with $\pi_{q,l=q}^{(k)}$
                        - it is as if $k$ has a biased coin for each true response (she may have higher sensitivity than specificity, i.e., may be better at detecting true-positives than at detecting true-negatives)
                - as all clinicians record values independently $\Pr(\boldsymbol{\pi}|j,l=q) \propto \prod_{k=1}^K\prod_{l=1}^J\left(\pi_{q,l}^{(k)}\right)^{n_{i,l}^{(k)}}$
                - *equ 2.1*: the probability of all the data on patient $i$ (true response $T_{ij}$ and recorded values L^{(i)}):
                $\Pr(T_{ij}, L^{(i)}|p_j,\boldsymbol{\pi},\mathbf{n}_{i,l}) \propto \prod_{j=1}^J\left\{p_j\prod_{k=1}^K\prod_{l=1}^J\left(\pi_{q,l}^{(k)}\right)^{n_{i,l}^{(k)}} \right\}^{T_{ij}}$
                    - note that the exponent, $T_{ij} \in \{0,1\}$ depending on $i$'s true response (here: $q$), and hence, only the product factors where $j = q$, the other $J-1$ products equal one
                    - for $j = q$ the product equals \Pr(L^{(i)} | T_{iq}=1)\Pr(T_{iq}=1)
                - *equ 2.2*: hence the likelihood function for the full data is $\prod_{i=1}^I\prod_{j=1}^J\left\{p_j\prod_{k=1}^K\prod_{l=1}^J\left(\pi_{q,l}^{(k)} \right)^{n_{i,l}^{(k)}} \right\}^{T_{ij}}$
                - the closed form solutions for the MLEs are
                    - *equ 2.3*: $\hat{\pi}_ {j,l}^{(k)} = \frac{\sum_i T_{ij}n_{il}^{(k)}}{\sum_l\sum_i T_{ij}n_{il}^{(k)}}$
                    - *equ 2.4*: $\hat{p}_j\sum_i T_{ij} / I$
        - **case 2: true responses are not available**
            - we don't know which value in $J$ is the true response (!!!)
            - hence, 
                - *equ 2.1* must be reformulated as *equ 2.6*: \Pr(L^{(i)}|p_j,\boldsymbol{\pi},\mathbf{n}_{i,l}) \propto \sum_{j=1}^J\left\{p_j\prod_{k=1}^K\prod_{l=1}^J\left(\pi_{q,l}^{(k)}\right)^{n_{i,l}^{(k)}} \right\}
                    - *equ 2.6* is a mixture of multinomial distributions with the marginal probabilities $p_j$ acting as weights
                - and the full likelihood of the data is (*equ 2.7*) \prod_{i=1}^I\left(\sum_{j=1}^J\left\{p_j\prod_{k=1}^K\prod_{l=1}^J\left(\pi_{q,l}^{(k)} \right)^{n_{i,l}^{(k)}} \right\}\right)$
            - treated the indicator variables ${T_{ij}: i \in I, j \in J}$ as missing data, and apply EM algorithm:
                1. take initial estimates/set initial values of the $T$s
                2. use equations 2.3 and 2.4 to obtain estimates of the $p$s and $\pi$s. 
                3. use equation 2.5 and the estimates of the $p$s and $\pi$s to calculate new estimates of the $T$s as $p(T_{ij} | L^{(i)}, p_j, \boldsymbol{\pi}, \mathbf{n}_ {i,l})$
                4. Repeat steps two and three until the results converge
            - final estimates obtained by the EM algorithm are those that
                - maximize *equ 2.7*,
                - and provided the ground for assigning patients to response classes (\hat{T} are "consensus probabilities", p. 24)  
    - *discussion* (pp. 24f.)
        - comments:
            1. in any application context at hand, "the concept of a true response must be meaningful" (e.g., as a hypothetically observable quantity, i.e., a latent trait)
            2. as a latent class model, the annotation model *lacks identifiability* 
                - given data, the likelihood defined in *equ 2.7* remains unchanged for any relabeling of indexes $j \in J$
                - there correspond $J!$ sets of estimates to the global maximum
                - potential remedies:
                    - verify that \pi_{jj}^{(k)} > \pi_{jl}^{(k)} j \neq l$ for most $k$ and $j$ for different runs of the algorithm
                    - choose a better *initial/starting values* for $T_{ij}$ than $1/J$ (for this is the 'center of symmetry' and will let the algorithm start at a saddle point of the likelihood surface, p. 24)
                        - e.g., the average value or the majority-winner assignment at the patient/item level
                - rerun the algorithm with different initial/starting values to verify the stability of estimates (especially when few data is used to use many parameters)
            3. the model does not provide for estimates of uncertainty in the parameters
                - use leave-observer(s)-out or leave-patient(s)-out validation to obtain a bootstrap estimate of error margins
            4. "unrealistically accurate estimates of the $T$'s are obtained" (p. 25)
                - again, use leave-one-out validation/bootstrapping to examine the level of certainty/stability in class assignment (i.e., "to establish which patients are well classified and for which patients the consensus is more uncertain", p. 25)
                - uncertain cases (i.e., shaky class assignments) will be sensitive to small changes in error-rate estimates $\pi_{jl}^{(k)}$), hence the rationale of leave-patient-out validation (cf. p. 27)



