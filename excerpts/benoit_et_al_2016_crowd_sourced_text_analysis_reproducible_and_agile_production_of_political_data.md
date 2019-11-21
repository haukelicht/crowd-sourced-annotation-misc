<!-- benoit_et_al_2016_crowd_sourced_text_analysis_reproducible_and_agile_production_of_political_data.md -->
- Benoit et al. (2016) "Crowd-sourced Text Analysis: Reproducible and Agile Production of Political Data"
    - introduce a methodological framework for the crowd-sourced annotation of text corpora and other labeling exercises in large amounts of  unstructured (e.g., textual) data in political science
    - they argue that their framework 
        - is both *reproducible* and *agile* (i.e., scalable and adaptable to different specific research interests)
        - facilitates the *reproducibility* of data (collection) rather than merely the replicability of data analyses (p. 278)
        - allows for imposing some structure in the iterative development of coding schemes in data annotation tasks (in a deploy-test/verify-(re)design cycle, p. 279)
    - provide a proof of concept: devise an estimator of the policy domain of sentences from political texts (economic: left vs. right, vs. social: liberal vs. conservative) that is based on crowd-sourced classification of this data
        - examine the validity of their approach externally (), as well as internally ()

    - *introduction* (pp. 278f.)
    - *harwesting the wisdom of crowds* (pp. 279f.)
        - crowd-sourcing:= using the Internet to distribute a large package of small tasks to a large number of anonymous workers
        - argue that, "[Provided] crowd workers are not systematically biased in relation to the 'true' value of the latent quantity of interest ..., the central tendency of even erratic workers will converge on this true value as the number of workers increases" (p. 279)
        - stress that "it is important to check for such bias" (p. 279)
        - stress the importance of having a 'gold standard' (an 'objective truth') to evaluate how crowd-workers perform
            - point to literature on Bayesean scaling/item-response models of annotation that allow determine such a gold standard and evaluate annotators performance characteristics simultaneously (p. 280)
    - *a method for replicable coding of political text* (pp. 280-2)  
        - portray the systematic extraction of meaning from a text corpus by human analysts as a canonical problem in many large-scale political science research projects (e.g., the CMP)
            - such methods depend on the reading and semantic interpretations of humans (i.e., "human analysts are employed to engage in *natural language processing* (NLP)", p. 280)
            - stress the pervasive shortage of human expert coders as "highly sophisticated natural language processors" (p. 280)
        - stress that many supervised learning methods depend on the availability of labeled training datasets
        - use the task of classifying sentences from political manifestos into classes covering economic or social policy and their positions on those dimensions as a *running example* and proof of concept (coding scheme and text corpus are described on pp. 281f.)
    - *scaling document policy positions from coded sentences* (pp. 282f.)
        - the goal of their running example "is to estimate the policy positions of entire documents" (i.e., manifestos)
            - when striving for this goal, they emphasize that aggregating the values (classes and ratings) of sentences to estimate manifestos' class membership and positions warrants to take into account "reader, sentence, and domain effects" (i.e., differences in annotators qualities, intrinsic sentence difficulty/ambiguity, and uncertainty in point estimates)
        - portray simple averaging (SA) as "the benchmark when there is no additional information on the quality of individual coders" (p. 282)
        - stress that using scaling methods based item response theory provide a feasible alternative to SA
        - model 
            - each sentence $j$ as a vector of four sentence attributes: *domain propensity* to be labeled (i) economic or (ii) social, and *latent position* of the sentence on (iii) the economic or (iv) the social dimensions
            - individual readers (annotators) $i$ have generic *biases* in each of these four dimensions (assumed to be independent from sentence), and four *sensitivities* (i.e., "their relative responsiveness to changes in the latent sentence attributes in each dimension", p. 282) 
            - the *latent coding* of a sentence $j$ by reader $i$ on dimension $d$ is then the sum of the latent sentence attribute (domain propensity or position) plus the reader's individual bias on this dimension, times the reader's sensitivity on this dimension
            - responses to domains are modeled as distributed multinominally on {"economic", "social", "neither"}, and choices of scale position as ordinal logit
        - find that the posterior means of the document level attributes strongly correlate with SA of readers' labels
    - *benchmarking a crowd of experts* (pp. 283-5)
        - external validity: examine whether it makes a difference if expert coders process documents sentence-by-sentence in sequential order versus random sentence order, and find that random-order values/ratings correlate strongly with sequential-order values/ratings
            - stress that serving sentences at random has the advantage of scalability
        - internal validity:
            - examine inter-expert agreement
                - find perfect inter-expert agreement on sentence domains and domain positions in 35% of all sentences
                - there even occurred some agreement on policy domain, but agreement on sentences' policy domain is high (Fleiss' kappa is .93)
                - stress, however, that aggregating rating at the document level yields externally valid results
            - scale reliability: Cronbach's alpha (measures of scale reliability across readers) is excellent (.95)
    - *deploying crowd sourced text coding* (pp. 286f.)
        - *HIT*:= human intelligence task
        - *CrowdFlower* (CF): a platform that consolidates access to dozens of channels such as Amazon's Mechanical Turk
            - regularly implements *quality control*:
                - filters out 'spammers' by  using gold-standard based or other *screening test* prior to assigning tasks and intermittently during task completion
                - "gold HITs" can be seeded into the data to continuously verify and measure workers' performances on tasks
                - "screeners" (sentences containing literal instructions, e.g. "select category X and click continue") can be used to check whether workers pay attention
            - the authors use a two-stage process of quality control in their running example (p. 286)
        - deployment:
            (i) oversampled sentences from specific documents (manifestos from 1987 to '97) to be later able to examine the effect of repeated labeling on estimates as a function of the number of annotators of sentences
            (ii) sampled remaining samples until had received five labels (also avoids tie-breaking in classification problems)
    - *crowd-sourced estimates of party policy positions* (pp. 287)   
        - high correlations between CMP-expert and crowd-coding estimates of (average) document positions on both policy domains 
        - stress that point estimates of scaling models and SA values are strongly correlated and go for SA because scale values are contingent on the datasets included (p. 287)
        - show that "uncertainty over the crowd-based estimates collapses as we increase the number of workers per sentence"
            - argue that this result hinges on the fact that estimates are derived by averaging values over sentences at the document level ("With five scores per sentence and about 1000 sentences per manifesto, we have about 5000 'little' estimates of the manifesto position ... This sample is big enough to achieve a reasonable level of precision, given the large number of sentences per manifesto", p. 289)
            - stress that for shorter documents (i.e., documents with fewer sentences), or analyses that do not aggregate at the document level, the number of workers pro sentence need to obtain solid estimates may be higher (p. 289) 
    - *crowd-sourcing data for specific projects: immigration policy* (pp. 289-)
        - use the same method to derive estimates of documents' positions in the domain of immigration policy
        - use an "adaptive sentence sampling strategy" that samples a minimum of three sentences and an additional two "unless the first three [annotations] were unanimous in judging a sentence *not* to concern immigration policy" (p. 290)  
            - argue that this approach focuses resources on the few sentences that are on immigration policy
                - [hli]: note however, that the error incurred by this approach increases as readers annotation quality/ability decreases ("3 choose 0" is the higher if p, readers accuracy, is lower, i.e., three coders may agree but they still may all be wrong with some none-zero probability)
    - *crowd sourced text analysis in other contexts and languages* (p. 290f.)
        - report promising results from examining inter-language correlations
            - use EU debate speeches, their corresponding official translations into different languages and voting outcomes to examine whether scores derived from reading in different languages correlate, and how well they predict voting outcomes
        





