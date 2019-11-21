<!-- snow_et_al_2008_cheap_and_fast_but_is_it_good_evaluating_non_expert_annotations_for_natural_language_tasks.md -->
- Snow et al. (2008) "Cheap and Fast—But is it Good? Evaluating Non-Expert Annotations for Natural Language Tasks"
    - consider crowd-sourced annotation a "promising alternative [to] ... the collection of non-expert annotations" for tasks involving human natural language processing (NLP)
    - examine whether crowd-sourced annotation of conventional NLP datasets in linguistics can produce data of quality that matches that of trained/expert coders
    - *related work* (p. 255)
    - *task design* (pp. 256)
        - describe the Amazon Mechanical Turk platform for distributing crowd-sourced human intelligence tasks (HITs)
    - *annotation tasks* (pp. 256-60)
        - compare crowd-worker to expert performances on five datasets where annotation task involve only multiple-choice responses or numeric input (p. 256):
            1. Affective Text Analysis (ATA) due to Strapparava and Mihalcea (2007)
                - annotators are asked to rank short headlines on an integer scale ranging from [-100,100] on six emotions: anger, disgust, fear, joy, sadness, and surprise 
                - examine interannotator agreement (ITA) among and between experts and non-experts (Table 1, p. 256)
                - demonstrate that by and large nine non-expert coders suffice to achieve the level of annotation quality produced by expert annotation (Table 2 and Figure 1, p. 257)
            2. Word Similarity (WS) due to Miller and Charles (1991)
                - annotators are asked for numeric judgments of word similarity for 30 word pairs on a continuous scale of [0,10]
                - demonstrate that ten annotators suffice to obtain an average per-item agreement of .93 (ITA in a study with 51/10 subjects produced an ITA of .97/.958, p. 257)
            3. Recognizing Textual Entailment (RTE) due to Dagan et al. (2006)
                - annotators are presented with a pair of sentences, and asked to assess whether the second sentence can be inferred from the information provided in the first sentence (i.e., by ignoring 'world knowledge')
                - collect 10 annotations 100 RTE tasks each and achieved ITA of .897 (expert annotation has achieved 91% and 96% ITA over various subsections of the corpus, p. 258)
            4. Event Annotation (EA) due to Pustejovsky et al. (2003)
                - experts are asked to annotate event pairs by temporal relation between, choosing labels from a set of fourteen relations (before, after, during, includes, etc.)
                - the tasks was simplified for crowd-workers: could select only the relations 'strictly before' or 'strictly after'
                - achieve an ITA of .94 by aggregating 10 non-expert annotations (compared to an ITA of .77 for the more complex task faced by experts)
            5. Word Sense Disambiguation (WSD) due to Pradhan et al. (2007)
                - annotators are presented with a sentence that contains a word with multiple synonyms, and annotators are asked to identify the 'correct' meaning of the word in this context
                - crowd-source 10 annotations for each of 177 examples of the noun "president" for three senses (executive officer of an organization, head of state, or head of the U.S.)
                - compare crowd-sourced ITA with the best automatic system performance for the disambiguation task (accuracy of .98)
                    - simple majority voting achieves an accuracy of .994 already with only 3 annotators
    - *bias correction for non-expert annotators* (pp. 260f.)
        - argue that, due to the variability in annotator quality and the prospect of identifying high-volume but low-quality workers, "controlling for individual worker quality could yield higher quality overall judgments" (p. 260)
        - propose to use a small amount of expert-labeled training data as a gold standard in order to correct for the individual biases of different non-expert annotators
            - crowd-worker response likelihoods $\Pr(y_w|x = Y )$ and $\Pr(y_w|x = N)$ are estimated from frequencies of worker performance on gold standard examples
            - the posterior probability of the true label for an item is $\sum_w \log\frac{\Pr(y_{i,w}|x_i = Y)}{\Pr(y_{i,w}|x_i = N)} + \log\frac{\Pr(x_i = Y)}{\Pr(x_i = N)}$
            - this amounts to a *weighted voting rule*: "each worker's vote is weighted by their log likelihood ratio for their given response" (p. 261)
        - demonstrate that bias correction when determining labels improves accuracy compared to 'naive' majority voting of the in the RTE task
    - *conclusion* (p. 262)
        - demonstrate that for many NLP tasks only a small number of non-expert annotations per item suffice necessary to match the performance of single expert annotators
        - argue that significant improvements in crowd-sourced annotation can be achieve by controlling for annotator biases when aggregating annotations into labels


