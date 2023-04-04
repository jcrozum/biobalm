# Meeting on 04/04/2023

## Issue 44: Single-node driver and terminal restriction space optimization

<https://github.com/jcrozum/nfvs-motifs/issues/44>

I also think that this improvement can dramatically reduce the size of the candidate set in some cases.

## Issue 46: Simplified succession diagram

<https://github.com/jcrozum/nfvs-motifs/pull/46>

### Question 1

I found that the full succession diagram construction (in the benchmarks) does not consider the motif-avoidant attractor checking. I just want to confirm with you.

### Question 2

Could you please explain briefly the idea of the symbolic pruning algorithm?

### Question 3

I understand that the symbolic pruning algorithm can lead to two cases. Case 1: the "proof" is succeeded, i.e., there is no motif-avoidant attractor. Case 2: the "proof" is failed.

Case 1: continue

Case 2: In this case, the succession diagram is conservatively expanded to cover all options. So, how the symbolic pruning algorithm can detect motif-avoidant attractors?

### Question 4

I think with the simplified succession diagram, our control algorithms can still work as usual. But, I am not sure.

### Question 5

I agree that this algorithm can be a baseline comparison for future faster methods, although I am afraid that it will be not efficient for complex models, for example, N-K models because it uses BDDs, whereas BDD-based approaches seems to be not efficient for N-K models (from my previous papers). We should test more and see the results.

### Question 6

I am not clear about your suggestion

> We can replace the symbolic pruning with a trappist check of minimal trap spaces. This should ensure we don't spend too much time in the symbolic pruning step, even with the complexity bound, but it means we no longer prove the absence of motif-avoidant attractors. So it may be worth it if you are only interested in the succession diagram itself, but not if you then also want to find all attractors.

Could you please explain more?

### Question 7

I think we have all necessary functions implemented. Now, we first need to specify what are the target features of our tool. For example,

+ Build the full succession diagram without attractor checking
+ Build the full succession diagram with attractor checking
+ Build the partial succession diagram with respect to the target minimal trap space (for control)
+ ...

We need to prepare a specification document for this.

### Question 8

For the 211 model, `mtsNFVS` can handle it within 324.35s. Hence, it is possible that `nfvs-motifs` can handle this model with much smaller time. We need to optimize `nfvs-motifs` and see the results.

Anyway, the symbolic pruning algorithm already appears very great. Now, I see that it is much better than both `AEON` and `mtsNFVS` in most cases.

## Control problem

I have some concerns on the target control problem.

First, in some papers (e.g., Cifuentes Fontanals, L., Tonello, E., & Siebert, H. (2020). Control strategy identification via trap spaces in Boolean networks. In _Computational Methods in Systems Biology: 18th International Conference, CMSB 2020, Konstanz, Germany, September 23–25, 2020, Proceedings 18_ (pp. 159-175). Springer International Publishing.), the authors claim that `pystablemotifs` is only applicable for the fully asynchronous update scheme. Let see an example in Figure 2 of this paper. The succession diagram of this Boolean networks is as follows. 

Initially, I thought the result for the target control should be `{x3 = 0, x1 = 1}`.  However, `pystablemotifs` returns one control strategy `{x3 = 0}`. I found that `pystablemotifs` detects there is only one path from node `{x3 = 0}` to the minimal trap space `110`, so it excludes node `{x1 = 1}` from the consideration. The problem is that this reduction depends on the update scheme. There may be a motif-avoidant attractor inside the subspace `{x3 = 0}`.

Hence, I think for the control of any update scheme, the result should be `{x3 = 0, x1 = 1}`. If we consider the fully asynchronous update scheme, the result should be `{x3 = 0}`.

From this example, it raises some issues that need to be thoroughly considered.

+ If we consider the fully asynchronous update scheme and there are some motif-avoidant attractors, how our control algorithms work?
+ `pystablemotifs` and the trap space-based method (in the case of fully asynchronous) may miss some control strategies. With the above example Boolean network, they miss `{x2 = 0}` that can be detected by `CAEABN`. My question is how we can use the information about the existence of motif-avoidant attractors to find more possible control strategies in our method (i.e., `nfvs-motifs`). Of course, we need to carefully consider two cases: no motif-avoidant attractors and some motif-avoidant attractors. Note that the position of a motif-avoidant attractor in the succession diagram may be important.

Now, I suggest that we should migrate the control features of `pystablemotifs` to `nfvs-motifs` first. Then, we will discuss in detail how to improve these features in `nfvs-motifs`. I have some ideas for this. I will show you some progress soon.

## Experiments

### Compared tools

#### On attractor detection

+ `CABEAN`

> Mizera, A., Pang, J., Qu, H., & Yuan, Q. (2018). Taming asynchrony for attractor detection in large Boolean networks. _IEEE/ACM transactions on computational biology and bioinformatics_, _16_(1), 31-42.

+ `AEON`

> Beneš, N., Brim, L., Pastva, S., & Šafránek, D. (2021). Computing bottom SCCs symbolically using transition guided reduction. In _Computer Aided Verification: 33rd International Conference, CAV 2021, Virtual Event, July 20–23, 2021, Proceedings, Part I 33_ (pp. 505-528). Springer International Publishing.

+ `PyBoolNet`

> Klarner, H., & Siebert, H. (2015). Approximating attractors of Boolean networks by iterative CTL model checking. _Frontiers in bioengineering and biotechnology_, _3_, 130.

+ `pystablemotifs`

> Rozum, J. C., Gómez Tejeda Zañudo, J., Gan, X., Deritei, D., & Albert, R. (2021). Parity and time reversal elucidate both decision-making in empirical models and attractor scaling in critical Boolean networks. _Science Advances_, _7_(29), eabf8124.

+ `iFVS-ABN`

> Van Giang, T., & Hiraishi, K. (2021, October). An Improved Method for Finding Attractors of Large-Scale Asynchronous Boolean Networks. In _2021 IEEE Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB)_ (pp. 1-9). IEEE.

+ `mtsNFVS`

> Trinh, V. G., Hiraishi, K., & Benhamou, B. (2022, August). Computing attractors of large-scale asynchronous boolean networks using minimal trap spaces. In _Proceedings of the 13th ACM International Conference on Bioinformatics, Computational Biology and Health Informatics_ (pp. 1-10).

+ `nfvs-motifs` (may have multiple variants)

#### On target control

It should be divided into two categories: any update scheme and fully asynchronous.

Any update scheme:

+ `PyBoolNet` (trap space-based method without considering attractors)

> Cifuentes Fontanals, L., Tonello, E., & Siebert, H. (2020). Control strategy identification via trap spaces in Boolean networks. In _Computational Methods in Systems Biology: 18th International Conference, CMSB 2020, Konstanz, Germany, September 23–25, 2020, Proceedings 18_ (pp. 159-175). Springer International Publishing.

> Fontanals Laura, C., Tonello, E., & Siebert, H. (2022, November). Computing trap space-based control strategies for Boolean networks using answer set programming. In _AIP Conference Proceedings_ (Vol. 2611, No. 1, p. 110002). AIP Publishing LLC.

+ `pystablemotifs`

> Rozum, J. C., Deritei, D., Park, K. H., Gómez Tejeda Zañudo, J., & Albert, R. (2022). pystablemotifs: Python library for attractor identification and control in Boolean networks. _Bioinformatics_, _38_(5), 1465-1466.

+ `DCGS`

> An, S., Jang, S. Y., Park, S. M., Lee, C. K., Kim, H. M., & Cho, K. H. (2023). Global stabilizing control of large-scale biomolecular regulatory networks. _Bioinformatics_, _39_(1), btad045.

+ `nfvs-motifs` (may have multiple variants)

Fully asynchronous:

+ `PyBoolNet` (trap space-based method with considering attractors)

> Cifuentes Fontanals, L., Tonello, E., & Siebert, H. (2020). Control strategy identification via trap spaces in Boolean networks. In _Computational Methods in Systems Biology: 18th International Conference, CMSB 2020, Konstanz, Germany, September 23–25, 2020, Proceedings 18_ (pp. 159-175). Springer International Publishing.

> Fontanals Laura, C., Tonello, E., & Siebert, H. (2022, November). Computing trap space-based control strategies for Boolean networks using answer set programming. In _AIP Conference Proceedings_ (Vol. 2611, No. 1, p. 110002). AIP Publishing LLC.

+ `PyBoolNet` (model checking-based method)

> Cifuentes-Fontanals, L., Tonello, E., & Siebert, H. (2022). Control in Boolean networks with model checking. _Frontiers in Applied Mathematics and Statistics_, _8_, 23.

+ `CABEAN`

> Su, C., & Pang, J. (2021). Target control of asynchronous Boolean networks. _IEEE/ACM Transactions on Computational Biology and Bioinformatics_.

+ `pystablemotifs`

> Rozum, J. C., Deritei, D., Park, K. H., Gómez Tejeda Zañudo, J., & Albert, R. (2022). pystablemotifs: Python library for attractor identification and control in Boolean networks. _Bioinformatics_, _38_(5), 1465-1466.

+ `nfvs-motifs` (may have multiple variants)

### Compared metrics

On attractor detection

+ ~~Accuracy~~
+ Running time

On target control

+ Minimum size of control strategies
+ Number of returned control strategies
+ Running time

The case of target control seems to be quite complicated. Needs to be considered thoroughly.

### Fairness

Regarding the comparison, the performance of `AEON`, `mtsNFVS`, and `nfvs-motifs` largely depends on the network structure. Hence, it should be fair if we perform network reductions on the original models before running the actual methods. For example,

+ leaf-node reduction
+ constant-node propagation

For attractor detection, there is no issue. But we need to be careful with the case of target control.

## GitHub issue

> I also wanted to give everyone a heads up that I was just notified that we are at 75% of the monthly GitHub actions usage limit for private repos (it's a cpu-time allotment as far as I can tell). It resets on April 13th. Otherwise we will have to pay for more time or take the repo public to continue using our CI. I think we should consider reducing the computational demands of our automated tests, at least on the server side, to avoid this limit.

## Next meeting

I will have to attend my lab's seminar next Tuesday. Hence, can we meet next next Tuesday?