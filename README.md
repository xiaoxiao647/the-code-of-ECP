# Elastic commensurate prior（ECP）：模拟实验代码文件说明
此文件夹包含用于对弹性相称先验（ECP）进行模拟研究的R代码，包括两个子文件夹，分别是正态终点（normal）和二分类终点（binary）的代码。每个子文件夹包含4份代码，分别用于使用ECP从单个外部数据借用信息，从多个外部数据借用信息，在纳入协变量信息共同度量数据异质性的情形下从外部数据借用信息的模拟研究，以及优化弹性函数g(T)=exp⁡(a+b∙log⁡(T))参数a和b的网格搜索。

*ECP_normal_model.R
该代码用于在正态终点下，使用ECP从单个外部数据借用信息的模拟研究，生成手稿中表4的模拟结果。它包含四个函数：
1）decide_para()用于确定弹性函数g(T)=exp⁡(a+b∙log⁡(T))中的参数a和b；
2）sample_poster()用于对对照组均值的后验进行抽样；
3）ECP_normal_model()用于获得每次试验里试验组疗效优于对照组疗效的概率；
4）rej_null_p()用于计算拒绝零假设，声称药物优效的概率。

* ECP_multiple_normal_model.R
该代码用于在正态终点下，使用ECP从多个外部数据借用信息的模拟研究，生成手稿中表6的模拟结果。它包含三个函数：

1）decide_para()用于确定弹性函数g(T)=exp⁡(a+b∙log⁡(T))中的参数a和b；

2）ECP_multiple_normal_model()用于获得每次试验里试验组疗效优于对照组疗效的概率；

3）rej_null_p()用于计算拒绝零假设，声称药物优效的概率。

* ECP_normal_covariates.R
该代码用于在正态终点下，使用ECP从单个外部数据借用信息的模拟研究，生成手稿中表8的模拟结果。它包含四个函数：

1）decide_para()用于确定弹性函数g(T)=exp⁡(a+b∙log⁡(T))中的参数a和b；

2）sample_poster()用于对对照组均值的后验进行抽样；

3）ECP_cov_normal()用于获得每次试验里试验组疗效优于对照组疗效的概率；

4）rej_null_p()用于计算拒绝零假设，声称药物优效的概率。

* grid_search_normal.R
该代码用于在正态终点下，使用网格搜索确定q1和q2的取值，从而优化弹性函数g(T)=exp⁡(a+b∙log⁡(T))参数a和b的取值。它包含四个函数：

1）decide_para()用于确定弹性函数g(T)=exp⁡(a+b∙log⁡(T))中的参数a和b；

2）sample_poster()用于对对照组均值的后验进行抽样；

3）ECP_normal_model()用于获得每次试验里试验组疗效优于对照组疗效的概率；

4）rej_null_p()用于计算拒绝零假设，声称药物优效的概率。

* ECP_binary_model.R
该代码用于在二分类终点下，使用ECP从单个外部数据借用信息的模拟研究，生成手稿中表5的模拟结果。它包含三个函数：

1）decide_para()用于确定弹性函数g(T)=exp⁡(a+b∙log⁡(T))中的参数a和b；

2）ECP_binary_model()用于获得每次试验里试验组疗效优于对照组疗效的概率；

3）rej_null_p()用于计算拒绝零假设，声称药物优效的概率。

* ECP_multiple_binary_model.R
该代码用于在二分类终点下，使用ECP从多个外部数据借用信息的模拟研究，生成手稿中表7的模拟结果。它包含三个函数：

1）decide_para()用于确定弹性函数g(T)=exp⁡(a+b∙log⁡(T))中的参数a和b；

2）ECP_multiple_binary_model()用于获得每次试验里试验组疗效优于对照组疗效的概率；

3）rej_null_p()用于计算拒绝零假设，声称药物优效的概率。

* ECP_binary_covariates.R
该代码用于在二分类终点下，使用ECP从单个外部数据借用信息的模拟研究，生成手稿中表9的模拟结果。它包含三个函数：

1）decide_para()用于确定弹性函数g(T)=exp⁡(a+b∙log⁡(T))中的参数a和b；

2）ECP_cov_binary()用于获得每次试验里试验组疗效优于对照组疗效的概率；

3）rej_null_p()用于计算拒绝零假设，声称药物优效的概率。

* grid_search_binary.R
该代码用于在二分类终点下，使用网格搜索确定q1和q2的取值，从而优化弹性函数g(T)=exp⁡(a+b∙log⁡(T))参数a和b的取值。它包含三个函数：

1）decide_para()用于确定弹性函数g(T)=exp⁡(a+b∙log⁡(T))中的参数a和b；

2）ECP_binary_model()用于获得每次试验里试验组疗效优于对照组疗效的概率；

3）rej_null_p()用于计算拒绝零假设，声称药物优效的概率。
