# Johnson–Lindenstrauss theorem(约翰逊-林登斯特劳斯定理)以及随机投影对其的利用

# 一 **Johnson–Lindenstrauss theorem(约翰逊-林登斯特劳斯定理)**

![Untitled](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-07-24-Johnson-Lindenstrauss_theorem/picture1.png)

![Untitled](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-07-24-Johnson-Lindenstrauss_theorem/picture2.png)

![Untitled](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-07-24-Johnson-Lindenstrauss_theorem/picture3.png)

[让人惊叹的Johnson-Lindenstrauss引理：理论篇 - 知乎 (zhihu.com)](https://zhuanlan.zhihu.com/p/413581747)  //作者证明地很好

- 需要的引理：
    - 马尔可夫不等式（切比雪夫不等式）
        
        ![Untitled](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-07-24-Johnson-Lindenstrauss_theorem/picture4.png)
        
        其中切比雪夫不等式不等号右边的分子就是方差
        
    - Cramér-Chernoff方法
        
        ![Untitled](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-07-24-Johnson-Lindenstrauss_theorem/picture5.png)
        
    - 单位模引理
        
        ![Untitled](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-07-24-Johnson-Lindenstrauss_theorem/picture6.png)
        
        注意上面推导中两个画圈的地方：
        
        - 第一个画圈的地方，作者大概是书写错误了，但是计算结果是正确的，简单拿正太分布的概率密度函数计算并带入$\int_{-\infin}^{+\infin} e^{-x^2}dx = \sqrt{\pi}$已知积分计算结论计算即可
        - 第二个画圈的地方，需要理解以下式子（容斥原理）：
            
            $$
            P(| \space ||u_i||^2 -1 \space | \ge \varepsilon) = P(||u_i||^2 -1 \ge \varepsilon 或 1-||u_i||^2 \ge \varepsilon)=P(||u_i||^2 -1 \ge \varepsilon)+P(1-||u_i||^2 \ge \varepsilon) - P(||u_i||^2 = \varepsilon + 1 ) \le P(||u_i||^2 -1 \ge \varepsilon)+P(1-||u_i||^2 \ge \varepsilon) \le 2 e^{-\frac{n\varepsilon^2}{8}}
            $$
            
    - 正交性引理（作者自己叫的）
        
        ![Untitled](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-07-24-Johnson-Lindenstrauss_theorem/picture7.png)
        
        注意那个两式相加，推导如下（补充高中知识）
        
        使用余弦定理：
        
        ![屏幕截图 2024-01-10 213343.png](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-07-24-Johnson-Lindenstrauss_theorem/picture8.png)
        
- 加下来证明JL定理本身：
    
    ![Untitled](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-07-24-Johnson-Lindenstrauss_theorem/picture9.png)
    
    ![Untitled](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-07-24-Johnson-Lindenstrauss_theorem/picture10.png)
    
    结论：
    
    - **如果有$N$个向量，那么只需要存储$logN$级别的维度就可以很好地表示其中任意两个向量之间的欧氏距离（夹角）。**
    - 如果是按照欧式距离来定义的，那么降维后两个向量相减的距离（欧氏距离）会在原欧氏距离的$(1-\varepsilon,1+\varepsilon)$倍之间
    - 如果是按照夹角来定义的，那么降维后两个向量的夹角大小的范围是$(-\varepsilon+\theta,\varepsilon+\theta)$，其中$\theta$是降维前原来两个向量的夹角大小。

Johnson, William B; Lindenstrauss, Joram. Extensions of Lipschitz mappings into a Hilbert space. Contemporary mathematics. 1984, **26** (1): 189-206.

# 二、随机投影

随机投影是一种基于Johnson-Lindenstrauss引理的降维技术。具体过程如下：

1. 利用上面的JL定理，假设向量数据库$D$中有$N$个数据点，并且维度大小为$m$，那么生成一个投影矩阵$A_{nm}$，它的每一个元素$a_{ij}$独立重复采样自高斯分布($a_{ij} \sim N(0,\frac{1}{n})$)，设置要降到的维度为$n$，并且$n$满足：$n\ge \frac{24logN}{\varepsilon^2}$，其中$\varepsilon$是可以调节的参数。
2. 利用上述已经设置好的变量，对于$\forall (v_i)_{1m},(v_j)_{1m} \in D$，将其通过投影矩阵降维：$(Av_i)_{1n},(Av_j)_{1n}$后，**则至少有$\frac{N-1}{N}$（这个可以根据$n$和$\varepsilon$的大小来调节，是一个放缩后的最小值）的概率，使得降维后的向量之间差值的模（$||Av_i - Av_j ||$）在原来向量差值模（$||v_i-v_j||$）的$(1-\varepsilon,1+\varepsilon)$倍之间。**
3. **通过调节$n$和$\varepsilon$的大小来得到我们想要的发生概率大小结果。**
4. 降维后，原先的维度灾难能得到一定程度上的缓解。

利用随即投影预处理数据的时间复杂度：

预设一个投影矩阵的时间复杂度大概为$O(nm)$，做投影变换的时间复杂度也是$O(nm)$，是一个与维度相关的平方复杂度。