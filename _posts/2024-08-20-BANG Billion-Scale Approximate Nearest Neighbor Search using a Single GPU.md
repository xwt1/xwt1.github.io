# BANG: Billion-Scale Approximate Nearest Neighbor Search using a Single GPU

# 1. Introduction

文章的主要贡献：

1. 作者的贡献是提出了BANG，这种架构能够使用各种之前已有的技术，如PQ等压缩数据。
2. 同时发明了多种GPU优化技术如指令预取，GPU与CPU之间的PCIE争抢，使得能够在不牺牲recall的前提下进行检索。

# 2. Background

## 2.1 GPU Architecture and Programming Model

简要的背景知识：

- GPU与CPU相互通信是依靠PCI-E总线进行通信
- **GPU架构**：
    - GPU是为高并行性和内存带宽而设计的专用加速器。
    - 它们由多个CUDA核心组成，这些核心被分组到流式多处理器（SM）中。
    - 每个SM共享片上共享内存，作为L1缓存使用。
- **编程模型**：
    - 使用的编程模型是CUDA，它是NVIDIA创建的并行计算平台和应用程序编程接口模型。
    - GPU的主内存称为全局内存。
    - GPU任务（称为内核）在SM上并行执行，使用线程网格。
    - 线程被组织成块，每个块分配给一个SM。
    - 在线程块内，线程可以使用共享内存进行通信和同步。
    - 线程块包含多个warps（每个warp包含32个线程），以SIMD（单指令多数据）方式执行。
- **任务并行性**：
    - 任务队列（称为流）用于内核调用和内存传输，从而实现任务并行性

## 2.2 A summary of DiskANN and Vamana Grap

老朋友不用介绍，直接看以往的记录

## 2.3 Vector Compression

老朋友不用介绍，直接看以往的记录

# 3. ANN SEARCH ON GPU

## 3.1 Challenges in Handling Billion-scale Data

文章中提到了在GPU上做检索最大的两个限制是：

- **有限的GPU存储限制：**在使用图索引时，会遭遇GPU内存有限的问题，图索引过大，无法将全部数据放入GPU显存中的问题
    - 解决方法：在CPU上利用图索引做搜索，在GPU利用乘积量化做距离计算。
    - 好处
        - **压缩数据传输**：将数据压缩后再传输到GPU，减少了数据传输量，从而降低了带宽需求。
        - **按需获取邻居信息**：在需要时才从CPU传输节点的邻居信息到GPU，而不是一次性传输整个图的数据。
- **优化的硬件使用：**使用多线程时现行办法无法充分利用系统资源，diskANN的方法无法直接被应用到CPU-GPU异构系统中，原因有两个：**第一个，基于GPU的检索，需要平衡CPU与GPU之间的工作流（这在CPU与DISK之间不会出现）。第二个：GPU的ANN检索有着更明显的内存占用量，因为它能够在GPU上并行解决大量的查询（GPU比起CPU更优良的并行计算能力），所以需要大量与这些查询有关的信息（如索引中的内容），而传输这些信息需要占用大量PCIe接口的带宽导致成为搜索速度的瓶颈**。
    
    更多的，假设CPU保存图索引本身，而GPU用于压缩向量之间的距离计算，在这些前提下，新系统需要确保CPU和GPU都能连续工作，避免出现空闲时间。
    
    并且，对于一系列需要执行并行查询的查询集合$Q$，与其有关的图上的结点集合$V$，在GPU上进行并行查询时，如何管理用于跟踪访问节点的数据结构的问题。具体来说，**作者提到需要决定这些数据结构应该放置在CPU还是GPU上**。为了设计一个高效的工作流，必须考虑几个关键因素，**包括与数据结构的内存大小相关的计算复杂度、在访问节点时的内存访问模式**（CPU和GPU的内存带宽不同导致访问时间不同），**以及在CPU和GPU之间通过PCIe互连传输数据所需的时间**。作者强调，设计工作流时要特别注意这些因素，以便在带宽有限的情况下实现最佳性能。
    
    为了解决上面的一系列挑战，作者提出了BANG的解决方案：
    
    - 解决方法：
        - **最大化并行处理**：BANG通过高效利用硬件资源（CPU、GPU和PCIe总线）来实现最大化的并行处理。
        - **负载均衡**：BANG在CPU和GPU之间有效地分配ANN搜索工作，以实现负载均衡。
        - **优化技术**：采用了预取和流水线技术，管理CPU和GPU的空闲时间。
        - **布隆过滤器**：集成了布隆过滤器，用于处理大图上多个查询的访问节点集合。
        - **GPU优化内核**：对GPU上的关键操作（如距离计算、排序和合并）进行了高度优化。

## 3.2  BANG Overview

![image.png](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-08-20-BANG Billion-Scale%20Approximate%20Nearest%20Neighbor%20Search%20using%20a%20Single%20GPU/image.png)

上图是BANG的整体架构图。

BANG使用了Vanama作为基础索引，由于其体积很大，将其存储在CPU内存中。距离比较计算在GPU中，并且保证使用压缩向量进行距离比较计算。

约定：设原始数据有$d$维度，乘积量化后，一共有$m$个子空间，每个子空间中有$c$个聚簇（向量），每个子空间中的code的类型是$uint8$，也就是一个字节。

BANG整体包含三个阶段，如架构图中黑色标号所标注的那样：

1. **Distance Table Construction**
    1. 在开始搜索之前，首先需要提前计算每一个数据点，对应每一个子空间中的对应聚簇与query对应子空间向量之间的距离，并且将这些距离制作成一张PQ Distance table。
        
        ![image.png](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-08-20-BANG Billion-Scale%20Approximate%20Nearest%20Neighbor%20Search%20using%20a%20Single%20GPU/image%201.png)
        
        上面的图是作者构造的例子，一共有12个数据点，$d=2$，并且子空间的大小$m=1$，可以看到在第一张图中每一个数据点所在位置以及构图。第二张图是每一个结点对应的聚簇号以及其与query（图中星星代表的点）对应子空间的距离。
        
2. **ANN Search**
    
    由于在一个批次中的所有查询都是独立的，它们就可以以一种过易并行的方式进行检索。在整体架构图中，每一个查询都会单独运行在一个individual CUDA thread-block中。
    
    搜索部分的逻辑如下（这个地方是我看图说话，它原文中没有像我一样解释）：
    
    1. 在搜索开始时，聚类重心点是起始搜索点，然后会将所有聚类重心点的邻居点从CPU中送到GPU中，在GPU中计算出所有聚类重心点的邻居与query point的距离。
    2. 接下来进入循环阶段，一共有A和B两个**CPU与GPU并行运行**的阶段和若干小步骤：
        1. 将下一个候选点从GPU传给CPU。
        2. 进入A：
            1. CPU获取一个候选点的邻居，并将候选点的邻居传给GPU
            2. GPU进行排序
            
            CPU与GPU并行这个过程，直到所有候选点的邻居被全部传送给GPU
            
        3. 进入B：
            1. GPU**过滤掉**A中传送过来的候选点的邻居中**已经访问过的结点。**
            2. GPU计算上一步剩余的结点与query点之间的距离。
            3. GPU挑选下一轮的候选结点。
3. **Re-ranking**：使用实际数据而非压缩数据对最终搜索结果进行排序，然后取最终最近的$k$个元素作为最终结果集合。

# 4. BANG: BILLION-SCALE ANN SEARCH ON A SINGLE GPU

## 4.1 Search Algorithm

算法1就是普通的图搜索算法

算法2描述了BANG的搜索方案，在Vamana graph上使用批次查询。由于同一个批次的所有查询是独立的，他们可以以一种过易并行的方式进行。每一个查询运行在一个独立的CUDA线程块中，这样就能保证最大化吞吐量。

![image.png](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-08-20-BANG Billion-Scale%20Approximate%20Nearest%20Neighbor%20Search%20using%20a%20Single%20GPU/image%202.png)

这里用语言描述一下算法2的流程：

1. 对于查询集合中的每一个查询向量，首先初始化一个起始向量$u_i^*$，并且接下来的流程对于每一个查询向量并行执行。（第二行）
2. 进入循环，如果所有在worklist中的点都已经被搜索过了，就退出循环到第11步，否则接下来的循环内容。（第四行）
3. 对于当前的$u_i^*$，找到它的所有邻居$N_i$。（第五行）
4. 将邻居集合$N_i$传给GPU。（第六行）
5. 用一个布隆过滤器**并行从$N_i$**过滤掉已经访问的点，形成$N_i'$。（第七行）
6. 对于在$N_i'$中的每一个点，并行计算每一个点与query的距离。（第十一行）
7. 用并行化归并排序将$N_i'$按照距离$D_i$从小到大进行排序。（第十三行）
8. 将排序后的$N_i'$与worklist的$L_i$进行合并。
9. 在GPU更新下一个要访问的结点，并将该点从GPU传输到CPU（15行）
10. 更新当前是否需要继续循环的标志位
11. 返回$k$个最近邻。

## 4.2 Construction of 𝑃𝑄𝐷𝑖𝑠𝑡𝑇𝑎𝑏𝑙𝑒 in Parallel

构建PQDistTable，用于搜索阶段的距离查找。PQDistTable包含每一个子空间中每一个聚簇向量（**centroid vector**）到query对应向量子空间的距离。这个PQDistTable会常驻在GPU中，直到对应的批次中的查询结束。

对于一个批次的查询中，PQDistTable负责提前计算并存储每个**centroid vector**到query对应向量子空间的距离。在这里，提到了使用ADC（Asymmetric Distance Computation，在之前的笔记中讲过）的方式在后续查询中计算查询点与数据点之间的距离。

作者维护使用一个大小为$\rho \cdot m \cdot 256$的连续数组PQDistTable，$\rho$是每一个query batch的大小，$m$是子空间的大小，$c=256$代表每一个子空间中聚簇的数量。作者从经验上，设置$m=74$，并且每一个线程块会处理一个查询，所以最后会有$\rho$个并行线程块在GPU中。

更多的，对于任意一个query，他的每一个子空间向量与聚簇的距离计算可以独立并行计算，并且注意到，这样的操作都是在一个线程中被完成的，这样可以达成每一个线程充足的工作负载并将每一个query的计算限制在一个线程块中。

构建整个PQDistTable的时间复杂度为$O((m \cdot subspace\_size)\cdot c \cdot \rho)$，原因在于，对于每一个查询（大小为$\rho$）中的每一个向量子空间（大小为$m$）中的每一个聚簇向量（大小为$c$），都要做距离计算（距离计算消耗$subspace\_size = \frac{d}{m}$）。由于并行化的缘故，算法的复杂度会消除掉与查询有关的复杂度，变为$O(m \cdot subspace\_size)$。

## 4.3 Handling Data Transfer Overheads

由于PCIe的带宽有限，在CPU与GPU之间传输内容是很耗时的相比起GPU惊人的处理速度，因此，BANG只会传输最少的必要的内容在CPU与GPU之间。特别地：

- 从GPU到CPU：
    - 会传输最终的搜索结果（AKNN最终结果 18行）
    - 会传输query的候选点（第十六行）
- 从CPU到GPU：
    - 将候选点的邻居传输过去（第6行）
    - 将非压缩的原始向量传送过去用于挑选最终结果。

更多地，在搜索进行中时，获取一个给定的结点的所有邻居，需要使用到所有的CPU线程。

利用图索引的结构，将结点与其邻居结点放置在相邻的CPU内存块中（这里仿照的是硬盘版本的局部性原理，硬盘那边是防止更多的I/O开销，这里也可以引申出这个道理）

## 4.4 Handling Visited Vertices: Bloom Filter

一个普通的布隆过滤器，只会有如下的两种形态：

- 对于每一个结点设置一个bit位记录其是否被访问过
    - 对于GPU-based的方法，这种方法对于十亿级别的数据以及上万的查询来说，需要难以承受的显存空间大小。故没有办法采用这种办法。
- 使用一个数据结构，如优先队列，哈希表之类的结构来记录每一个被标记为找到的结点。
    - 虽然这种办法只会存储被标记的点，故只会有较少的点被记录保存在数据结构中。然而，一个哈希表，优先队列，之类的结构是一个动态数据结构，对于使用GPU来维护它们是一种挑战，并且效率也不好。

鉴于之前所提到的各种挑战，作者采用了一种好的布隆过滤器的方法，这种方法是一种**近似集合成员测试**的方法，并且有最小的内存占用量，非常好的利用了GPU的并发性。

作者对于每一个query使用了一个拥有$z$个bool变量的数组，$z$是固定的，由于采用的是近似集合成员测试的方法，所以会允许一个名叫false-positive rate的小错误率存在（也就是预测错误的情况），并且采用多个哈希函数来做Bloom Filter。

这里多写一点哈希函数如何做Bloom Filter：

- **初始化Bloom Filter**: 创建一个大小为`m`的位数组（通常用0来初始化），长度`m`是根据期望的误报率和插入的元素数量确定的。
- **定义哈希函数**: 选择`k`个哈希函数。每个哈希函数将输入的元素映射到一个范围在`0`到`m-1`之间的整数。哈希函数需要尽可能独立，以避免冲突。
- **插入元素**:
    - 对于要插入的每个元素，依次通过这`k`个哈希函数进行哈希处理，得到`k`个哈希值。
    - 将这`k`个哈希值对应的位数组中的位置设为`1`。
- **检查元素是否存在**:
    - 对于要查询的元素，同样使用这`k`个哈希函数计算哈希值。
    - 检查这些哈希值对应的位数组位置是否都为`1`。如果所有位置都为`1`，则该元素**可能存在**；如果有任何一个位置为`0`，则该元素**一定不存在**。

## 4.5 Parallel Neighbor Distance Computation

- **非对称距离计算**：对于批处理中每个查询点，算法会在每次迭代中计算该查询点与当前邻居列表之间的非对称距离。这是通过使用之前构建的`PQDistTable`（产品量化距离表）数据结构来实现的，该表保存了质心与查询向量之间的距离。计算是并行进行的，每个查询点被分配到一个线程块中。
- **并行化策略**：给每一个query独立分配一个线程块。
- **汇总策略**：~~当使用GPU-based的方法进行最近邻检索计算，汇总各个并行单元的结果并使用英伟达的库是很常见的。这种策略可能会在块级别，warp-level（SIMT：单指令多线程）发挥作用。作者采用了一种独特的方式，其将一个线程块拆成了$g$个组，每个组有$g_{size}=\frac{t_b}{g}$个线程。这些组合作计算**一个**query距离**一个**邻居的距离。~~
    
    ~~所以对于计算一个query到$m$个聚类中心的距离来说，每个线程组可以被分配到$\frac{m}{g_{size}}$，每一个线程**顺序**地计算在一个段中的总和值使用一个局部线程的寄存器，避免同步操作。~~
    
    - **线程块划分**：
        - 传统方法通常使用Nvidia的CUB库进行归约操作，通常在线程块或warp级别进行。而在这个新方法中，一个线程块（大小为$t_b$，即线程总数）被划分为多个较小的组。
        - 这个线程块被分为$g$个组，每个组有$g_{size}$个线程，其中$g_{size}=\frac{t_b}{g}$。这些组共同计算查询点与邻居之间的距离。
    - **距离计算**：
        - 每个线程组（大小为$g_{size}$）负责累加$m$个质心（centroid）的部分距离。这些质心被分成多个段，每个段的大小为$\frac{m}{g_{size}}$。（这个地方有点疑惑，这个式子代表的意思是对于一个线程块，其中任意一个线程需要计算多少个聚类簇。感觉好像这样分配不合理？怎么会特别针对某个线程块做这种运算？没看懂）
        - 每个线程在它们各自负责的段内顺序地计算数值之和，使用的是线程本地寄存器。这意味着每个线程在计算过程中不需要与其他线程进行同步，从而避免了同步带来的开销。
    - **归约策略**：
        - 当线程计算出每个段的和后，会采用两种策略之一来聚合这些结果：
            1. **原子操作（Atomics）**：使用`atomicAdd()`函数将结果累加到最终的和中。
            2. **子warp级别归约**：使用CUB的`WarpReduce`在部分线程（子warp）内进行归约操作。
        - 结果显示，第二种方法（子warp级别归约）比第一种原子操作略微更好，从而在实验中达到了最佳性能。
    - **性能考量**：
        - 这种分段方法比传统的CUB的`WarpReduce`和`BlockReduce`等标准方法更有效，尤其是在处理大规模数据集时，此内核占据了显著的运行时间。

## 4.6  Prefetching Candidate Nodes

- **候选节点的传统传送方式**：
    - 在GPU更新工作列表 $L_i$ 后，候选节点 $u_i^∗$ 会被传送给CPU。
    - 在此期间，GPU会处于空闲状态，等待CPU完成邻居节点的收集并返回结果。
- **优化策略：提前预测候选节点**：
    - 优化策略在计算邻居距离后，立即预测下一个候选节点，而不是等待邻居列表排序和工作列表更新完成后再传送。
    - 这一优化选择新邻居列表 $N_i'$ 中距离最近的节点和工作列表 $L_i$ 中未访问节点的首个节点，择优选择最佳节点作为候选节点。
- **并行任务执行**：
    - GPU在继续执行当前迭代的任务（如邻居列表排序和工作列表更新）的同时，将提前选定的候选节点传送给CPU。
    - CPU则并行执行邻居节点的收集工作，从而在GPU完成任务的同时，完成邻居节点的收集，并立即返回给GPU。
- **性能提升**：
    - 这种优化策略有效地消除了GPU的空闲时间，实验结果显示，吞吐量提高了大约10%。
- **具体实施步骤**：
    - 在算法第13行之前插入提前选择候选节点的操作。
    - 将提前选择的候选节点发送至CPU，并行执行后续任务。

## 4.7 Parallel Merge Sort

## 4.8 Parallel Merge

- 作者采用了并行归并排序的方式来做排序操作
    - 传统的归并排序方法是：
        
        **算法通常采用自底向上的方式，逐步合并小的已排序子列表来形成更大的已排序列表。具体来说，一个线程通常负责将两个已排序的子列表合并成一个更大的已排序列表**。在这种方法中，随着排序进程的推进，列表的数量减少，而每个线程需要处理的工作量增加。这会导致两个问题：
        
        1. **并行度降低**：由于待合并的列表数量减少，能够并行工作的线程数量也会减少，导致GPU的利用率降低。
        2. **线程空闲**：在合并较大的列表时，排序过程中的大部分线程会处于闲置状态，等待其他线程完成其合并任务。
    - 本文采用的归并排序的方法：
        - **小型列表的并行排序**：BANG方法将待排序的邻居节点列表划分为较小的子列表，通常最多包含64个邻居节点。**每个子列表由一个线程块来处理，每个线程负责一个邻居节点的排序**。这样可以充分利用GPU的并行计算能力，因为所有的线程都在同时工作，并且列表较小，减少了线程等待的情况。
        - **并行合并例程**：在BANG中，两个已排序的列表并不是由单个线程合并，而是通过并行的方式，多个线程同时处理不同部分的列表。例如，假设有两个已排序的列表A和B，BANG方法会为每个元素分配一个线程，利用二分查找的方式确定该元素在合并后的列表中的位置。这样，每个线程都能独立地工作，最大化GPU的利用率。
            
            ![image.png](https://raw.githubusercontent.com/xwt1/xwt1.github.io/main/_misc/picture/2024-08-20-BANG Billion-Scale%20Approximate%20Nearest%20Neighbor%20Search%20using%20a%20Single%20GPU/image%203.png)
            
        - 总的来说：一次排序的总时间复杂度计算如下：
            - 大体流程和传统方法一致，采取自底向上的方法，故可写出时间复杂度递推式：$T(n)=2*T(\frac{n}{2})+nlog(n) = O(nlog^2n)$，后半部分为合并的复杂度，前半部分为合并次数，$n$是排序数组的长度。
            - 合并如上面的图所示，每个线程负责一个元素，若当前需要合并$L_{in1}$和$L_{in2}$两个列表，每一个线程会查找自身负责的那个元素在另一个列表中的位置，如图中$T_1$负责的值为$3$的元素在$L_{in1}$的1号位置，那么它就需要利用二分查找在$L_{in2}$中找到对应位置，为1号位置，则其最终在合并后的列表中就应该处于$1+1=2$号位置。
            - 由于待排序的元素实际是候选邻居的数量（成为候选点），故$n=NumNbrs$，$NumNbrs$是每一轮排序中候选邻居的数量，一般和最终的$k$（$k$个最近邻）是同一个数量级的，故复杂度为$O(NumNbrs \cdot log^2(NumNbrs) \cdot \rho)$，而由于并行化的原因，每一个查询（$\rho$）都会分配一个线程块，每一个候选邻居结点（$NumNbrs$）都会安排一个线程块中的一个线程，故这两项都可以省略，最终复杂度为$O(log^2(NumNbrs))$，比起传统复杂度要快很多$O(NumNbrs \cdot log(NumNbrs))$。

## 4.9 Re-ranking

BANG系统引入了最终的重新排序步骤以提高整体的召回率。在搜索过程中，使用的是近似距离（通过PQ距离计算得到），而在搜索收敛之后，会进行重新排序步骤，以确保最终结果的准确性。

重新排序步骤具体包括以下内容：

1. **计算精确的L2距离**：对于每个查询向量和其对应的候选节点，重新计算它们之间的精确L2距离。这个步骤是并行执行的，每个查询点的距离计算是独立的。
2. **排序候选节点**：根据计算出的精确距离，对候选节点进行排序，最终报告距离查询向量最近的前k个节点作为最终的近邻。
3. **性能提升**：通过重新排序步骤，实验表明，BANG系统在某些数据集上的召回率提高了10-15%。这表明重新排序在处理高维数据时对提高检索精度非常有效。

工作跨度分析显示，重新排序的总工作量为$O((d × |C| + |C| × log²(|C|)) × ρ)$，其中$|C|$是每个查询的候选节点数量，$d$是向量维度。通过并行化方案，重新排序的跨度为$O(d + log²(|C|))$。

# 5 IMPLEMENTATION

除了之前处理的大型数据集的BANG版本，作者还提出了针对小型和中型数据集的优化版本，这些数据集可以完全存储在GPU内存中。

## 5.1  In-memory Version

在这个版本中，当图结构可以完全放入GPU内存时，BANG系统的CPU与GPU之间的数据传输可以显著减少，从而提高整体吞吐量。具体优化包括：

- **消除CPU-GPU通信开销**：将整个图结构保存在GPU内存中，消除在搜索过程中CPU从主内存中获取邻居节点并将其传输给GPU的步骤。
- **提高吞吐量**：通过减少PCIe总线上的数据传输开销，吞吐量得到了显著提升，实验中提升幅度可达50%。

## 5.2 Exact-distance Version

在这个版本中，BANG进一步优化了In-memory版本，直接使用精确的L2距离计算，而不是在搜索过程中使用PQ距离。具体优化包括：

- **移除近似距离计算**：在搜索过程中直接计算精确的L2距离，省去了对压缩向量的近似距离计算以及后续的重新排序步骤。
- **提升精度**：通过直接使用精确距离计算，消除了因压缩带来的距离误差，使得检索结果更加精确。

# 6. EXPERIMENTAL SETUP

实验部分自然是大吹特吹，暂且不提。