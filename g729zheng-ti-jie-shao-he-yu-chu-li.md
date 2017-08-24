# 整体介绍

本文处理的信号是由电话带宽的模拟信号经过G.712推荐的带通滤波，然后经过8Khz采样，并将样值做16bit PCM量化得到的。  
CS-ACELP编码基于码激励线性预测（CELP）编码模型。编码器处理单位为8khz采样率下10ms共80个样点的数据。对每一个10ms的数据帧进行分析，可以得到CELP模型的参数（线性预测滤波器系数、自适应码书和固定码书索引和增益），这些参数经过编码后进行传输。编码结果的Bit分配如下图所示。在解码端，这些参数被用于检索激励和合成滤波器参数。

![G729编码bit分配](http://note.youdao.com/yws/public/resource/fedd85d892cd55c62ed9d238f41e4155/xmlnote/A992B1BAEBD549A38E9ACF958D2A7E2A/103744)

解码段通过一个基于10阶线性预测滤波器的短时合成滤波器对激励进行滤波；长时滤波器或称为基音合成滤波器通过自适应码书方式实现。重建得到的语音信号，最后还要经过一个后置滤波器的增强。

语音信号在预处理模块经过高通滤波和量化然后输出，作为后续其他模块分析的输入信号。LP分析每10ms一帧进行一次计算LP系数，然后将这些LP系数抓换为LSP系数，转换后对LSP系数通过一个两段矢量量化器量化为18bit。激励信号的选择是通过一个分析合成搜索过程，通过使原始语音信号和合成语音信号的感知加权误差最小化得到的。实际实现是通过将误差信号通过一个有LP系数得到的感知加权滤波器进行滤波完成的，感知加权滤波器是自适应的，目的是为了输入语音信号有一个平坦的频率响应。

激励参数（固定码书参数和自适应码书参数）是通过5ms子帧进行计算的的。第1个子帧使用内插的LP滤波器参数包括量化和非量化的，第2个子帧使用计算得到的LP滤波器参数（包括量化和非量化的）。前面利用感知家全国的语音信号10ms子帧得到开环基音延迟的一个估计值。

然后对每一个子帧执行如下步骤的处理。通过加权的合成滤波器`$\tfrac{W(z)}{\hat{A(z)}}$`对LP估计残差进行滤波，滤波器的内部状态通过对LP残差和激励信号进行滤波来更新的。这种方式等价于传统的从加权语音信号中减去加权合成滤波器零输入响应的方法。接下来计加权合成滤波器的冲击响应h\(n\)计算。  
然后是搜索自适应码书延迟和增益的闭环基因分析，该过程是利用目标信号x\(n\)及h\(n\)，在开环基音分析结果附近进行搜索完成的，搜索采用了1/3份数倍的精度。对第一个子帧，基音延迟采用8bit编码，第二个子帧使用5bit编码和第一个子帧的差值。目标信号x\(n\)减去自适应码书的对合成结果的贡献，作为固定码书搜索的目标信号，搜索固定码书表示的最佳激励。  
一个基于代数码书的17bit编码用来表示固定码书搜索的结果，自适应码书和固定码书的增益使用7bit进行适量量化统一编码（固定码书的增益计算采用的MA预测）。  
最后，用确定的激励信号更新滤波器的存储。

本文将详细介绍下图中各个模块的功能和实现方式。  
![G729编码框图](http://note.youdao.com/yws/public/resource/fedd85d892cd55c62ed9d238f41e4155/xmlnote/16BFEB69E5494E2A884D73578CE92346/103736)

# 基本概念

## DPF

DPF \(Double Precision Format \)  
A是一个正常的32位整数，DPF的表示方式会将A用2个16位有符号数表示，他们关系如下。

```cpp
A = AH<<16 + AL<<1
```

由于AL当中包含了移位符号位，所以他们和标准32位整数的表示方式是不同的。使用DPF的原因是可以方便做快速的乘法和除法计算。

## DPF的运算

## 综合滤波

Syn\_filt\(Aq, ai\_zero, h1, L\_SUBFR, zero, 0\);  
将输入信号做`$\frac{1}{A(z)}$`综合滤波。

* Aq，Q12，\(m+1\)阶预测滤波器A\(z\)的系数
* ai\_zero，输入信号
* h1，输出信号
* L\_SUBFR，滤波长度
* zero，滤波器初始状态（如果更新滤波器状态，则调用返回保存结束状态）
* 0，不更新滤波器状态

## 预测滤波误差

```c
void Residu(
  Word16 a[],    /* (i) Q12 : prediction coefficients                     */
  Word16 x[],    /* (i)     : speech (values x[-m..-1] are needed         */
  Word16 y[],    /* (o)     : residual signal                             */
  Word16 lg      /* (i)     : size of filtering                           */
)
```

# 代码整体介绍

```
graph TD
    A[打开wav文件] --> B[初始化]
    B --> G[读一帧数据]
    G --> C[预处理]
    C --> D[ld8k编码]
    D --> E[bit流->Byte流]
    E --> F[字节流输出到文件]
    F --> |继续读取| G
```

如果开启了时间对齐，则会先预读40个到从speech开始的内存当中，作为上一帧的内容，并进行单独的预处理过程。

# 预处理

二阶高通滤波器，截止频率为140Hz，滤除低频噪声。

```cpp
#ifdef SYNC
  /* Read L_NEXT first speech data */
  fread(&new_speech[-L_NEXT], sizeof(Word16), L_NEXT, f_speech);
  Pre_Process(&new_speech[-L_NEXT], L_NEXT);
#endif
```

充传输函数得到滤波计算公式：

```math
y[i]=\tfrac{1}{2}b_0*x[i]+\tfrac{1}{2}b_1*x[i-1]+\tfrac{1}{2}b_2*x[i-2]+a_1*y[i-1]+a_2*y[i-2]
```

`Coder_ld8k()`函数，处理的数据来源是读取到全局数组的原始语音信号，长度为80字节。输出为编码结果和合成的语音信号（用于调试）。

# 开环基音分析

## 计算开环延迟

```
graph LR
A(防止计算溢出) --> B(分阶段寻找基因时延)
B --> C(基因延迟选定)
```

对于不同区间计算得到的延迟机器对应的自相关结果，并不是直接对比自相关的大小选取延迟，而是以0.85倍递减。

### 防止计算溢出

```cpp
// Pitch_ol()
  /*--------------------------------------------------------*
   *  Verification for risk of overflow.                    *
   *--------------------------------------------------------*/

   Overflow = 0;
   t0 = 0;

   for(i= -pit_max; i< L_frame; i++)
     t0 = L_mac(t0, signal[i], signal[i]);

  /*--------------------------------------------------------*
   * Scaling of input signal.                               *
   *                                                        *
   *   if Overflow        -> scal_sig[i] = signal[i]>>3     *
   *   else if t0 < 2^20  -> scal_sig[i] = signal[i]<<3     *
   *   else               -> scal_sig[i] = signal[i]        *
   *--------------------------------------------------------*/
   if(Overflow == 1)
   {
     for(i=-pit_max; i<L_frame; i++)
       scal_sig[i] = shr(signal[i], 3);
   }
   else 
   {
     L_temp = L_sub(t0, (Word32)1048576L); // 1048576 = 2^20，源代码注释错误
     if ( L_temp < (Word32)0 )  /* if (t0 < 2^20) */
     {
        for(i=-pit_max; i<L_frame; i++)
          scal_sig[i] = shl(signal[i], 3);
     }
     else
     {
       for(i=-pit_max; i<L_frame; i++)
         scal_sig[i] = signal[i];
     }
   }
```

### 分阶段寻找基音延迟

大的延迟区间到小的延迟区间，区间级差乘以0.85递减。

```
graph TD
    A[T_op=t1, R=R_1] --> B{R_2 >= 0.85R ?}
    B -->|yes| C[R=R_2, T_op=t2]
    B -->|No| D{R_3 >= 0.85R ?}
    C --> D
    D --> |yes| E[R=R_3, T_op=t3]
    D --> |No| F[return T_op, R]
    E --> F
```

```cpp
  /*--------------------------------------------------------------------*
   *  The pitch lag search is divided in three sections.                *
   *  Each section cannot have a pitch multiple.                        *
   *  We find a maximum for each section.                               *
   *  We compare the maximum of each section by favoring small lag.     *
   *                                                                    *
   *  First section:  lag delay = pit_max     downto 4*pit_min          *
   *  Second section: lag delay = 4*pit_min-1 downto 2*pit_min          *
   *  Third section:  lag delay = 2*pit_min-1 downto pit_min            *
   *--------------------------------------------------------------------*/
   // 阶段1：[pit_max, 4*pit_min]
   j = shl(pit_min, 2);
   p_max1 = Lag_max(scal_sig, L_frame, pit_max, j, &max1);
   // 阶段2：(4*pit_min, 2*pit_min]
   i = sub(j, 1); j = shl(pit_min, 1);
   p_max2 = Lag_max(scal_sig, L_frame, i, j, &max2);
   // 阶段2：(2*pit_min, pit_min]
   i = sub(j, 1);
   p_max3 = Lag_max(scal_sig, L_frame, i, pit_min , &max3);

  /*--------------------------------------------------------------------*
   * Compare the 3 sections maximum, and favor small lag.               *
   *--------------------------------------------------------------------*/
  if( sub(mult(max1, THRESHPIT), max2)  < 0) // Threshold to favor small pitch; 0.85,Q15
  {
    max1 = max2;        // max2 较大，代表选定的lag相关性较好
    p_max1 = p_max2;    // 只要符合条件，就倾向选择更小的基音延迟
  }
  if( sub(mult(max1, THRESHPIT), max3)  < 0)
  {
    p_max1 = p_max3;
  }
```

最终，将p\_max1，即基因延迟返回。

### 计算基音延迟

```cpp
// Lag_max()
  max = MIN_32; 

  p_max = lag_max;
  for (i = lag_max; i >= lag_min; i--)
  {
    p  = signal;
    p1 = &signal[-i];
    t0 = 0;

    for (j=0; j<L_frame; j++, p++, p1++) // 计算
      t0 = L_mac(t0, *p, *p1);

    L_temp = L_sub(t0,max);
    if (L_temp >= 0L)   // 寻找最大的自相关对应的延迟
    {
      max    = t0;
      p_max = i;
    }
  }

  /* compute energy */
  t0 = 0;
  p = &signal[-p_max];
  for(i=0; i<L_frame; i++, p++)
    t0 = L_mac(t0, *p, *p);

  /* 1/sqrt(energy),    result in Q30 */
  t0 = Inv_sqrt(t0);

  /* max = max/sqrt(energy)                   */
  /* This result will always be on 16 bits !! */

  L_Extract(max, &max_h, &max_l);
  L_Extract(t0, &ener_h, &ener_l);

  t0 = Mpy_32(max_h, max_l, ener_h, ener_l);    // 实际为 max / sqrt(t0)
  *cor_max = extract_l(t0);     // 最大的自相关值
```

# 自适应码书搜索

自适应码书：声音激励模型将可能的各种激励预先存储在存储器内，通过某种判据决定哪一种激励是当前信号的最佳激励，并且把最佳激励的存储地址作为激励的表征，例如码激励模型或矢量激励模型等。存储器内容随时间变化的部分称为自适应码书。**自适应码书的搜索等价于基音检测**。

相比于g723，可以看出g729对自适应激励作了化简，当基音周期较短时，引入分数基音延迟，这个对应g723的五阶加权闭环基音预测。当基音周期较大时，则简化了处理。当取自适应激励时，我们仍然看到做了升抽样，这在一定程度上，仍然是对应g723的五阶加权。  
闭环基音分析每个子帧（40个符号）进行一次。

## 闭环基音分析

在整数基音周期附近\(左右3个单位\)，搜索分数基音周期，解析度为1/3。

```
graph LR
A(上采样) --> B(低通滤波)
```

首先构造出感知加权滤波器`$W(z)$`和综合滤波器`$A(z)$`，并将它们串联。

计算知觉加权滤波器的LPC系数。

```cpp
/*---------------------------------------------------------------*
 * Find the weighted LPC coefficients for the weighting filter.  *
 *---------------------------------------------------------------*/
Weight_Az(A, gamma1[i_gamma], M, Ap1);
Weight_Az(A, gamma2[i_gamma], M, Ap2);
i_gamma = add(i_gamma,1);
```

计算感知加权滤波器的冲击响应。  
用单位冲激响应去激励这样的一个系统`$\frac{1}{Aq(z)}*\frac{Ap1(z)}{Ap2(z)}$`，即将量化后的预测系数A与感知加权滤波器做一个串联，得到一个冲激响应。

```cpp
/*---------------------------------------------------------------*
 * Compute impulse response, h1[], of weighted synthesis filter  *
 *---------------------------------------------------------------*/
for (i = 0; i <= M; i++) 
{
    ai_zero[i] = Ap1[i];
}
```

这里得到完整的冲激响应，加入感知加权滤波，用未量化的Az系数生成感知加权，然后用经由量化的Az系数滤波。

```
Syn_filt(Aq, ai_zero, h1, L_SUBFR, zero, 0);
Syn_filt(Ap2, h1, h1, L_SUBFR, zero, 0);
```

根据这个冲激响应,得到残差信号\(激励信号\)，并推算出扣除零输入响应的目标语音信号xn。

## 激励信号

```cpp
//
    *              |------|  res[n]                                          *
    *  speech[n]---| A(z) |--------                                          *
    *              |------|       |   |--------| error[n]  |------|          *
    *                    zero -- (-)--| 1/A(z) |-----------| W(z) |-- target *
    *                    exc          |--------|           |------|          *
```

计算目标向量，而没有从感知加权滤波器中减去零输入响应。这种配置在使用定点实现时有更好的性能。  
`$\frac{1}{A(z)}$`的存储使用`$res[n]-exec[n]$`通过`$\frac{1}{A(z)}$`滤波后的结果更新，或者简单的使用输入语音减去综合的语音：

```math
error[n] = speech[n] - syn[n]
```

`$W(z)$`的存储使用`$error[n]$`通过`$W(z)$`滤波后的结果来更新，或简单使用目标矢量减去滤波后的自适应和固定码矢量：

```math
target[n] - gain\_pit*y1[n] - gain\_code*y2[n]
```

这里计算目标向量，有可能做了去除零输入响应之类的操作，先用包含零输入响应的信号滤波\(量化后的Az系数\)，得到残差信号，再用残差信号还原，就得到去除零输入响应的语音信息，对语音信号进行感知加权\(采用的是未量化的Az系数\)。

```cpp
Residu(Aq, &speech[i_subfr], &exc[i_subfr], L_SUBFR);   /* LPC residual ，保存在exc当中*/
Syn_filt(Aq, &exc[i_subfr], error, L_SUBFR, mem_err, 0); // exc 经过 1/Aq(z)滤波
Residu(Ap1, error, xn, L_SUBFR);
Syn_filt(Ap2, xn, xn, L_SUBFR, mem_w0, 0);    /* target signal xn[]*/
```

### 闭环即音延迟搜索

首先，限定第一个子帧基音延迟搜索的范围。这个在搜索到开环基音延迟后就进行了。

```cpp
/* Range for closed loop pitch search in 1st subframe */

T0_min = sub(T_op, 3);  // T0_min = T_op - 3;  T_op开环基音延迟搜索结果
if (sub(T0_min,PIT_MIN)<0)  // PIT_MIN最小基音延迟
{
    T0_min = PIT_MIN;
}

T0_max = add(T0_min, 6);     // T0_max = Tmin + 6
if (sub(T0_max ,PIT_MAX)>0)  // 如果max越界，则重新设定max和min
{
 T0_max = PIT_MAX;
 T0_min = sub(T0_max, 6);
}
```

TO\_min TO\_max 是开环基音周期TO-3、Tmin+6得到的，Pitch\_fr3 在这个范围内进一步搜索自相关最大的基音延迟。

```cpp
    /*----------------------------------------------------------------------*
     *                 Closed-loop fractional pitch search                  *
     *----------------------------------------------------------------------*/
    T0 = Pitch_fr3(&exc[i_subfr], xn, h1, L_SUBFR, T0_min, T0_max,
                               i_subfr, &T0_frac);
    index = Enc_lag3(T0, T0_frac, &T0_min, &T0_max,PIT_MIN,PIT_MAX,i_subfr);    // 基音周期参数打包

    *ana++ = index; // 对P1高6位进行奇偶校验
    if (i_subfr == 0) 
    {
      *ana++ = Parity_Pitch(index);
    }
```

Pitch\_fr3负责闭环基音周期搜索,如必要再进一步做更精细的分数基音延迟搜索。  
\(当基音周基小于84就需要搜索分数延迟的基音周期，保证在基音周期较短的情况下，得到的自适应激励更精确一些\)

## 闭环搜索

与开环用信号的自相关最大作为搜索依据不同，闭环搜索的依据是自适应激励与h1卷积与目标信号的相关最大\(这点与g723的五阶闭环基音周期搜索类似\)。

x\[n\]是目标语音信号  yk\[n\]是自适应激励与h1卷积  自然yk\[n\]应该有九组  
而我们就是使找出corr\[k\]最大的那个yk\[n\]  对应的,也就是找到了最佳的自适应激励

### 节省卷积计算量

Norm\_Corr中有一小段代码,与g723的类似,为了节省计算量,  
由于每个循环我们都在不断更新yk\[n\],而这个可以根据卷积的性质做一些化简,代码片段如下:

```cpp
    k=sub(k,1);//lsc 这里往 t_max方向移动一格,然后再更新卷积,参照723的作法,节省一些运算量
    for (j = L_subfr-(Word16)1; j > 0; j--)
    {
        s = L_mult(exc[k], h[j]);
        s = L_shl(s, h_fac);             /* h is in Q(12-scaling) */
        s_excf[j] = add(extract_h(s), s_excf[j-1]);
    }
    s_excf[0] = shr(exc[k], scaling);
```

这个类似的技巧在分析g723时分析过了,这里不再详细说明,读者只需要注意一下即可。

### 分数倍基音延迟搜索

sinc函数\(broadcom公司的标志\)，升抽样就是借助原始信号与sinc的卷积来完成的。原理很简单：采样的过程\(假设满足香农定理\),那么采样后信号的频域曲线为原始信号频域函数周期重复，只要在频域上做一个低通滤波即可还原。对应时域的处理就是原始信号与sinc函数卷积。

然后就是分数基音延迟搜索了,因为搜索涉及到升抽样,  
g729采用了一些取巧的做法,即直接对互相关的结果corr\(k\)这个序列进行升抽样.  
Interpol\_3 这个函数就是负责升抽样的,原来在上两节讲过了,就是将序列与sinc函数进行卷积。  
应该注意 corr\(k\)这个序列因为要进行升抽样,而间隔要被拉开,所以对应的inter\_3在循环跳变时的间隔为UP\_SAMP而本身sinc函数是关于y轴对称的，以及一些非因果系统的原因，卷积的过程会比较怪异，但它终究只是在计算卷积。  
小技巧就不分析了

```cpp
Word16 Interpol_3(      /* (o)  : interpolated value  */
  Word16 *x,            /* (i)  : input vector        */
  Word16 frac           /* (i)  : fraction            */
)
{//lsc 这个是升抽样, 通信信息与sinc函数卷积来达到升抽样的目的
//lsc 抽样的过程,在频域的表现,是将原始信号的频域按抽样频率复制多份,要还原信号,
//lsc 我们可以通过在频域上还原信号,即把抽样后的信号做一个低通滤波,去掉重复的高频成份,就达到还原的目的了
//lsc 这段代码就是在做这些事情,而频域上的矩形滤波器对应时域的信号就是sinc
//lsc 在时域上,只要将原始信号与sinc进行卷积,就完成了升抽样的工作
```

到这里已经得到基音延迟了。

## 得到自适应激励

至此基音延迟正数k和分数t已经确定，自适应码矢量v\(n\)可以通过内插历史激励信号u\(n\)得到。

获取自适应激励，这里也涉及到升抽样，原理与Interpol\_3一样。Pred\_lt\_3就是对历史激励源进行升抽样,然后取出一组激励作为自适应激励源。

gain\_pit 其实是自适应码本增益的一个估值\(自适应码本的增益会在后继章节接着讲\)，利用这个估值,得到固定码本的搜索目标信号

G\_pitch 不但完成了这个估值计算,还计算了一些后继增益量化时会用到的一些项  
保存在g\_coeff,这些会在后继章节分析\(大体也是求偏导之类的\)

有了这个估值 gain\_pit,就可以从目标语音信号xn中扣除自适应激励成份,而得到固定码本的目标向量,

就是做一个卷积,然后从xn中扣除目标向量保存在xn2数组里头.

```cpp
//-- 
Pred_lt_3(&exc[i_subfr], T0, T0_frac, L_SUBFR); // 对最佳激励进行升抽样
Convolve(&exc[i_subfr], h1, y1, L_SUBFR);       // 升抽样后的历史激励与冲激响应卷积
gain_pit = G_pitch(xn, y1, g_coeff, L_SUBFR);

/* clip pitch gain if taming is necessary */
temp = test_err(T0, T0_frac);

if( temp == 1){
  if (sub(gain_pit, GPCLIP) > 0) {
    gain_pit = GPCLIP;
  }
}

/* xn2[i]   = xn[i] - y1[i] * gain_pit  */
// 扣除自适应激励,搜索固定码本
for (i = 0; i < L_SUBFR; i++)
{
  L_temp = L_mult(y1[i], gain_pit);
  L_temp = L_shl(L_temp, 1);               /* gain_pit in Q14 */
  xn2[i] = sub(xn[i], extract_h(L_temp));
}
```

# 参考文献

1. [g729源码分析-2-共振锋感知加权](http://blog.csdn.net/lsccsl/article/details/7449361)
2. [深蓝怒火的CSDN博客](http://blog.csdn.net/lsccsl/article/list/3)



