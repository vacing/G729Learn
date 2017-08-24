# 代数码书搜索（固定码书搜索）
经过自适应码本搜索后，语音信号的时域相关性被极大地去除了，剩下的残差信号接近于随机信号。
残差信号是与一个伪随机的固定码本进行"匹配"，算出最佳的码本索引与增益。

```cpp
// ---
   /*-----------------------------------------------------*
    * - Innovative codebook search.                       *
    *-----------------------------------------------------*/

    index = ACELP_Codebook(xn2, h1, T0, sharp, i_subfr, code, y2, &i);
    *ana++ = index;        /* Positions index */
    *ana++ = i;            /* Signs index     */
```

## 固定码书搜索
固定码本搜索与g723的固定码本搜索算法几乎是一模一样的。g729与g723另一个不同的是，增益的估值与量化，这些笔者将在后继的章节分析。   
ACELP_Codebook和g723的ACELP_LBC_code函数算法是一模一样的。
```cpp
Word16  ACELP_Codebook(  /* (o)     :index of pulses positions    */ // 激励信号索引
  Word16 x[],            /* (i)     :Target vector                */ // 目标向量
  Word16 h[],            /* (i) Q12 :Impulse response of filters  */ // 滤波器冲击响应
  Word16 T0,             /* (i)     :Pitch lag                    */ // 基音延迟
  Word16 pitch_sharp,    /* (i) Q14 :Last quantized pitch gain    */ // 上一次量化的基音增益
  Word16 i_subfr,        /* (i)     :Indicator of 1st subframe,   */ // 是否是第一个子帧
  Word16 code[],         /* (o) Q13 :Innovative codebook          */ // 新的码书
  Word16 y[],            /* (o) Q12 :Filtered innovative codebook */ // 滤波后的新码书
  Word16 *sign           /* (o)     :Signs of 4 pulses            */ // 冲击的符号
)
```

```cpp
//--
  // 如果基音延迟T小于子帧长度，则将固定增益基音贡献算入冲击响应，并计算自相关
 /*-----------------------------------------------------------------*
  * Include fixed-gain pitch contribution into impulse resp. h[]    *
  * Find correlations of h[] needed for the codebook search.        *
  *-----------------------------------------------------------------*/
  sharp = shl(pitch_sharp, 1);          /* From Q14 to Q15 */
  if (sub(T0, L_SUBFR)<0) // 前置预滤波
     for (i = T0; i < L_SUBFR; i++){    /* h[i] += pitch_sharp*h[i-T0] */
       h[i] = add(h[i], mult(h[i-T0], sharp));
     }   
  Cor_h(h, rr); // 计算自相关

  // 计算目标矢量和冲击响应的互相关
 /*-----------------------------------------------------------------*
  * Compute correlation of target vector with impulse response.     *
  *-----------------------------------------------------------------*/
  Cor_h_X(h, x, Dn);

  // 查找码书
 /*-----------------------------------------------------------------*
  * Find innovative codebook.                                       *
  *-----------------------------------------------------------------*/
  index = D4i40_17(Dn, rr, h, code, y, sign, i_subfr);

  // 将固定增益基音贡献算入码矢量
 /*-----------------------------------------------------------------*
  * Compute innovation vector gain.                                 *
  * Include fixed-gain pitch contribution into code[].              *
  *-----------------------------------------------------------------*/
  if(sub(T0 ,L_SUBFR) <0)
     for (i = T0; i < L_SUBFR; i++) {  /* code[i] += pitch_sharp*code[i-T0] */
       code[i] = add(code[i], mult(code[i-T0], sharp));
     }

  return index; // 返回码书索引值
```

## 代数码书搜索

需要计算的公式如下：
```math
H=\begin{bmatrix}
h(0) & 0    & 0     & ... & 0           \\
h(1) & h(0) & 0     & ... & 0           \\
h(2) & h(1) & h(0)  & ... & 0           \\
...  & ...  & ...   & ... & ...         \\
h(39)& h(38)& h(37) & ... & h(0)        
\end{bmatrix}

H^T=\begin{bmatrix}
h(0)    & h(1)  & h(2)  & ... & h(39)           \\
0       & h(0)  & h(1)  & ... & h(38)           \\
0       & 0     & h(0)  & ... & h(37)           \\
...     & ...   & ...   & ... & ...         \\
0       & 0     & 0     & ... & h(0)        
\end{bmatrix}

x=\begin{bmatrix}x(39) & x(38) & x(37) & ... & x(0)\end{bmatrix}

d=xH^T
```
```math
d(i)=\sum_{n=i}^{39}x\prime(n)h(n-i)
```
```math
C=|d(m_0)|+|d(m_1)|+|d(m_2)|+|d(m_3)|
```

```math
\phi=H^TH 

\phi\prime(i,j)=\frac{1}{2}*sign[d(i)]*sign[d(j)]*\phi(i,j)=\frac{1}{2}*sign[d(i)]*sign[d(j)]*\sum_{n=j}^{39}h(n-i)h(n-j)
```
```math
\frac{E}{2}=\phi\prime(m_0, m_0)
+\phi\prime(m_1, m_1)+\phi\prime(m_0, m_1)

+\phi\prime(m_2, m_2)+\phi\prime(m_0, m_2)+\phi\prime(m_1, m_2)

+\phi\prime(m_3, m_3)+\phi\prime(m_0, m_3)+\phi\prime(m_1, m_3)+\phi\prime(m_2, m_3)
```

```cpp
static Word16 D4i40_17( /* (o)    : Index of pulses positions.               */ // 冲击索引
  Word16 Dn[],          /* (i)    : Correlations between h[] and Xn[].       */ // 冲击响应和目标矢量户相关
  Word16 rr[],          /* (i)    : Correlations of impulse response h[].    */ // 冲击响应自相关
  Word16 h[],           /* (i) Q12: Impulse response of filters.             */ // 滤波器冲击响应
  Word16 cod[],         /* (o) Q13: Selected algebraic codeword.             */ // 选定的代数码书码字
  Word16 y[],           /* (o) Q12: Filtered algebraic codeword.             */ // 滤波后代数码书码字
  Word16 *sign,         /* (o)    : Signs of 4 pulses.                       */ // 冲击的符号
  Word16 i_subfr        /* (i)    : subframe flag                            */ // 子帧标志
)
```

### 计算目标矢量和滤波器的互相关`$d(i)$`



```cpp
// acelp_co.c -> Cor_h_X()
   for (i = 0; i < L_SUBFR; i++)
   {
     s = 0;
     for (j = i; j <  L_SUBFR; j++)
       s = L_mac(s, X[j], h[j-i]);  // x(i)h(0)+x(i+1)h(1)+...

     y32[i] = s;        // 互相关值

     s = L_abs(s);
     L_temp =L_sub(s,max);
     if(L_temp>0L) {
        max = s;        // 得到max，做归一化使用
     }
   }
```



### 计算目标矢量的自相关和互相关

|名称   |符号  |位置 |
| ---   | ---  | --- |
|m0     |(+-1) | 0, 5, 10, 15, 20, 25, 30, 35  |
|m1     |(+-1) | 1, 6, 11, 16, 21, 26, 31, 36  |
|m2     |(+-1) | 2, 7, 12, 17, 22, 27, 32, 37  |
|m3     |(+-1) | 3, 8, 13, 18, 23, 28, 33, 38  |
|m4     |(+-1) | 4, 9, 14, 19, 24, 29, 34, 39  |

在[JBoss - org.mobicents.media.server.impl.dsp.audio.g729.LD8KConstants]中的对应部分有如下注释
```cpp
public static final short DIM_RR = 616; /* size of correlation matrix                            */
public static final short NB_POS = 8  ; /* Number of positions for each pulse                    */
public static final short STEP   = 5 ;  /* Step betweem position of the same pulse.              */
public static final short MSIZE  = 64 ; /* Size of vectors for cross-correlation between 2 pulses*/
```
从注释可已看出，NB_POS表示每个冲击可能的位置个数；STEP表示对于某个冲击，相邻位置的间隔；MSIZE表示两个冲击互相关可能取值的个数。
- 自相关`$R_{ii}$`个数 = 3 * NB_POS + (2*NB_POS) = 3 * 8 + (2*16) = 40   
    `$|\phi_{xx\_11}$` `$\phi_{xx\_22}$` `$\phi_{xx\_33}$` `$\phi_{xx\_44}$` `$\phi_{xx\_55}$` `$\phi_{xx\_66}$` `$\phi_{xx\_77}$` `$\phi_{xx\_88}|$`
- 互相关`$R_{ij}$`个数 = 3 * 64 + 3 * (8*16) = 192 + 3 * 128 = 576   
    `$|\phi_{xy\_11}$` `$\phi_{xy\_12}$` `$\phi_{xy\_13}$` `$\phi_{xy\_14}$` `$\phi_{xy\_15}$` `$\phi_{xy\_16}$` `$\phi_{xy\_17}$` `$\phi_{xy\_18}|$`
    
    `$|\phi_{yx\_12}$` `$\phi_{xy\_22}$` `$\phi_{xy\_23}$` `$\phi_{xy\_24}$` `$\phi_{xy\_25}$` `$\phi_{xy\_26}$` `$\phi_{xy\_27}$` `$\phi_{xy\_28}|$`
    
    `$|\phi_{yx\_13}$` `$\phi_{yx\_23}$` `$\phi_{xy\_33}$` `$\phi_{xy\_34}$` `$\phi_{xy\_35}$` `$\phi_{xy\_36}$` `$\phi_{xy\_37}$` `$\phi_{xy\_38}|$`
    
    `$|\phi_{yx\_14}$` `$\phi_{yx\_24}$` `$\phi_{yx\_34}$` `$\phi_{xy\_44}$` `$\phi_{xy\_45}$` `$\phi_{xy\_46}$` `$\phi_{xy\_47}$` `$\phi_{xy\_48}|$`
    
    `$|\phi_{yx\_15}$` `$\phi_{yx\_25}$` `$\phi_{yx\_35}$` `$\phi_{yx\_45}$` `$\phi_{xy\_55}$` `$\phi_{xy\_56}$` `$\phi_{xy\_57}$` `$\phi_{xy\_58}|$`
    
    `$|\phi_{yx\_16}$` `$\phi_{yx\_26}$` `$\phi_{yx\_36}$` `$\phi_{yx\_46}$` `$\phi_{yx\_56}$` `$\phi_{xy\_66}$` `$\phi_{xy\_67}$` `$\phi_{xy\_68}|$`
    
    `$|\phi_{yx\_17}$` `$\phi_{yx\_27}$` `$\phi_{yx\_37}$` `$\phi_{yx\_47}$` `$\phi_{yx\_57}$` `$\phi_{yx\_67}$` `$\phi_{xy\_77}$` `$\phi_{xy\_78}|$`
    
    `$|\phi_{yx\_18}$` `$\phi_{yx\_28}$` `$\phi_{yx\_38}$` `$\phi_{yx\_48}$` `$\phi_{yx\_58}$` `$\phi_{yx\_68}$` `$\phi_{yx\_78}$` `$\phi_{xy\_88}|$`
- 总个数 = 自相关个数 + 互相关个数 = **616**

代码当中处理时，位置3的16个可能取值被分为两组进行处理，而且计算时不计算这两组间的互相关，结果是一样的。


首先计算自相关，示例如下   
```math
\phi(i,i)=\sum_{n=i}^{39}h(n-i)h(n-i)

\phi(39,39)=h(0)h(0)

\phi(38,38)=h(0)h(0)+h(1)h(1)

\phi(37,37)=h(0)h(0)+h(1)h(1)+h(2)h(2)


...
```
代码如下：
```cpp
/*
 * i0 (+-1) : 0, 5, 10, 15, 20, 25, 30, 35                                *
 * i1 (+-1) : 1, 6, 11, 16, 21, 26, 31, 36                                *
 * i2 (+-1) : 2, 7, 12, 17, 22, 27, 32, 37                                *
 * i3 (+-1) : 3, 8, 13, 18, 23, 28, 33, 38                                *
 *            4, 9, 14, 19, 24, 29, 34, 39                                *
*/
  ptr_h1 = h;
  cor    = 0;
  for(i=0;  i<NB_POS; i++)
  {
    // fai_hihi, i=0,1,2,3,4....
    cor = L_mac(cor, *ptr_h1, *ptr_h1); ptr_h1++; // 0, 5, 10, 15, 20, 25, 30, 35
    *p4-- = extract_h(cor);

    cor = L_mac(cor, *ptr_h1, *ptr_h1); ptr_h1++;
    *p3-- = extract_h(cor);

    cor = L_mac(cor, *ptr_h1, *ptr_h1); ptr_h1++;
    *p2-- = extract_h(cor);

    cor = L_mac(cor, *ptr_h1, *ptr_h1); ptr_h1++;
    *p1-- = extract_h(cor);

    cor = L_mac(cor, *ptr_h1, *ptr_h1); ptr_h1++;
    *p0-- = extract_h(cor);
  }
```

计算互相关时，由`$\phi(i,j)$`的计算公式 
`$\phi(i,j)=\sum_{n=j}^{39}h(n-i)h(n-j)$` `$i=0,1,...,39$`,`$j=i,i+1,...,39$` 可知，
计算`$m_0$`和`$m_1$`的互相关性时，由于对索引的大小有要求，因此不是二者所有的组合都存在相关结果。
`$m_0$`和`$m_1$`存在两种关系：
- 关系1：`$m_0[i]<m_1[i]<m_1[i+1]<...<m_1[8]$`，`$R_{01}$`对应图中红色位置。
- 关系2：`$m_0[8]>...>m_0[i+1]>m_1[i]$`，`$R_{10}$`对应图中绿色位置。

![image](http://note.youdao.com/yws/public/resource/fedd85d892cd55c62ed9d238f41e4155/xmlnote/9003CB614B72435A8C5D981279424B4C/103970)

其他任意两行的互相关都存在类似关系。

以`$R_{01}$`、`$R_{12}$`、`$R_{23}$`关系1和`$R_{04}$`关系2为例：
```math
R_{01\_88}=\phi(37,38)=\sum_{n=38}^{39}h(n-37)h(n-38)=h(1)h(0)+h(2)h(1)

R_{12\_88}\phi(36,37)=\sum_{n=37}^{39}h(n-36)h(n-37)=h(1)h(0)+h(2)h(1)+h(3)h(2)

R_{23\_88}\phi(35,36)=\sum_{n=36}^{39}h(n-35)h(n-36)=h(1)h(0)+h(2)h(1)+h(3)h(2)+h(4)h(3)

R_{40\_78}\phi(34,35)=\sum_{n=35}^{39}h(n-34)h(n-35)=h(1)h(0)+h(2)h(1)+h(3)h(2)+h(4)h(3)+h(5)h(4)

...
```
对应代码如下：
```cpp
/*
 * i4 (+-1) :    4, 9,  14, 19, 24, 29, 34, 39                            *
 * i0 (+-1) : 0, 5, 10, 15, 20, 25, 30, 35                                *
 * i1 (+-1) : 1, 6, 11, 16, 21, 26, 31, 36                                *
 * i2 (+-1) : 2, 7, 12, 17, 22, 27, 32, 37                                *
 * i3 (+-1) : 3, 8, 13, 18, 23, 28, 33, 38                                *
*/
 /*-----------------------------------------------------------------*
  * Compute elements of: rri2i3[上], rri1i2[上], rri0i1[上] and rri0i4[下]  *
  *-----------------------------------------------------------------*/

  l_fin_sup = MSIZE-1;
  l_fin_inf = l_fin_sup-(Word16)1;
  ldec = NB_POS+1;

  ptr_hd = h;
  ptr_hf = ptr_hd + 1;

  for(k=0; k<NB_POS; k++) {
          // ROUND 1
          p3 = rri2i3 + l_fin_sup;
          p2 = rri1i2 + l_fin_sup;
          p1 = rri0i1 + l_fin_sup;
          p0 = rri0i4 + l_fin_inf;  // rri4i0
          cor = 0;
          ptr_h1 = ptr_hd;
          ptr_h2 =  ptr_hf;

          for(i=k+(Word16)1; i<NB_POS; i++ ) {
                  // round 1
                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // h0*h1
                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // R_37_38=h0*h1+h1*h2
                  *p3 = extract_h(cor);

                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // R_36_37=h0*h1+h1*h2+h2*h3
                  *p2 = extract_h(cor);

                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // R_35_36=h0*h1+h1*h2+h2*h3+h3*h4
                  *p1 = extract_h(cor);

                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // R_34_35=h0*h1+h1*h2+h2*h3+h3*h4+h4*h5
                  *p0 = extract_h(cor);

                  p3 -= ldec;   // -=9
                  p2 -= ldec;
                  p1 -= ldec;
                  p0 -= ldec;
                  
                  // round 2
                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // h5*h6
                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // R_32_33=h0*h1+h1*h2+h2*h3+h3*h4+h4*h5+h6*h7+
                  *p3 = extract_h(cor);

                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // R_31_32=h0*h1+h1*h2+h2*h3+h3*h4+h4*h5+h6*h7+h7*h8+
                  *p2 = extract_h(cor);

                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // R_30_31=h0*h1+h1*h2+h2*h3+h3*h4+h4*h5+h6*h7+h7*h8+h8*h9
                  *p1 = extract_h(cor);

                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // R_29_30=h0*h1+h1*h2+h2*h3+h3*h4+h4*h5+h6*h7+h7*h8+h8*h9+h9*h10
                  *p0 = extract_h(cor);

                  p3 -= ldec;   // -=9
                  p2 -= ldec;
                  p1 -= ldec;
                  p0 -= ldec;
                  
                  // round 3
                  ...
          }
          cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;           
          cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;           // R_2_3
          *p3 = extract_h(cor);

          cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;           // R_1_2
          *p2 = extract_h(cor);

          cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;           // R_0_1
          *p1 = extract_h(cor);

          l_fin_sup -= NB_POS;
          l_fin_inf--;
          ptr_hf += STEP;
          
          // ROUND 2
          p3 = rri2i3 + l_fin_sup;
          p2 = rri1i2 + l_fin_sup;
          p1 = rri0i1 + l_fin_sup;
          p0 = rri0i4 + l_fin_inf;  // rri4i0
          cor = 0;
          ptr_h1 = ptr_hd;
          ptr_h2 =  ptr_hf;

          for(i=k+(Word16)1; i<NB_POS; i++ ) {
                  // round 1
                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // h0*h6
                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // R_32_38=h0*h6+h1*h7
                  *p3 = extract_h(cor);

                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // R_31_37=h0*h6+h1*h7+h2*h8
                  *p2 = extract_h(cor);

                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // R_30_36=h0*h6+h1*h7+h2*h8+h3*h9
                  *p1 = extract_h(cor);

                  cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;   // R_29_35=h0*h6+h1*h7+h2*h8+h3*h9+h4*h10
                  *p0 = extract_h(cor);

                  p3 -= ldec;   // -=9
                  p2 -= ldec;
                  p1 -= ldec;
                  p0 -= ldec;
                  
                  // round 2
                  ...
          }
          cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;           
          cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;           // R_2_8
          *p3 = extract_h(cor);

          cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;           // R_1_7
          *p2 = extract_h(cor);

          cor = L_mac(cor, *ptr_h1, *ptr_h2); ptr_h1++; ptr_h2++;           // R_0_6
          *p1 = extract_h(cor);

          l_fin_sup -= NB_POS;
          l_fin_inf--;
          ptr_hf += STEP;
          
          // ROUND 3
          ...
  }
```

### 实际搜索

- 得到互相关值的符号
- 通过前3个冲击和滤波器的互相关，得到搜索门限
    - 搜索找出最大的三个值的和`$max_3$`
    - 搜索找出三个值和的平均值`$avg_3$`
    - 计算门限`$thr_3=avg_3+k*(max_3-avg_3)$`，`$k=0.4$`
- 计算互相关`$\phi$`，将冲击的符号考虑进去（只考虑互相关部分即可）


通过最大化`$\frac{correlation^2}{energy}$`的值，搜索4个冲击所在的最佳位置。实现时，搜索是通过4层嵌套实现的。每层嵌套当中，一个冲击会分别影响自相关和能量。   
只有当前3个冲击得到的自相关超过预设门限时，才会进入第4个循环。


```math
ps_0=d[i_0]

ps_1=ps_0+d[i_1]=d[i_0]+d[i_1]

ps_2=ps_1+d[i_2]=d[i_0]+d[i_1]+d[i_2]

ps_3=ps_2+d[i_3]=d[i_0]+d[i_1]+d[i_2]+d[i_3]

ps_3\_c=ps_3*ps_3
```

```math
\alpha_0=R_{0\_0}[i]

\alpha_1=\alpha_0 + R_{1\_1}[i] + 2 * R_{0\_1}[i]

\alpha_2=\alpha_1 + R_{2\_2}[i] + 2 * \left(R_{0\_2}[i]+R_{1\_2}[i]\right)

\alpha_3=\alpha_2 + R_{3\_3}[i] + 2 * \left(R_{0\_3}[i]+R_{1\_3}[i]+R_{2\_3}[i]\right)

```

把求商的最值，变换成乘积差的比较，得到最大值对应的冲击位置。
```math
ps_3*ps_3*\alpha - psc*\alpha_3 > 0
==>>
\frac{ps_3*ps_3}{alpha_3}>\frac{psc}{\alpha}

```


```cpp
// --
 time   = add(MAX_TIME, extra);     // 首帧搜索(75+30=105)次，其他帧等于(75+上次搜索剩余次数)。
 for (i0 = 0; i0 < L_SUBFR; i0 += STEP)        /* first pulse loop  */
 {
   ps0  = Dn[i0];
   alp0 = *ptr_ri0i0++;

   ptr_ri1i1 = rri1i1;    /* Init. pointers that depend on second loop */
   ptr_ri1i2 = rri1i2;
   ptr_ri1i3 = rri1i3;
   ptr_ri1i4 = rri1i4;
   for (i1 = 1; i1 < L_SUBFR; i1 += STEP)      /* second pulse loop */
   {
     ps1  = add(ps0, Dn[i1]);   // Dn[i0]+Dn[i1]

     /* alp1 = alp0 + *ptr_ri1i1++ + 2.0 * ( *ptr_ri0i1++); */
     alp1 = L_mult(alp0, 1);
     alp1 = L_mac(alp1, *ptr_ri1i1++, 1);
     alp1 = L_mac(alp1, *ptr_ri0i1++, 2);

     ptr_ri2i2 = rri2i2;     /* Init. pointers that depend on third loop */
     ptr_ri2i3 = rri2i3;
     ptr_ri2i4 = rri2i4;
     for (i2 = 2; i2 < L_SUBFR; i2 += STEP)    /* third pulse loop */
     {
       ps2  = add(ps1, Dn[i2]);

       /* alp2 = alp1 + *ptr_ri2i2++ + 2.0 * (*ptr_ri0i2++ + *ptr_ri1i2++); */
       alp2 = L_mac(alp1, *ptr_ri2i2++, 1);
       alp2 = L_mac(alp2, *ptr_ri0i2++, 2);
       alp2 = L_mac(alp2, *ptr_ri1i2++, 2);

       /* Test threshold */
       if ( sub(ps2, thres) > 0)
       {
         ptr_ri3i3 = rri3i3;    /* Init. pointers that depend on 4th loop */
         for (i3 = 3; i3 < L_SUBFR; i3 += STEP)      /* 4th pulse loop */
         {
           ps3 = add(ps2, Dn[i3]);

           /* alp3 = alp2 + *ptr_ri3i3++                                */
           /*       + 2.0*( *ptr_ri0i3++ + *ptr_ri1i3++ + *ptr_ri2i3++); */
           alp3 = L_mac(alp2, *ptr_ri3i3++, 1);
           alp3 = L_mac(alp3, *ptr_ri0i3++, 2);
           alp3 = L_mac(alp3, *ptr_ri1i3++, 2);
           alp3 = L_mac(alp3, *ptr_ri2i3++, 2);
           alp  = extract_l(L_shr(alp3, 5));

           ps3c = mult(ps3, ps3);   // L_temp = ps3*ps3*alpha - psc*alp
           L_temp = L_mult(ps3c, alpha);
           L_temp = L_msu(L_temp, psc, alp);
           if( L_temp > 0L )
           {
             psc = ps3c;
             alpha = alp;
             ip0 = i0;
             ip1 = i1;
             ip2 = i2;
             ip3 = i3;
           }
         }  /*  end of for i3 = */
         ptr_ri0i3 -= NB_POS;
         ptr_ri1i3 -= NB_POS;

         ptr_ri4i4 = rri4i4;    /* Init. pointers that depend on 4th loop */
         for (i3 = 4; i3 < L_SUBFR; i3 += STEP)      /* 4th pulse loop */
         {
           ps3 = add(ps2, Dn[i3]);

           /* alp3 = alp2 + *ptr_ri4i4++                                */
           /*       + 2.0*( *ptr_ri0i4++ + *ptr_ri1i4++ + *ptr_ri2i4++); */
           alp3 = L_mac(alp2, *ptr_ri4i4++, 1);
           alp3 = L_mac(alp3, *ptr_ri0i4++, 2);
           alp3 = L_mac(alp3, *ptr_ri1i4++, 2);
           alp3 = L_mac(alp3, *ptr_ri2i4++, 2);
           alp  = extract_l(L_shr(alp3, 5));

           ps3c = mult(ps3, ps3);
           L_temp = L_mult(ps3c, alpha);
           L_temp = L_msu(L_temp, psc, alp);
           if( L_temp > 0L )
           {
             psc = ps3c;
             alpha = alp;
             ip0 = i0;
             ip1 = i1;
             ip2 = i2;
             ip3 = i3;
           }
         }  /*  end of for i3 = */
         ptr_ri0i4 -= NB_POS;
         ptr_ri1i4 -= NB_POS;

         time = sub(time, 1);
         if(time <= 0 ) goto end_search;     // 限定搜索次数，超出后退出 /* Maximum time finish */

       }  /* end of if >thres */
       else
       {
         ptr_ri2i3 += NB_POS;
         ptr_ri2i4 += NB_POS;
       }

     } /* end of for i2 = */

     ptr_ri0i2 -= NB_POS;
     ptr_ri1i3 += NB_POS;
     ptr_ri1i4 += NB_POS;
   } /* end of for i1 = */

   ptr_ri0i2 += NB_POS;
   ptr_ri0i3 += NB_POS;
   ptr_ri0i4 += NB_POS;
 } /* end of for i0 = */
```



```math
y[i]=\sum_{i=0}^{n}{cod[i]h[n-i]}=\sum_{k=1}^{4}{cod[i_k]h[n-i_k]}
```

**i3 (+-1):**
```
3, 8, 13, 18, 23, 28, 33, 38                                
4, 9, 14, 19, 24, 29, 34, 39  
```
上下两行除以5得到的结果相同，因而必须额外增加一个信息表示是哪一行的位置。

```cpp
end_search:
 extra = time;  // 下次额外搜索次数，等于本次剩余次数（循环正常结束时））

 /* Set the sign of impulses */
 i0 = p_sign[ip0];  // 得到对应冲击的符号
 i1 = p_sign[ip1];
 i2 = p_sign[ip2];
 i3 = p_sign[ip3];

 /* Find the codeword corresponding to the selected positions */
 for(i=0; i<L_SUBFR; i++) {cod[i] = 0; }    // 得到搜索结果对应的码字
 cod[ip0] = shr(i0, 2);         /* From Q15 to Q13 */
 cod[ip1] = shr(i1, 2);
 cod[ip2] = shr(i2, 2);
 cod[ip3] = shr(i3, 2);

 /* find the filtered codeword */
 // 只需计算非0位置的冲击即可。
 for (i = 0; i < L_SUBFR; i++) {y[i] = 0;  }

 if(i0 > 0)
   for(i=ip0, j=0; i<L_SUBFR; i++, j++) {
       y[i] = add(y[i], h[j]); }
 else
   for(i=ip0, j=0; i<L_SUBFR; i++, j++) {
       y[i] = sub(y[i], h[j]); }

 if(i1 > 0)
   for(i=ip1, j=0; i<L_SUBFR; i++, j++) {
       y[i] = add(y[i], h[j]); }
 else
   for(i=ip1, j=0; i<L_SUBFR; i++, j++) {
       y[i] = sub(y[i], h[j]); }

 if(i2 > 0)
   for(i=ip2, j=0; i<L_SUBFR; i++, j++) {
       y[i] = add(y[i], h[j]); }
 else
   for(i=ip2, j=0; i<L_SUBFR; i++, j++) {
       y[i] = sub(y[i], h[j]); }

 if(i3 > 0)
   for(i=ip3, j=0; i<L_SUBFR; i++, j++) {
       y[i] = add(y[i], h[j]); }
 else
   for(i=ip3, j=0; i<L_SUBFR; i++, j++) {
       y[i] = sub(y[i], h[j]); }

 /* find codebook index;  17-bit address */
 // 保存符号位
 i = 0;
 if(i0 > 0) i = add(i, 1);
 if(i1 > 0) i = add(i, 2);
 if(i2 > 0) i = add(i, 4);
 if(i3 > 0) i = add(i, 8);
 *sign = i;
 // 保存码字对应位置索引
 ip0 = mult(ip0, 6554);         /* ip0/5 */
 ip1 = mult(ip1, 6554);         /* ip1/5 */
 ip2 = mult(ip2, 6554);         /* ip2/5 */
 i   = mult(ip3, 6554);         /* ip3/5 */
 j   = add(i, shl(i, 2));       /* j = i*5 */
 j   = sub(ip3, add(j, 3));     /* j= ip3%5 -3 */
 ip3 = add(shl(i, 1), j);       // 行信息，保存在ip3的最后一个bit位
 // 保存到一个Word16当中
 i = add(ip0, shl(ip1, 3));
 i = add(i  , shl(ip2, 6));
 i = add(i  , shl(ip3, 9));

 return i; // 保存着冲击位置的信息
```


计算固定码字的增益，计算时考虑自适应码书的结果。

如果基音延迟`$T_0<40$`，则将`$[T_0,40)$`之间的的值重新用以下公式计算。
```math
code[i] = code[i]+pitch_sharp*cod[i-T_0],i=[T_0,40)
```

```cpp
 /*-----------------------------------------------------------------*
  * Compute innovation vector gain.                                 *
  * Include fixed-gain pitch contribution into code[].              *
  *-----------------------------------------------------------------*/

  if(sub(T0 ,L_SUBFR) <0)
     for (i = T0; i < L_SUBFR; i++) {  /* code[i] += pitch_sharp*code[i-T0] */
       code[i] = add(code[i], mult(code[i-T0], sharp));
     }
```
