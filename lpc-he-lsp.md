# 整体介绍

线性预测的基本原理和语音信号的数字模型密切相关。

语音信号存在两种类型的相关性，即在样点之间的短时相关性（声道特性），和相邻基音周期之间的长时相关性（激励信号特性）。相邻样本之间的预测利用了比较相邻的样本值，所以称为“短时预测”，它实际上是频谱包络的预测；而为了区别于短时预测，将基于**基音**周期的预测称为“长时预测”，它实际上是基于频谱的细微结构的预测。

由于语音样点之间存在短时相关性，所以可以用过去的样点值来预测现在或未来的样点值，即一个语音的抽样能够用过去若干个语音抽样或它们的线性组合来逼近。通过使实际语音抽样和线性预测抽样之间的误差在某个准则下达到最小值来决定唯一的一组预测系数。而这组预测系数就反应了语音信号的特性，可以作为语音信号特征参数用于语音识别、语音合成等。

LPC即为预测语音信号的短时相关性，即发声模型的声道特性。清音信号由随机噪声通过声道激励产，更适用与LPC分析；浊音信号具有准周期性，所以相邻周期的样本之间具有很大的相关性，LPC分析结果不仅体现了声道特性，还体现了激励的特性。

```
graph TD
    A[海明窗处理] --> B[计算自相关]
    B --> C[Lag窗]
    C --> |LevinsionDurbin|D[LPC]
    D --> |LPC->LSP|E[LSP]
    E --> |Chepseve多项式|H[LPC_R和LPC_RQ]
    E --> F1[LSP量化]
    F1 --> |LSP->LPC| H
```

整体代码：

```c
  /* LP analysis */
  Autocorr(p_window, M, r_h, r_l);              /* Autocorrelations */
  Lag_window(M, r_h, r_l);                      /* Lag windowing    */
  Levinson(r_h, r_l, &A_t[MP1],rc);             /* Levinson Durbin  */
  Az_lsp(&A_t[MP1], lsp_new, lsp_old);          /* From A(z) to lsp */

  /* LSP quantization */
  Qua_lsp(lsp_new, lsp_new_q, ana);
  ana += 2;                         /* Advance analysis parameters pointer */

 /*--------------------------------------------------------------------*
   * Find interpolated LPC parameters in all subframes (both quantized  *
   * and unquantized).                                                  *
   * The interpolated parameters are in array A_t[] of size (M+1)*4     *
   * and the quantized interpolated parameters are in array Aq_t[]      *
   *--------------------------------------------------------------------*/
  Int_lpc(lsp_old, lsp_new, lsf_int, lsf_new,  A_t);
  Int_qlpc(lsp_old_q, lsp_new_q, Aq_t);

  /* update the LSPs for the next frame */
  for(i=0; i<M; i++)
  {
    lsp_old[i]   = lsp_new[i];
    lsp_old_q[i] = lsp_new_q[i];
  }
```

# 海明窗处理

语音信号的数字化一般包括放大及增益控制、防混叠滤波、采样、AD转换及编码（一般是PCM编码）。  
1、预滤波的目的有两个：抑制输入信号各频率中所有超过采样率一半的分量，防止混叠；防止50Hz电源工频干扰。  
2、AD采样一般存在量化误差，或称为量化噪声。

预处理一般包括：预加重（提升高频）、加窗、分帧等。此外，在分析处理之前必须要把要分析的语音信号部分从输入信号中照出来，这项工作叫做语音信号的端点检测。  
1、预加重：由于语音信号的平均功率谱收到声门激励和口鼻辐射的影响，高频大约在800Hz以上按6dB/倍频程跌落，及6dB/2倍频或20dB/10倍频。所以求语音信号的频谱时，频率越高相应的成分越小，高频部分的频谱比低频部分的频谱难求，因此要在预处理中进行预加重处理。

预加重处理的目的是提升高频部分，是信号的频谱变得平坦，保持在低频到高频的整个频带中，能用同样的信噪比求频谱，以便于频谱分析或声道参数分析。预加重可以在防混叠之前进行，这样不仅可以进行预加重，还可以压缩信号的动态范围，有效提高信噪比。  
    在恢复原信号时，要进行相应的去加重处理。

2、加窗分帧。分帧一般采用交叠分段的方法，这是为了使帧与帧之间平滑过渡，保持其连续性。分帧时，窗口的选择（形状和长度）对分析参数的特性影响很大，为此应选择合适的窗口，使短时参数更好地反应语音信号的特性变化。

窗口长度的选择，更重要的是考虑语音信号的基音周期。通常认为在一个语音帧内应包含1-7个基音周期，然而不同人的基音周期变化很大，从女性和儿童的2ms到老年男人的14ms（500-70Hz），所以N的选择比较困难。通常在10Khz取样频率下，N折衷选择为100-200点为宜（10-20ms的持续时间）。  
    采样周期TS=1/fs，窗口长度N，频率分辨率△f之间存在如下关系：  △f=1/\(N \* TS\)

**窗函数的作用**：截取信号，并尽量不影响信号的频域特性，所以窗函数的频域特性，应该是主瓣窄，并且旁瓣低。由于同时满足两者有冲突，比如频域上是个冲激\(这是最理想的，不会改变信号的任何频域特性\)，但在时域则为无穷个点，这种窗显然没有任何实际意义。于是数学家们设计了一些窗函数，分别适用于某种特定应用场合，海明窗是其中的一个。

语音信号的数字化和预处理是一个很重要的环节，在对一个语音信号处理系统进行性能评价时，作为语音参数分析条件，采样频率和精度、采用了什么样的预加重、窗函数、帧长、和帧移是多少都必须交代清楚以供参考。

```c
  for(i=0; i<L_WINDOW; i++)
  {
    y[i] = mult_r(x[i], hamwindow[i]);
  }
```

# 计算自相关

**防止溢出**:  
首先计算自相关r\[0\]，如果计算过程中发生溢出，则先将输入右移2位（除以4）再重新计算r\[0\]，以此类推，直到不在溢出为止。

```c
  do {
    Overflow = 0;
    sum = 1;                   /* Avoid case of all zeros */
    for(i=0; i<L_WINDOW; i++)
      sum = L_mac(sum, y[i], y[i]); // y[i]为加窗后的数据，sum + y[i]*y[i]
                                    // sum至少为1，如果有溢出L_mac中置位全局变量Overflow

    /* If overflow divide y[] by 4 */
    if(Overflow != 0)
    {
      for(i=0; i<L_WINDOW; i++)
      {
        y[i] = shr(y[i], 2); // 右移2位
      }
    }
  }while (Overflow != 0);   // 直到不再溢出为止
```

**计算归一化参数**:  
计算可以将一个32位变量L\_var1通过左移进行归一化需要的位数。  
即将L\_var1看做小数，通过移位的方式使其尽可能接近1（0x7fffffff），达到除法归一化的目的。

```c
Word16 norm_l(Word32 L_var1)
{
  Word16 var_out;

  if (L_var1 == 0){ // 0 << 0 = 0
    var_out = 0;
  }  else  {
    if (L_var1 == (Word32)0xffffffffL) { // -1，直接取反则会出现全0的情况
      var_out = 31;
    }    else    {
      if (L_var1 < 0)      {
        L_var1 = ~L_var1; // 取反表示对应的整数
      }

      for(var_out = 0;L_var1 < (Word32)0x40000000L;var_out++) { // var_out 次高位为1 停止（最高位为符号位）
        L_var1 <<= 1;
      }
    }
  }

   return(var_out);
}
```

**计算自相关，并归一化**:  
由于前面已经保证了自相关最大值r\[0\]不出现溢出，所以互相关一定不会溢出。同时，利用前面计算得到的归一化参数，通过移位的方式实现了自相关函数的归一化。

```c
  for (i = 1; i <= m; i++)
  {
    sum = 0;
    for(j=0; j<L_WINDOW-i; j++)
      sum = L_mac(sum, y[j], y[j+i]);

    sum = L_shl(sum, norm);
    L_Extract(sum, &r_h[i], &r_l[i]);
  }
```

# LAG window

1. 带宽扩展，避免可能会干扰基音和共振峰参数的不正常剧变谐波。
2. 增强Levinson-Durbin算法的稳定性。
   ```c
   for(i=1; i<=m; i++)
   {
     x  = Mpy_32(r_h[i], r_l[i], lag_h[i-1], lag_l[i-1]); // 
     L_Extract(x, &r_h[i], &r_l[i]);
   }
   ```

# Levinson-Durbin算法计算LPC 【待续】

LPC是对自相关的LPC？

```
//--
 |        R[i]    autocorrelations.                                          |
 |       A[i]    filter coefficients.                                        |
 |       K       reflection coefficients.                                    |
 |       Alpha   prediction gain.                                            |
 |                                                                           |
 |       Initialization:                                                     |
 |               A[0] = 1                                                    |
 |               K    = -R[1]/R[0]                                           |
 |               A[1] = K                                                    |
 |               Alpha = R[0] * (1-K**2]                                     |
 |                                                                           |
 |       Do for  i = 2 to M                                                  |
 |            S =  SUM ( R[j]*A[i-j] ,j=1,i-1 ) +  R[i]                      |
 |            K = -S / Alpha                                                 |
 |            An[j] = A[j] + K*A[i-j]   for j=1 to i-1                       |
 |                                      where   An[i] = new A[i]             |
 |            An[i]=K                                                        |
 |            Alpha=Alpha * (1-K**2)                                         |
 |                                                                           |
 |       END                                                                 |
```

# LSP分析

由于LPC的系数`$a_i$`是在时域的，其较小的量化误差就会影响到LPC滤波器响应的每个频段，因此不适宜直接进行量化传输。  
将其转化为隐式分解的形式，然后分解后的因式求根，每个根的误差只会影响某频段的能量，这样就适合进行矢量量化了。  
LSF系数实际是A\(z\)系数的变形，并且LSF所有的根都在单位圆上，求其系数可以有很多简便的方法。

## A\(z\) To LSP

g729 lpc与lsp的转换用的是**切比雪夫多项式**来逼近求单位圆上的根。

### 切比雪夫多项式系数计算f\(i\)

A\(z\)是10阶LPC滤波器，a\(i\)是A\(z\)的系数，fn\(i\)是辅助多项式n对应的切比雪夫多项式多项式的的参数。  
f\(i\)计算公式公式如下：

```cpp
/*-------------------------------------------------------------*
 * f1[0] = 1.0;                                                *
 * f2[0] = 1.0;                                                *
 *                                                             *
 * for (i = 0; i< NC; i++)                                     *
 * {                                                           *
 *   f1[i+1] = a[i+1] + a[M-i] - f1[i] ;                       *
 *   f2[i+1] = a[i+1] - a[M-i] + f2[i] ;                       *
 * }                                                           *
 *-------------------------------------------------------------*/
```

Chebps\_10或Chebps\_11是计算切比雪夫多项式系数的函数，首先选择Q11的定点，如果出现溢出，则采用Q10定点。  
代码实现如下：

```cpp
 for (i = 0; i< NC; i++)
 {
   Overflow = 0;
   t0 = L_mult(a[i+1], 16384);          /* x = (a[i+1] + a[M-i]) >> 1        */
   t0 = L_mac(t0, a[M-i], 16384);       /*    -> From Q12 to Q11             */
   x  = extract_h(t0);
   if ( Overflow ) 
   {
     ovf_coef = 1;      
   }

   Overflow = 0;
   f1[i+1] = sub(x, f1[i]);             /* f1[i+1] = a[i+1] + a[M-i] - f1[i] */
   if ( Overflow ) 
   {
     ovf_coef = 1;      
   }

   Overflow = 0;
   t0 = L_mult(a[i+1], 16384);          /* x = (a[i+1] - a[M-i]) >> 1        */
   t0 = L_msu(t0, a[M-i], 16384);       /*    -> From Q12 to Q11             */
   x  = extract_h(t0);
   if ( Overflow ) 
   {
     ovf_coef = 1;      
   }

   Overflow = 0;
   f2[i+1] = add(x, f2[i]);             /* f2[i+1] = a[i+1] - a[M-i] + f2[i] */
   if ( Overflow ) 
   {
     ovf_coef = 1;      
   }
 }

 if ( ovf_coef ) // Q11溢出，采用Q10
 {
   /*printf("===== OVF ovf_coef =====\n");*/

   pChebps = Chebps_10;     // 默认是Chebps_11

   f1[0] = 1024;            /* f1[0] = 1.0 is in Q10 */
   f2[0] = 1024;             /* f2[0] = 1.0 is in Q10 */
   for (i = 0; i< NC; i++)
   {
     t0 = L_mult(a[i+1], 8192);         /* x = (a[i+1] + a[M-i]) >> 1        */
     t0 = L_mac(t0, a[M-i], 8192);      /*    -> From Q11 to Q10             */
     x  = extract_h(t0);
     f1[i+1] = sub(x, f1[i]);           /* f1[i+1] = a[i+1] + a[M-i] - f1[i] */

     t0 = L_mult(a[i+1], 8192);         /* x = (a[i+1] - a[M-i]) >> 1        */
     t0 = L_msu(t0, a[M-i], 8192);      /*    -> From Q11 to Q10             */
     x  = extract_h(t0);
     f2[i+1] = add(x, f2[i]);           /* f2[i+1] = a[i+1] - a[M-i] + f2[i] */
   }
 }
```

## 计算LSP参数

Az\_lsp函数最后,会判断搜索出来的根有几个,如果个数不足，仍然沿用旧的lsp系数

利用直线两点式方程，线性内插，估计过零点。

```math
(x-x_l)/(y-y_l)=(x_h-x_l)/(y_h-y_l)

y_l =0

x=x_l - y_l*(x_h-x_l)/(y_h-y_l)
```

代码如下：

```cpp
 nf=0;           /* number of found frequencies */
 ip=0;          /* indicator for f1 or f2      */
 coef = f1;
 xlow = grid[0];    // [0-pi)三角函数表，Q15, grid[i] = cos((6.283185307*i)/(2.0*grid_points));
 ylow = (*pChebps)(xlow, coef, NC); // 区间结束点

 j = 0;
 while ( (nf < M) && (j < GRID_POINTS) )
 {
   j =add(j,1);
   xhigh = xlow;
   yhigh = ylow;    // 区间起始点
   xlow  = grid[j];
   ylow  = (*pChebps)(xlow,coef,NC);    // 区间结束点

   L_temp = L_mult(ylow ,yhigh);        
   if ( L_temp <= (Word32)0)            // 区间中含有零点
   {
        /* divide 4 times the interval */
        for (i = 0; i < 4; i++)
        {
           xmid = add( shr(xlow, 1) , shr(xhigh, 1)); /* xmid = (xlow + xhigh)/2 */
           ymid = (*pChebps)(xmid,coef,NC);
           L_temp = L_mult(ylow,ymid);
           if ( L_temp <= (Word32)0)        // 零点在低半区间
           {
             yhigh = ymid;
             xhigh = xmid;
           }
           else                             // 零点在高半区间
           {
             ylow = ymid;
             xlow = xmid;
           }
        }

        /*-------------------------------------------------------------*
         * Linear interpolation                                        *
         *    xint = xlow - ylow*(xhigh-xlow)/(yhigh-ylow);            *
         *-------------------------------------------------------------*/
        x = sub(xhigh, xlow);
        y = sub(yhigh, ylow);

        if(y == 0)
        {
            xint = xlow;
        }
        else    // 线性内插近似过零点
        {
            sign= y;
            y   = abs_s(y);
            exp = norm_s(y);
            y   = shl(y, exp);
            y   = div_s( (Word16)16383, y);
            t0  = L_mult(x, y);
            t0  = L_shr(t0, sub(20, exp) );
            y   = extract_l(t0);            /* y= (xhigh-xlow)/(yhigh-ylow) in Q11 */

            if(sign < 0) 
            {
                y = negate(y);
            }

            t0   = L_mult(ylow, y);                  /* result in Q26 */
            t0   = L_shr(t0, 11);                    /* result in Q15 */
            xint = sub(xlow, extract_l(t0));         /* xint = xlow - ylow*y */
        }

        lsp[nf] = xint;
        xlow    = xint;
        nf =add(nf,1);  // 已经找到的频率个数

        if(ip == 0)     // 较差计算另一个辅助多项式的根
        {
            ip = 1;
            coef = f2;
        }
        else
        {
            ip = 0;
            coef = f1;
        }

        ylow = (*pChebps)(xlow,coef,NC);
    }
 }
```

## LSP参数量化

LSP和LSF之间的转换通过查找表+线性内插完成。

## 将LSP转换为LSF

```cpp
  ind = 63;           /* begin at end of table2 -1 */
  for(i= m-(Word16)1; i >= 0; i--)
  {
    /* find value in table2 that is just greater than lsp[i] */
    while( sub(table2[ind], lsp[i]) < 0 )   // table2为Q15的cos(x)表，x=(ind<<8)*2pi
    {
      ind = sub(ind,1);
      if ( ind <= 0 )
        break;
    }
    offset = sub(lsp[i], table2[ind]);  // 实际lsp值与查表所得之间的差值

    // 不能精确预测的差值，使用线性插值估计得到
    /* acos(lsp[i])= ind*512 + (slope_acos[ind]*offset >> 11) */
    L_tmp  = L_mult( slope_acos[ind], offset );   /* L_tmp in Q28 */ 
    freq = add(shl(ind, 9), extract_l(L_shr(L_tmp, 12)));
    lsf[i] = mult(freq, 25736);           /* 25736: 2.0*PI in Q12 */
  }
```

### 计算预测器权重矢量

```cpp
    buf[0] = sub( flsp[1], (PI04+8192) );  // PIO4即pi*0.04, buf Q13         /* 8192:1.0(Q13) */
    for ( i = 1 ; i < M-1 ; i++ ) 
    {
        tmp = sub( flsp[i+1], flsp[i-1] );
        buf[i] = sub( tmp, 8192 );
    }

    buf[M-1] = sub( (PI92-8192), flsp[M-2] );   // 

    // 计算 wi
    for ( i = 0 ; i < M ; i++ ) 
    {
        if ( buf[i] > 0 )
        {
            wegt[i] = 2048;                    /* 2048:1.0(Q11) */
        }
        else 
        {
            L_acc = L_mult( buf[i], buf[i] );           /* L_acc in Q27 */
            tmp = extract_h( L_shl( L_acc, 2 ) );       /* tmp in Q13 */

            L_acc = L_mult( tmp, CONST10 );             // CONST10 10.0, /* L_acc in Q25 */
            tmp = extract_h( L_shl( L_acc, 2 ) );       /* tmp in Q11 */

            wegt[i] = add( tmp, 2048 );                 /* wegt in Q11 */
        }
    }

    // 特殊处理
    L_acc = L_mult( wegt[4], CONST12 );             // CONST12 1.2, /* L_acc in Q26 */
    wegt[4] = extract_h( L_shl( L_acc, 1 ) );       /* wegt in Q11 */

    L_acc = L_mult( wegt[5], CONST12 );             /* L_acc in Q26 */
    wegt[5] = extract_h( L_shl( L_acc, 1 ) );       /* wegt in Q11 */

    /* wegt: Q11 -> normalized */
    tmp = 0;
    for ( i = 0; i < M; i++ )       // 寻找最大值，寻求归一化
    {
        if ( sub(wegt[i], tmp) > 0 ) 
        {
            tmp = wegt[i];
        }
    }
    sft = norm_s(tmp);
    for ( i = 0; i < M; i++ ) 
    {
        wegt[i] = shl(wegt[i], sft);                  /* wegt in Q(11+sft) */
    }
```

### 预测

```cpp
for(mode = 0; mode<MODE; mode++)  // 预测的模式，2种
{
    Lsp_prev_extract(lsp, rbuf, fg[mode], freq_prev, fg_sum_inv[mode]);
    Lsp_pre_select(rbuf, lspcb1, &cand_cur );
    cand[mode] = cand_cur;
    Lsp_select_1(rbuf, lspcb1[cand_cur], wegt, lspcb2, &index);

    tindex1[mode] = index;
    for( j = 0 ; j < NC ; j++ )
    {
        buf[j] = add( lspcb1[cand_cur][j], lspcb2[index][j] );
    }

    Lsp_expand_1(buf, GAP1);
    Lsp_select_2(rbuf, lspcb1[cand_cur], wegt, lspcb2, &index);
    tindex2[mode] = index;

    for( j = NC ; j < M ; j++ )
    {
        buf[j] = add( lspcb1[cand_cur][j], lspcb2[index][j] );
    }
    Lsp_expand_2(buf, GAP1);
    Lsp_expand_1_2(buf, GAP2);
    Lsp_get_tdist(wegt, buf, &L_tdist[mode], rbuf, fg_sum[mode]);
}

Lsp_last_select(L_tdist, &mode_index);

code_ana[0] = shl( mode_index,NC0_B ) | cand[mode_index];
code_ana[1] = shl( tindex1[mode_index],NC1_B ) | tindex2[mode_index];

Lsp_get_quant(lspcb1, lspcb2, cand[mode_index],
      tindex1[mode_index], tindex2[mode_index],
      fg[mode_index], freq_prev, lspq, fg_sum[mode_index]) ;
```

### 将LSF转换为LSP

`ind    = b8-b15 of freq`  
`offset = b0-b7  of freq`

```cpp
  for(i=0; i<m; i++)
  {
    /*    freq = abs_s(freq);*/
    freq   = mult(lsf[i], 20861);        /* 20861: 1.0/(2.0*PI) in Q17 */
    ind    = shr(freq, 8);               /* ind    = b8-b15 of freq */
    offset = freq & (Word16)0x00ff;      /* offset = b0-b7  of freq */

    if ( sub(ind, 63)>0 ){
      ind = 63;                 /* 0 <= ind <= 63 */
    }

    /* lsp[i] = table2[ind]+ (slope_cos[ind]*offset >> 12) */
    // 差值对应的弧度差值，用线性估计得到
    L_tmp   = L_mult(slope_cos[ind], offset);   /* L_tmp in Q28 */
    lsp[i] = add(table2[ind], extract_l(L_shr(L_tmp, 13)));
  }
```

### 非量化LPC内插

上采样2\(M+1\)倍  
同时计算得到内插后的A\(z\)的参数和LSF的参数，以及内插钱非量化的LSF参数。

```math
lsp[i] = lsp\_new[i] * 0.5 + lsp\_old[i] * 0.5
```

```cpp
  /*  lsp[i] = lsp_new[i] * 0.5 + lsp_old[i] * 0.5 */

  for (i = 0; i < M; i++) 
  {
    lsp[i] = add(shr(lsp_new[i], 1), shr(lsp_old[i], 1));
  }

  Lsp_Az(lsp, Az);      // LSP 计算 Az的系数

  // 根据内插得到的LSP计算LSF参数（查表+线性估计）
  Lsp_lsf(lsp, lsf_int, M);      /* transformation from LSP to LSF (freq.domain) */
  Lsp_lsf(lsp_new, lsf_new, M);  /* transformation from LSP to LSF (freq.domain) */
```

### Lsp\_lsf/Lsf\_lsp与Lsp\_lsf2/Lsf\_lsp2的差别

Lsp\_lsf2\(\)和Lsf\_lsp2\(\)用于处理量化的LSP和LSF之间的转换，存储的LSF值都转化为了实际的弧度值（pi）为单位。Lsp\_lsf2\(\)在将LSP转换为LSF时，将最后的转换结果乘以`${2\pi}$`。

```cpp
...
lsf[i] = mult(freq, 25736);           /* 25736: 2.0*PI in Q12 */
```

而Lsf\_lsp2\(\)在将LSF转换为LSP时，将最后的转换结果乘以`$\frac{1}{2\pi}$`。

```cpp
freq = mult(lsf[i], 20861);          /* 20861: 1.0/(2.0*PI) in Q17 */
...
```

而它们对应的Lsp\_lsf\(\)和Lsf\_lsp\(\)则没有这两步操作。

### 量化LPC内插

上采样2\(M+1\)倍，后续按照2个子帧处理。  
量化LPC的LSF参数在LPC量化时已经进行过计算。

```cpp
  /*  lsp[i] = lsp_new[i] * 0.5 + lsp_old[i] * 0.5 */
  for (i = 0; i < M; i++) 
  {
    lsp[i] = add(shr(lsp_new[i], 1), shr(lsp_old[i], 1));
  }

  Lsp_Az(lsp, Az);              /* Subframe 1 */
  Lsp_Az(lsp_new, &Az[MP1]);    /* Subframe 2 */    //MP1=M+1, M是LPC滤波器阶数
```

## 更新LSP参数

更新LSP参数处理下一帧是使用。

```cpp
for(i=0; i<M; i++)
{
    lsp_old[i]   = lsp_new[i];
    lsp_old_q[i] = lsp_new_q[i];
}
```



