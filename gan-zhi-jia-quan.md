# 感知加权

首先构造出感知加权滤波器`$W(z)$`和总和滤波器`$A(z)$`，并将它们串联。

```cpp
 /*----------------------------------------------------------------------*
  * - Find the weighting factors                                         *
  *----------------------------------------------------------------------*/

  perc_var(gamma1, gamma2, lsf_int, lsf_new, rc);

/*----------------------------------------------------------------------*
  * - Find the weighted input speech w_sp[] for the whole speech frame   *
  * - Find the open-loop pitch delay                                     *
  *----------------------------------------------------------------------*/

  Weight_Az(&A_t[0], gamma1[0], M, Ap1);
  Weight_Az(&A_t[0], gamma2[0], M, Ap2);
  Residu(Ap1, &speech[0], &wsp[0], L_SUBFR);
  Syn_filt(Ap2, &wsp[0], &wsp[0], L_SUBFR, mem_w, 1);

  Weight_Az(&A_t[MP1], gamma1[1], M, Ap1);
  Weight_Az(&A_t[MP1], gamma2[1], M, Ap2);
  Residu(Ap1, &speech[L_SUBFR], &wsp[L_SUBFR], L_SUBFR);
  Syn_filt(Ap2, &wsp[L_SUBFR], &wsp[L_SUBFR], L_SUBFR, mem_w, 1);

  /* Find open loop pitch lag */
  T_op = Pitch_ol(wsp, PIT_MIN, PIT_MAX, L_FRAME);
```

## 计算Gamma值

基于非量化的LP滤波器系数`$a_i$`  
G729的共振峰感知加权是自适应的，perc\_var这个函数来对共振峰感加权进行估值。ITU文档3.3节公式30是判断语音普是否平坦的一个条件。人类语音谱有个特点，越高频的共振峰能量会越弱，共振峰感知加权要注意这个现象，如果频谱倾斜了（高频共振峰能量弱），要加强加权系数。

### 计算对数面积比

分段线性化，求解近似值。  
对数面积比趋向于负无穷的时候，此时的对应的`$\frac{R(1)}{R(0)}$`是接近于1的，趋近于正无穷时，对应的`$\frac{R(1)}{R(0)}$`接近于-1。即对数面积比系数越小，可以估计出高频分量越低，对数面积比系数越大，高频分量就越高。

当条件1满足时,我们其实可以计算出 k1差不多为`$$`-55/56，即 R\(1\)/R\(0\) = 55/56 即高频分量低R\(2\)/R\(0\)其实对应着次高频分量,再来看k2,R\(2\)/R\(0\)大概可以由k2的取值推断出R\(2\)/R\(0\)也是一个非常接近1的值,即,次高频分量还是低的.\(此时k2大约为3/5,由于分母极小,所以分子也不可能太大,推断于R\(2\)/R\(0\)是接近于1的值\)

结合前一帧是平的,由于高频分量都低了,就可以推断出当前帧应该是倾斜了

同理可以推断于条件2是一个相反的过程,由倾斜,而高频分量高了,认为当前帧的频谱是平坦的  
注:人类语音能量集中在前两个共振峰,后面的共振峰能量会依次降低,共振峰对听觉心理的影响最大,  
所以感知加权要加强共振峰的强度与带宽引入的感知加权.  
我们可以因式分解成:

```
          1
 -----------------------------------------------------------------------
 (z^-1 + r*cos(b1) + i*r*sin(b1)) ... (z^-1 + r*cos(b10) + i*r*sin(b10))
我们将 z=(z/a)  0<a<1代入,看其中的一个因子
             1
  ------------------------------------
  ((z/a)^-1 + r*cos(b1) + i*r*sin(b1))

可化为
            1/a
 ------------------------------------------
   (z^-1 + (r/a)*cos(b1) + i*(r/a)*sin(b1))
将 z = e^jw代入,则整个分式的绝对值(对应频域的振幅)
             1/a
 ----------------------------
  1+(r/a)^2+ 2*(r/a)cos(w+b1)
```

我们画一下这个函数的图,观察a变化,引起的幅度变化,自然能得出a越小,共振峰带宽扩展,共振峰加强的结论,b1则表示共振峰的位置。  
对比723,729引入了感知加权自适应机制,即,两个共振峰的位置如果太接近了,带宽扩展有可能导致两个共振峰出现重合了  
也就是最近的两个共振峰越接近,则带宽扩展的加权越低\(即a越大\)

以上就是分析出来的共振峰加权系数取值的一些推导依据,代码就相应简单了,基本上照本宣科  
最终求得两个共振峰感知加权系数。

得到加权系数后,对信号进行滤波,比较简单,不详述了

```math
Lar(i) = \log{\frac{1+rc}{1-rc}},(i=1.2)
```

其中rc为反射系数，在Levinsion-Durbin递推时得到。

```cpp
// perc_var()
  LarNew = &Lar[2];
  /* ---------------------------------------- */
  /* Reflection coefficients ---> Lar         */
  /* Lar(i) = log10( (1+rc) / (1-rc) )        */
  /* Approximated by                          */
  /* x <= SEG1            y = x               */
  /* SEG1 < x <= SEG2     y = A1 x - B1_L     */
  /* SEG2 < x <= SEG3     y = A2 x - B2_L     */
  /* x > SEG3             y = A3 x - B3_L     */
  /* ---------------------------------------- */
  for (i=0; i<2; i++) 
  {
    cur_rc = abs_s(r_c[i]);     // 计算绝对值
    cur_rc = shr(cur_rc, 4);

    if (sub(cur_rc ,SEG1) <= 0)     // [0, 0.6341], SEG1, 0.6341 in Q11
    {
        LarNew[i] = cur_rc; 
    }
    else 
    {
      if (sub(cur_rc,SEG2)<= 0)     // (0.6341, 0.8864], SEG2, 0.8864 in Q11 
      {
        cur_rc = shr(cur_rc, 1);
        L_temp = L_mult(cur_rc, A1);    // A1, 2.23   in Q11
        L_temp = L_sub(L_temp, L_B1);   // L_B1, 0.78   in Q22 
        L_temp = L_shr(L_temp, 11); 
        LarNew[i] = extract_l(L_temp);
      }
      else 
      {
        if (sub(cur_rc ,SEG3)<= 0)  // (0.8864, 0.9490], SEG3, 0.9490 in Q11 
        {
          cur_rc = shr(cur_rc, 1);
          L_temp = L_mult(cur_rc, A2);  // A2, 5.75   in Q11 
          L_temp = L_sub(L_temp, L_B2); // L_B2, 3.90   in Q22
          L_temp = L_shr(L_temp, 11);
          LarNew[i] = extract_l(L_temp);
        }
        else                        // (0.9490, 1]
        {
          cur_rc = shr(cur_rc, 1);
          L_temp = L_mult(cur_rc, A3);  // A3, 13.40   in Q11
          L_temp = L_sub(L_temp, L_B3); // L_B3, 11.16   in Q22
          L_temp = L_shr(L_temp, 11);
          LarNew[i] = extract_l(L_temp);
        }
      }
    }
    if (r_c[i] < 0)     // 和反射系数同号
    {
        LarNew[i] = sub(0, LarNew[i]);
    }
  }
```

### 对数面积比内插

子帧1为帧和前一帧的平均值，子帧2直接使用当前帧值（即LarNew）。

```CPP
// perc_var()
  /* Interpolation of Lar for the 1st subframe */
  temp = add(LarNew[0], LarOld[0]);
  Lar[0] = shr(temp, 1);
  LarOld[0] = LarNew[0];

  temp = add(LarNew[1], LarOld[1]);
  Lar[1] = shr(temp, 1);
  LarOld[1] = LarNew[1];
```

### 滤波器平滑

```cpp
// perc_var()
  for (k=0; k<2; k++)  /* LOOP : gamma2 for 1st to 2nd subframes */
  {
      /* ---------------------------------------------------------- */
      /* First criterion based on the first two Lars                */
      /* smooth == 1  ==>  gamma2 can vary from 0.4 to 0.7          */
      /* smooth == 0  ==>  gamma2 is set to 0.6                     */
      /*                                                            */
      /* Double threshold + hysteresis :                            */
      /* if smooth = 1                                              */
      /*  if (CritLar0 < THRESH_L1) and (CritLar1 > THRESH_H1)      */
      /*                                                 smooth = 0 */
      /* if smooth = 0                                              */
      /*  if (CritLar0 > THRESH_L2) or (CritLar1 < THRESH_H2)       */
      /*                                                 smooth = 1 */
      /* ---------------------------------------------------------- */

      CritLar0 = Lar[2*k];
      CritLar1 = Lar[2*k+1];

      if (smooth != 0) // smooth为全局变量
      {
        if ((sub(CritLar0,THRESH_L1)<0)&&( sub(CritLar1,THRESH_H1)>0)) 
        {   // THRESH_L1, -1.74 in Q11; THRESH_H1, -1.52 in Q11;
            smooth = 0;
        }
      }
      else  // smooth == 0
      {
        if ( (sub(CritLar0 ,THRESH_L2)>0) || (sub(CritLar1,THRESH_H2) <0) ) 
        {   // THRESH_L2, 0.65 in Q11; THRESH_H2, 0.43 in Q11;
            smooth = 1;
        }
      }

    if (smooth == 0) 
    {
      /* ------------------------------------------------------ */
      /* Second criterion based on the minimum distance between */
      /*                two successives LSPs                    */
      /*                                                        */
      /*           gamma2[k] = -6.0 * pi * d_min + 1.0          */
      /*                                                        */
      /*       with Lsfs normalized range 0.0 <= val <= 1.0     */
      /* ------------------------------------------------------ */
      gamma1[k] = GAMMA1_0;
      if (k == 0) 
      {
        Lsf = LsfInt;
      }
      else 
      {
        Lsf = LsfNew;
      }
      d_min = sub(Lsf[1], Lsf[0]);
      for (i=1; i<M-1; i++) 
      {
        temp = sub(Lsf[i+1],Lsf[i]); // Lsf[i+1]-Lsf[i]
        if (sub(temp,d_min)<0)       // min {Lsf[i+1]-Lsf[i]}
        {
            d_min = temp;
        }
      }
      temp = mult(ALPHA, d_min);    // ALPHA, 6*pi in Q10
      temp = sub(BETA, temp);       // BETA, 1 in Q10
      temp = shl(temp, 5);
      gamma2[k] = temp;

      // GAMMA2_0_L, 0.40 in Q15; GAMMA2_0_H, 0.70 in Q15
      if (sub(gamma2[k] , GAMMA2_0_H)>0) 
      {
        gamma2[k] = GAMMA2_0_H;
      }
      if (sub(gamma2[k] ,GAMMA2_0_L)<0) 
      {
        gamma2[k] = GAMMA2_0_L;
      }

    }
    else    // smooth == 1
    {       // GAMMA1_1, 0.94 in Q15; GAMMA2_1, 0.60 in Q15
      gamma1[k] = GAMMA1_1;
      gamma2[k] = GAMMA2_1;
    }
  }
```

## 语音信号感知加权

729感知加权的问题，同723一样，仍然对高频分量做一个估值，高频低，则认为频谱是倾斜的，要感知加权的时候应做适当的补偿。 因为人类的语音信息，在前两个共振峰能量较大，高频共振峰幅度则会依次回落。   
对比723，729引入了感知加权自适应机制。即：  
两个共振峰的位置如果太接近了，带宽扩展有可能导致两个共振峰出现重合了，也就是最近的两个共振峰越接近，则带宽扩展的加权越低\(即a越大\) 。

```math
W(z) = A *\frac{\dfrac{z}{gamma1}}{\dfrac{z}{gamma2}} 

sw(n)=s(n)+\sum_{i=1}^{10}{a_i \gamma_{1}^{i}\,s(n-i)}-
           \sum_{i=1}^{10}{a_i \gamma_{2}^{i}\, sw(n-i)},\qquad
           (n=0,\ldots,39)
```

两个子帧分别处理。每个子帧计算一次。

```cpp
// 子帧1, 0 - 39
Weight_Az(&A_t[0], gamma1[0], M, Ap1);
Weight_Az(&A_t[0], gamma2[0], M, Ap2);
Residu(Ap1, &speech[0], &wsp[0], L_SUBFR);
Syn_filt(Ap2, &wsp[0], &wsp[0], L_SUBFR, mem_w, 1);
// 子帧2, 40 - 79
Weight_Az(&A_t[MP1], gamma1[1], M, Ap1);
Weight_Az(&A_t[MP1], gamma2[1], M, Ap2);
Residu(Ap1, &speech[L_SUBFR], &wsp[L_SUBFR], L_SUBFR);
Syn_filt(Ap2, &wsp[L_SUBFR], &wsp[L_SUBFR], L_SUBFR, mem_w, 1);
```

### 计算加权滤波器系数

```math
ap[i]  =  a[i] * (gamma ^ i)
```

计算代码：

```cpp
// Weight_Az()
ap[0] = a[0];
fac   = gamma;
for(i=1; i<m; i++)
{
    ap[i] = round( L_mult(a[i], fac) );
    fac   = round( L_mult(fac, gamma) );
}
ap[m] = round( L_mult(a[m], fac) );
```

## 语音信号感知加权

滑动平均项：

```cpp
// Residu()
  for (i = 0; i < lg; i++)  // 误差滤波器，卷积滤波
  {
    s = L_mult(x[i], a[0]);
    for (j = 1; j <= M; j++) 
    {
      s = L_mac(s, a[j], x[i-j]);   // a[j]*x[i-j]+s
    }

    s = L_shl(s, 3);
    y[i] = round(s);
  }
```

自回归项：

```cpp
// Syn_filt()
  for(i=0; i<M; i++)
  {
    *yy++ = mem[i];     // filter initial state
  }

  /* Do the filtering. */   // 1/A(z), 自回归滤波
  for (i = 0; i < lg; i++)
  {
    s = L_mult(x[i], a[0]);
    for (j = 1; j <= M; j++)
    {
      s = L_msu(s, a[j], yy[-j]);   // a[j]*yy[-j]-s
    }

    s = L_shl(s, 3);
    *yy++ = round(s);
  }

  for(i=0; i<lg; i++)
  {
    y[i] = tmp[i+M];    // 
  }

  /* Update of memory if update==1 */
  if(update != 0)
  {
     for (i = 0; i < M; i++)
    {
       mem[i] = y[lg-M+i];
    }
  }
```



