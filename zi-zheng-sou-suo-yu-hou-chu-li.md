# 增益量化

寻找xn、y1、y2之间的相关性，并记录相关结果的量化Q值。
```cpp
   /*-----------------------------------------------------*
    * - Quantization of gains.                            *
    *-----------------------------------------------------*/

    g_coeff_cs[0]     = g_coeff[0];                   /* <y1,y1> */
    exp_g_coeff_cs[0] = negate(g_coeff[1]);           /* Q-Format:XXX -> JPN  */
    g_coeff_cs[1]     = negate(g_coeff[2]);           /* (xn,y1) -> -2<xn,y1> */       
    exp_g_coeff_cs[1] = negate(add(g_coeff[3], 1));   /* Q-Format:XXX -> JPN  */

    // xn：基音搜索目标向量
    // y1: 滤波后自适应码书激励
    // y2：滤波后固定码书激励
    Corr_xy2( xn, y1, y2, g_coeff_cs, exp_g_coeff_cs );  /* Q0 Q0 Q12 ^Qx ^Q0 */
                         /* g_coeff_cs[3]:exp_g_coeff_cs[3] = <y2,y2>   */
                         /* g_coeff_cs[4]:exp_g_coeff_cs[4] = -2<xn,y2> */
                         /* g_coeff_cs[5]:exp_g_coeff_cs[5] = 2<y1,y2>  */

    *ana++ = Qua_gain(code, g_coeff_cs, exp_g_coeff_cs,
                      L_SUBFR, &gain_pit, &gain_code, temp);
```

搜索下式均方误差的最小值
```math
|(x-g_p y-g_c z)|^2=(x-g_p y-g_c z)^T*(x-g_p y-g_c z)
```

## 增益预测   
计算固定码书的增益，预测能量和增益

因为已经两输入除以2，这里`$\bar{E}$`从36修改为30。
27为`$E_{sum}$`的`$E_{fix}$`形式表示的定点，即`$E_{sum}=\tfrac{E_{fix}}{2^{27}}$`。
```math
T=\bar{E}-E
=\bar{E}-10\lg{\tfrac{E_{sum}}{40}}

= 30 - \left[10\lg{{E_{fix}}} - 10\lg{2^{27}} - 10\lg{40} \right]

= 30 - 10\lg{{E_{fix}}} + 10\lg{40} + 10\lg{2^{27}}

= 30 - 10\lg{{E_{fix}}} + 16.02 + 81.278

= 127.298 - 10\lg{{E_{fix}}}

= 127.298 - 3.0103*10\log_2{{E_{fix}}}

log_{10}{x}=\lg{x}=\lg2*\log_{2}{x}=0.30103*\log_{2}{x}


```

```cpp
/* MA gain prediction coeff ={0.68, 0.58, 0.34, 0.19} in Q13 */
Word16 pred[4] = { 5571, 4751, 2785, 1556 };
```


```math
P=\bar{E}-E+E^{m}=T+E^{m}=T+\sum_{i=0}^{3}b_iU^{m-i}

g_c=10^{\tfrac{\bar{E}-E+E^{m}}{20}}
=2^{3.3219*\tfrac{\bar{E}-E+E^{m}}{20}}
=2^{0.166*P}
```

```cpp
// gainpred.c->Gain_predict

  /*-------------------------------*
   * Energy coming from code       *
   *-------------------------------*/
   L_tmp = 0;
   for(i=0; i<L_subfr; i++) // 固定码书能量
     L_tmp = L_mac(L_tmp, code[i], code[i]);

  /*-----------------------------------------------------------------*
   *  Compute: means_ener - 10log10(ener_code/ L_sufr)               *
   *  Note: mean_ener change from 36 dB to 30 dB because input/2     *
   *                                                                 *
   * = 30.0 - 10 log10( ener_code / lcode)  + 10log10(2^27)          *
   *                                          !!ener_code in Q27!!   *
   * = 30.0 - 3.0103 * log2(ener_code) + 10log10(40) + 10log10(2^27) *
   * = 30.0 - 3.0103 * log2(ener_code) + 16.02  + 81.278             *
   * = 127.298 - 3.0103 * log2(ener_code)                            *
   *-----------------------------------------------------------------*/
   // exp包含求log结果的整数部分，frac包含小数部分。
   Log2(L_tmp, &exp, &frac);               /* Q27->Q0 ^Q0 ^Q15       */
   L_tmp = Mpy_32_16(exp, frac, -24660);   /* Q0 Q15 Q13 -> ^Q14     */
                                           /* hi:Q0+Q13+1            */
                                           /* lo:Q15+Q13-15+1        */
                                           /* -24660[Q13]=-3.0103    */
   L_tmp = L_mac(L_tmp, 32588, 32);        /* 32588*32[Q14]=127.298  */

  /*-----------------------------------------------------------------*
   * Compute gcode0.                                                 *
   *  = Sum(i=0,3) pred[i]*past_qua_en[i] - ener_code + mean_ener    *
   *-----------------------------------------------------------------*/

   L_tmp = L_shl(L_tmp, 10);                      /* From Q14 to Q24 */
   for(i=0; i<4; i++)
     L_tmp = L_mac(L_tmp, pred[i], past_qua_en[i]); /* Q13*Q10 ->Q24 */

   *gcode0 = extract_h(L_tmp);                    /* From Q24 to Q8  */

  /*-----------------------------------------------------------------*
   * gcode0 = pow(10.0, gcode0/20)                                   *
   *        = pow(2, 3.3219*gcode0/20)                               *
   *        = pow(2, 0.166*gcode0)                                   *
   *-----------------------------------------------------------------*/
   L_tmp = L_mult(*gcode0, 5439);       /* *0.166 in Q15, result in Q24*/
   L_tmp = L_shr(L_tmp, 8);             /* From Q24 to Q16             */
   L_Extract(L_tmp, &exp, &frac);       /* Extract exponent of gcode0  */

   *gcode0 = extract_l(Pow2(14, frac)); /* Put 14 as exponent so that  */
                                        /* output of Pow2() will be:   */
                                        /* 16768 < Pow2() <= 32767     */
   *exp_gcode0 = sub(14,exp);
```





```cpp
   /*------------------------------------------------------------*
    * - Update pitch sharpening "sharp" with quantized gain_pit  *
    *------------------------------------------------------------*/

    sharp = gain_pit;
    if (sub(sharp, SHARPMAX) > 0) { sharp = SHARPMAX;         }
    if (sub(sharp, SHARPMIN) < 0) { sharp = SHARPMIN;         }

   /*------------------------------------------------------*
    * - Find the total excitation                          *
    * - find synthesis speech corresponding to exc[]       *
    * - update filters memories for finding the target     *
    *   vector in the next subframe                        *
    *   (update error[-m..-1] and mem_w_err[])             *
    *   update error function for taming process           *
    *------------------------------------------------------*/

    for (i = 0; i < L_SUBFR;  i++)
    {
      /* exc[i] = gain_pit*exc[i] + gain_code*code[i]; */
      /* exc[i]  in Q0   gain_pit in Q14               */
      /* code[i] in Q13  gain_cod in Q1                */

      L_temp = L_mult(exc[i+i_subfr], gain_pit);
      L_temp = L_mac(L_temp, code[i], gain_code);
      L_temp = L_shl(L_temp, 1);
      exc[i+i_subfr] = round(L_temp);
    }

    update_exc_err(gain_pit, T0);

    Syn_filt(Aq, &exc[i_subfr], &synth[i_subfr], L_SUBFR, mem_syn, 1);

    for (i = L_SUBFR-M, j = 0; i < L_SUBFR; i++, j++)
    {
      mem_err[j] = sub(speech[i_subfr+i], synth[i_subfr+i]);
      temp       = extract_h(L_shl( L_mult(y1[i], gain_pit),  1) );
      k          = extract_h(L_shl( L_mult(y2[i], gain_code), 2) );
      mem_w0[j]  = sub(xn[i], add(temp, k));
    }

    A  += MP1;           /* interpolated LPC parameters for next subframe */
    Aq += MP1;

```

## 预筛选

一个定点数保存为了两部分，二进制表示和小数位的bit数，即16Qx中的x。

```math


```

计算自适应码书基音增益，未量化。
```cpp
  /*-----------------------------------------------------------------*
   *  pre-selection                                                  *
   *-----------------------------------------------------------------*/
  /*-----------------------------------------------------------------*
   *  calculate best gain                                            *
   *                                                                 *
   *  tmp = -1./(4.*coeff[0]*coeff[2]-coeff[4]*coeff[4]) ;           *
   *  best_gain[0] = (2.*coeff[2]*coeff[1]-coeff[3]*coeff[4])*tmp ;  *
   *  best_gain[1] = (2.*coeff[0]*coeff[3]-coeff[1]*coeff[4])*tmp ;  *
   *  gbk_presel(best_gain,&cand1,&cand2,gcode0) ;                   *
   *                                                                 *
   *-----------------------------------------------------------------*/

  /*-----------------------------------------------------------------*
   *  tmp = -1./(4.*coeff[0]*coeff[2]-coeff[4]*coeff[4]) ;           *
   *-----------------------------------------------------------------*/
   L_tmp1 = L_mult( g_coeff[0], g_coeff[2] );
   exp1   = add( add( exp_coeff[0], exp_coeff[2] ), 1-2 ); // 相当于乘以 4
   L_tmp2 = L_mult( g_coeff[4], g_coeff[4] );
   exp2   = add( add( exp_coeff[4], exp_coeff[4] ), 1 );

   if( sub(exp1, exp2)>0 ){     // Q(L_tmp1) < Q(L_tmp2)，将L_tmp1化成和L_tmp2同定点后做减法
      L_tmp = L_sub( L_shr( L_tmp1, sub(exp1,exp2) ), L_tmp2 );     // L_tmp1>>(exp1-exp2) - L_tmp2
      exp = exp2; 
   }
   else{
      L_tmp = L_sub( L_tmp1, L_shr( L_tmp2, sub(exp2,exp1) ) );     // Q(L_tmp1) > Q(L_tmp2)，将L_tmp2化成和L_tmp1同定点后做减法
      exp = exp1;
   }
   sft = norm_l( L_tmp ); 
   denom = extract_h( L_shl(L_tmp, sft) );      // 归一化 
   exp_denom = sub( add( exp, sft ), 16 );      // 重新计算Q值

   inv_denom = div_s(16384,denom);              
   inv_denom = negate( inv_denom );             // -1 / denom
   exp_inv_denom = sub( 14+15, exp_denom );     // 计算Q值

  /*-----------------------------------------------------------------*
   *  best_gain[0] = (2.*coeff[2]*coeff[1]-coeff[3]*coeff[4])*tmp ;  *
   *-----------------------------------------------------------------*/
   L_tmp1 = L_mult( g_coeff[2], g_coeff[1] );
   exp1   = add( exp_coeff[2], exp_coeff[1] );
   L_tmp2 = L_mult( g_coeff[3], g_coeff[4] );
   exp2   = add( add( exp_coeff[3], exp_coeff[4] ), 1 );

   if( sub(exp1, exp2)>0 ){
      L_tmp = L_sub( L_shr( L_tmp1, add(sub(exp1,exp2),1 )), L_shr( L_tmp2,1 ) );
      exp = sub(exp2,1);
   }
   else{
      L_tmp = L_sub( L_shr( L_tmp1,1 ), L_shr( L_tmp2, add(sub(exp2,exp1),1 )) );
      exp = sub(exp1,1);
   }
   sft = norm_l( L_tmp );
   nume = extract_h( L_shl(L_tmp, sft) );
   exp_nume = sub( add( exp, sft ), 16 );

   sft = sub( add( exp_nume, exp_inv_denom ), (9+16-1) );
   L_acc = L_shr( L_mult( nume,inv_denom ), sft );      // 乘以 tmp
   best_gain[0] = extract_h( L_acc );             /*-- best_gain[0]:Q9 --*/

   if (tameflag == 1){                  // 需要柔和，即自适应搜索基音延迟误差较大
        if(sub(best_gain[0], GPCLIP2) > 0) 
            best_gain[0] = GPCLIP2;     // 如果需要柔和，设置基音增益最大值
   }

  /*-----------------------------------------------------------------*
   *  best_gain[1] = (2.*coeff[0]*coeff[3]-coeff[1]*coeff[4])*tmp ;  *
   *-----------------------------------------------------------------*/
   L_tmp1 = L_mult( g_coeff[0], g_coeff[3] );
   exp1   = add( exp_coeff[0], exp_coeff[3] ) ;
   L_tmp2 = L_mult( g_coeff[1], g_coeff[4] );
   exp2   = add( add( exp_coeff[1], exp_coeff[4] ), 1 );

   if( sub(exp1, exp2)>0 ){
      L_tmp = L_sub( L_shr( L_tmp1, add(sub(exp1,exp2),1) ), L_shr( L_tmp2,1 ) );
      exp = sub(exp2,1);
   }
   else{
      L_tmp = L_sub( L_shr( L_tmp1,1 ), L_shr( L_tmp2, add(sub(exp2,exp1),1) ) );
      exp = sub(exp1,1);
   }
   sft = norm_l( L_tmp );
   nume = extract_h( L_shl(L_tmp, sft) );
   exp_nume = sub( add( exp, sft ), 16 );

   sft = sub( add( exp_nume, exp_inv_denom ), (2+16-1) );
   L_acc = L_shr( L_mult( nume,inv_denom ), sft );      // // 乘以 tmp
   best_gain[1] = extract_h( L_acc );             /*-- best_gain[1]:Q2 --*/

  /*--- Change Q-format of gcode0 ( Q[exp_gcode0] -> Q4 ) ---*/
   if( sub(exp_gcode0,4) >= 0 ){
      gcode0_org = shr( gcode0, sub(exp_gcode0,4) );
   }
   else{
      L_acc = L_deposit_l( gcode0 );        // 16->32位扩展，扩展符号位
      L_acc = L_shl( L_acc, sub( (4+16), exp_gcode0 ) );    // 20 - exp_gcode0
      gcode0_org = extract_h( L_acc );              /*-- gcode0_org:Q4 --*/
   }
```

预筛选共轭码书，得到搜索起始位置。


```cpp
// Gbk_presel，返回后续搜索的起始位置
 /*--------------------------------------------------------------------------*
   x = (best_gain[1] - (coef[0][0]*best_gain[0] + coef[1][1])*gcode0) * inv_coef;
  *--------------------------------------------------------------------------*/
   L_cfbg = L_mult( coef[0][0], best_gain[0] );        /* L_cfbg:Q20 -> !!y */
   L_acc = L_shr( L_coef[1][1], 15 );                  /* L_acc:Q20     */
   L_acc = L_add( L_cfbg , L_acc );
   acc_h = extract_h( L_acc );                         /* acc_h:Q4      */
   L_preg = L_mult( acc_h, gcode0 );                   /* L_preg:Q9     */
   L_acc = L_shl( L_deposit_l( best_gain[1] ), 7 );    /* L_acc:Q9      */
   L_acc = L_sub( L_acc, L_preg );
   acc_h = extract_h( L_shl( L_acc,2 ) );              /* L_acc_h:Q[-5] */
   L_tmp_x = L_mult( acc_h, INV_COEF );                /* L_tmp_x:Q15   */

 /*--------------------------------------------------------------------------*
   y = (coef[1][0]*(-coef[0][1]+best_gain[0]*coef[0][0])*gcode0
                                      -coef[0][0]*best_gain[1]) * inv_coef;
  *--------------------------------------------------------------------------*/
   L_acc = L_shr( L_coef[0][1], 10 );                  /* L_acc:Q20   */
   L_acc = L_sub( L_cfbg, L_acc );                     /* !!x -> L_cfbg:Q20 */
   acc_h = extract_h( L_acc );                         /* acc_h:Q4    */
   acc_h = mult( acc_h, gcode0 );                      /* acc_h:Q[-7] */
   L_tmp = L_mult( acc_h, coef[1][0] );                /* L_tmp:Q10   */

   L_preg = L_mult( coef[0][0], best_gain[1] );        /* L_preg:Q13  */
   L_acc = L_sub( L_tmp, L_shr(L_preg,3) );            /* L_acc:Q10   */

   acc_h = extract_h( L_shl( L_acc,2 ) );              /* acc_h:Q[-4] */
   L_tmp_y = L_mult( acc_h, INV_COEF );                /* L_tmp_y:Q16 */
   
   // 搜索
   sft_y = (14+4+1)-16;         /* (Q[thr1]+Q[gcode0]+1)-Q[L_tmp_y] */
   sft_x = (15+4+1)-15;         /* (Q[thr2]+Q[gcode0]+1)-Q[L_tmp_x] */

    if(gcode0 > 0){
      /*-- pre select codebook #1 --*/
      *cand1 = 0 ;
        do{ // thr1，门限，递增，共(NCODE1-NCAN1=4)个
            L_temp = L_sub( L_tmp_y, L_shr(L_mult(thr1[*cand1],gcode0),sft_y));     // L_tmp_y - [(thr1[*cand1]*gcode0) >> sft_y]
            if(L_temp >0L ) // 大于0，则继续 
                (*cand1) =add(*cand1,1);
            else                        
                break ;      // 从这里退出表示符合要求的个数大于要求个数，但也只取前面的 
        } while(sub((*cand1),(NCODE1-NCAN1))<0) ;   // NCODE1码字总个数，NCAN1选择的码字个数，保证候选码字个数足够
        
        /*-- pre select codebook #2 --*/
        *cand2 = 0 ;
        do{ // thr2，门限，递增，共(NCODE2-NCAN2=8)个
            L_temp = L_sub( L_tmp_x , L_shr(L_mult(thr2[*cand2],gcode0),sft_x));    // L_tmp_x - [(thr2[*cand2]*gcode0) >> sft_x]
            if( L_temp >0L) // 大于0，则继续 
                (*cand2) =add(*cand2,1);
            else               
                break ;     // 从这里退出表示符合要求的个数大于要求个数，但也只取前面的 
        } while(sub((*cand2),(NCODE2-NCAN2))<0);    // NCODE2码字总个数，NCAN2选择的码字个数，保证候选码字个数足够
    }else{
        /*-- pre select codebook #1 --*/
        *cand1 = 0 ;
        do{ // thr1，门限，递增，共(NCODE1-NCAN1=4)个
            L_temp = L_sub(L_tmp_y ,L_shr(L_mult(thr1[*cand1],gcode0),sft_y));      // L_tmp_y - [(thr1[*cand1]*gcode0) >> sft_y]
            if( L_temp <0L) // 小于0，则继续搜索
                (*cand1) =add(*cand1,1);
            else               
                break ;     // 从这里退出表示符合要求的个数大于要求个数，但也只取前面的   
        } while(sub((*cand1),(NCODE1-NCAN1))) ;
        
        /*-- pre select codebook #2 --*/
        *cand2 = 0 ;
        do{ // thr2，门限，递增，共(NCODE2-NCAN2=8)个
            L_temp =L_sub(L_tmp_x ,L_shr(L_mult(thr2[*cand2],gcode0),sft_x));       // L_tmp_x - [(thr2[*cand2]*gcode0) >> sft_x]
            if( L_temp <0L) // 小于0，则继续搜索
                (*cand2) =add(*cand2,1);
            else               
                break ;     // 从这里退出表示符合要求的个数大于要求个数，但也只取前面的   
        } while(sub( (*cand2),(NCODE2-NCAN2))) ;
    }
```

### 执行搜索过程
```cpp
/*---------------------------------------------------------------------------*
 *                                                                           *
 * Find the best quantizer.                                                  *
 *                                                                           *
 *  dist_min = MAX_32;                                                       *
 *  for ( i=0 ; i<NCAN1 ; i++ ){                                             *
 *    for ( j=0 ; j<NCAN2 ; j++ ){                                           *
 *      g_pitch = gbk1[cand1+i][0] + gbk2[cand2+j][0];                       *

 *      g_code = gcode0 * (gbk1[cand1+i][1] + gbk2[cand2+j][1]);             *
 *      dist = g_pitch*g_pitch * coeff[0]                                    *
 *           + g_pitch         * coeff[1]                                    *
 *           + g_code*g_code   * coeff[2]                                    *
 *           + g_code          * coeff[3]                                    *
 *           + g_pitch*g_code  * coeff[4] ;                                  *
 *                                                                           *
 *      if (dist < dist_min){                                                *
 *        dist_min = dist;                                                   *
 *        indice1 = cand1 + i ;                                              *
 *        indice2 = cand2 + j ;                                              *
 *      }                                                                    *
 *    }                                                                      *
 *  }                                                                        *
 *                                                                           *
 * g_pitch         = Q13                                                     *
 * g_pitch*g_pitch = Q11:(13+13+1-16)                                        *
 * g_code          = Q[exp_gcode0-3]:(exp_gcode0+(13-1)+1-16)                *
 * g_code*g_code   = Q[2*exp_gcode0-21]:(exp_gcode0-3+exp_gcode0-3+1-16)     *
 * g_pitch*g_code  = Q[exp_gcode0-5]:(13+exp_gcode0-3+1-16)                  *
 *                                                                           *
 * term 0: g_pitch*g_pitch*coeff[0] ;exp_min0 = 13             +exp_coeff[0] *
 * term 1: g_pitch        *coeff[1] ;exp_min1 = 14             +exp_coeff[1] *
 * term 2: g_code*g_code  *coeff[2] ;exp_min2 = 2*exp_gcode0-21+exp_coeff[2] *
 * term 3: g_code         *coeff[3] ;exp_min3 = exp_gcode0  - 3+exp_coeff[3] *
 * term 4: g_pitch*g_code *coeff[4] ;exp_min4 = exp_gcode0  - 4+exp_coeff[4] *
 *                                                                           *
 *---------------------------------------------------------------------------*/

   exp_min[0] = add( exp_coeff[0], 13 );
   exp_min[1] = add( exp_coeff[1], 14 );
   exp_min[2] = add( exp_coeff[2], sub( shl( exp_gcode0, 1 ), 21 ) );
   exp_min[3] = add( exp_coeff[3], sub( exp_gcode0, 3 ) );
   exp_min[4] = add( exp_coeff[4], sub( exp_gcode0, 4 ) );

   e_min = exp_min[0];
   for(i=1; i<5; i++){
      if( sub(exp_min[i], e_min) < 0 ){
         e_min = exp_min[i];
      }
   }

   /* align coeff[] and save in special 32 bit double precision */

   for(i=0; i<5; i++){
     j = sub( exp_min[i], e_min );
     L_tmp = L_deposit_h( g_coeff[i] );
     L_tmp = L_shr( L_tmp, j );          /* L_tmp:Q[exp_g_coeff[i]+16-j] */
     L_Extract( L_tmp, &coeff[i], &coeff_lsf[i] );          /* DPF */
   }

   /* Codebook search */

   L_dist_min = MAX_32;

   /* initialization used only to suppress Microsoft Visual C++  warnings */
   index1 = cand1;
   index2 = cand2;

if(tameflag == 1){
   for(i=0; i<NCAN1; i++){
      for(j=0; j<NCAN2; j++){
         g_pitch = add( gbk1[cand1+i][0], gbk2[cand2+j][0] );     /* Q14 */
         if(g_pitch < GP0999) {
         L_acc = L_deposit_l( gbk1[cand1+i][1] );
         L_accb = L_deposit_l( gbk2[cand2+j][1] );                /* Q13 */
         L_tmp = L_add( L_acc,L_accb );
         tmp = extract_l( L_shr( L_tmp,1 ) );                     /* Q12 */

         g_code   = mult( gcode0, tmp );         /*  Q[exp_gcode0+12-15] */
         g2_pitch = mult(g_pitch, g_pitch);                       /* Q13 */
         g2_code  = mult(g_code,  g_code);       /* Q[2*exp_gcode0-6-15] */
         g_pit_cod= mult(g_code,  g_pitch);      /* Q[exp_gcode0-3+14-15] */

         L_tmp = Mpy_32_16(coeff[0], coeff_lsf[0], g2_pitch);
         L_tmp = L_add(L_tmp, Mpy_32_16(coeff[1], coeff_lsf[1], g_pitch) );
         L_tmp = L_add(L_tmp, Mpy_32_16(coeff[2], coeff_lsf[2], g2_code) );
         L_tmp = L_add(L_tmp, Mpy_32_16(coeff[3], coeff_lsf[3], g_code) );
         L_tmp = L_add(L_tmp, Mpy_32_16(coeff[4], coeff_lsf[4], g_pit_cod) );

         L_temp = L_sub(L_tmp, L_dist_min);

         if( L_temp < 0L ){
            L_dist_min = L_tmp;
            index1 = add(cand1,i);
            index2 = add(cand2,j);
         }
        }
      }
   }

}
else{
   for(i=0; i<NCAN1; i++){
      for(j=0; j<NCAN2; j++){
         g_pitch = add( gbk1[cand1+i][0], gbk2[cand2+j][0] );     /* Q14 */
         L_acc = L_deposit_l( gbk1[cand1+i][1] );
         L_accb = L_deposit_l( gbk2[cand2+j][1] );                /* Q13 */
         L_tmp = L_add( L_acc,L_accb );
         tmp = extract_l( L_shr( L_tmp,1 ) );                     /* Q12 */

         g_code   = mult( gcode0, tmp );         /*  Q[exp_gcode0+12-15] */
         g2_pitch = mult(g_pitch, g_pitch);                       /* Q13 */
         g2_code  = mult(g_code,  g_code);       /* Q[2*exp_gcode0-6-15] */
         g_pit_cod= mult(g_code,  g_pitch);      /* Q[exp_gcode0-3+14-15] */

         L_tmp = Mpy_32_16(coeff[0], coeff_lsf[0], g2_pitch);
         L_tmp = L_add(L_tmp, Mpy_32_16(coeff[1], coeff_lsf[1], g_pitch) );
         L_tmp = L_add(L_tmp, Mpy_32_16(coeff[2], coeff_lsf[2], g2_code) );
         L_tmp = L_add(L_tmp, Mpy_32_16(coeff[3], coeff_lsf[3], g_code) );
         L_tmp = L_add(L_tmp, Mpy_32_16(coeff[4], coeff_lsf[4], g_pit_cod) );

         L_temp = L_sub(L_tmp, L_dist_min);

         if( L_temp < 0L ){
            L_dist_min = L_tmp;
            index1 = add(cand1,i);
            index2 = add(cand2,j);
         }

      }
   }
}

```

### 得到搜索结果


```cpp
// 
   /* Read the quantized gain */

  /*-----------------------------------------------------------------*
   * *gain_pit = gbk1[indice1][0] + gbk2[indice2][0];                *
   *-----------------------------------------------------------------*/
   *gain_pit = add( gbk1[index1][0], gbk2[index2][0] );      /* Q14 */

  /*-----------------------------------------------------------------*
   * *gain_code = (gbk1[indice1][1]+gbk2[indice2][1]) * gcode0;      *
   *-----------------------------------------------------------------*/
   L_acc = L_deposit_l( gbk1[index1][1] );
   L_accb = L_deposit_l( gbk2[index2][1] );
   L_gbk12 = L_add( L_acc, L_accb );                          /* Q13 */
   tmp = extract_l( L_shr( L_gbk12,1 ) );                     /* Q12 */
   L_acc = L_mult(tmp, gcode0);                /* Q[exp_gcode0+12+1] */

   L_acc = L_shl(L_acc, add( negate(exp_gcode0),(-12-1+1+16) ));
   *gain_cod = extract_h( L_acc );   

```

### 更新量化器参数

```math
U^{(m)}=E^{(m)}-\tilde{E}^{(m)}=20\log{(\gamma)}

\gamma=L\_gbk12=gbk1[index1][1]+gbk2[index2][1]
```

```cpp
// gainpred.c->Gain_update()

   for(i=3; i>0; i--){  // 更新量化能量误差记录
      past_qua_en[i] = past_qua_en[i-1];         /* Q10 */
   }
  /*----------------------------------------------------------------------*
   * -- past_qua_en[0] = 20*log10(gbk1[index1][1]+gbk2[index2][1]); --    *
   *    2 * 10 log10( gbk1[index1][1]+gbk2[index2][1] )                   *
   *  = 2 * 3.0103 log2( gbk1[index1][1]+gbk2[index2][1] )                *
   *  = 2 * 3.0103 log2( gbk1[index1][1]+gbk2[index2][1] )                *
   *                                                 24660:Q12(6.0205)    *
   *----------------------------------------------------------------------*/
   
   Log2( L_gbk12, &exp, &frac );               /* L_gbk12:Q13       */
   L_acc = L_Comp( sub(exp,13), frac);         /* L_acc:Q16           */ // DPF->WORD32
   tmp = extract_h( L_shl( L_acc,13 ) );       /* tmp:Q13           */
   past_qua_en[0] = mult( tmp, 24660 );        /* past_qua_en[]:Q10 */   
```

### 更新sharp

```cpp
   /*------------------------------------------------------------*
    * - Update pitch sharpening "sharp" with quantized gain_pit  *
    *------------------------------------------------------------*/

    sharp = gain_pit;
    if (sub(sharp, SHARPMAX) > 0) { sharp = SHARPMAX;         }
    if (sub(sharp, SHARPMIN) < 0) { sharp = SHARPMIN;         }
```

# 子帧搜索后处理
为下一子帧处理做准备。
```cpp
//--
   /*------------------------------------------------------*
    * - Find the total excitation                          *
    * - find synthesis speech corresponding to exc[]       *
    * - update filters memories for finding the target     *
    *   vector in the next subframe                        *
    *   (update error[-m..-1] and mem_w_err[])             *
    *   update error function for taming process           *
    *------------------------------------------------------*/

    for (i = 0; i < L_SUBFR;  i++)
    {
      /* exc[i] = gain_pit*exc[i] + gain_code*code[i]; */
      /* exc[i]  in Q0   gain_pit in Q14               */
      /* code[i] in Q13  gain_cod in Q1                */

      L_temp = L_mult(exc[i+i_subfr], gain_pit);
      L_temp = L_mac(L_temp, code[i], gain_code);
      L_temp = L_shl(L_temp, 1);
      exc[i+i_subfr] = round(L_temp);
    }

    update_exc_err(gain_pit, T0);

    Syn_filt(Aq, &exc[i_subfr], &synth[i_subfr], L_SUBFR, mem_syn, 1);

    for (i = L_SUBFR-M, j = 0; i < L_SUBFR; i++, j++)
    {
      mem_err[j] = sub(speech[i_subfr+i], synth[i_subfr+i]);
      temp       = extract_h(L_shl( L_mult(y1[i], gain_pit),  1) );
      k          = extract_h(L_shl( L_mult(y2[i], gain_code), 2) );
      mem_w0[j]  = sub(xn[i], add(temp, k));
    }

    // MP1=M+1，M为LP滤波器阶数
    A  += MP1;           /* interpolated LPC parameters for next subframe */
    Aq += MP1;

```

## 计算总的激励

## 计算自适应码字误差
计算编码端和解码段自适应码字的误差。

```cpp
    L_worst = -1L;
    n = sub(T0, L_SUBFR);

    if(n < 0) {     // 基音延迟小于子帧长度
        L_Extract(L_exc_err[0], &hi, &lo);
        L_temp = Mpy_32_16(hi, lo, gain_pit);   // L_exc_err[0]*gain_pit（基音增益）
        L_temp = L_shl(L_temp, 1);
        L_temp = L_add(0x00004000L, L_temp);    // L_exc_err[0]*gain_pit*2 - L_worst
        L_acc = L_sub(L_temp, L_worst);
        if(L_acc > 0L) {
                L_worst = L_temp;
        }
        L_Extract(L_temp, &hi, &lo);
        L_temp = Mpy_32_16(hi, lo, gain_pit);
        L_temp = L_shl(L_temp, 1);
        L_temp = L_add(0x00004000L, L_temp);    // L_temp[0]*gain_pit*2 - L_worst
        L_acc = L_sub(L_temp, L_worst);
        if(L_acc > 0L) {
                L_worst = L_temp;
        }
    } else {  // 基音延迟大于子帧长度
        zone1 = tab_zone[n];    // Table for taming procedure test_err.

        i = sub(T0, 1);
        zone2 = tab_zone[i];

        for(i = zone1; i <= zone2; i++) {
            L_Extract(L_exc_err[i], &hi, &lo);
            L_temp = Mpy_32_16(hi, lo, gain_pit);
            L_temp = L_shl(L_temp, 1);
            L_temp = L_add(0x00004000L, L_temp);
            L_acc = L_sub(L_temp, L_worst);
            if(L_acc > 0L) {
                L_worst = L_temp;
            }
        }
    }

    for(i=3; i>=1; i--) {
        L_exc_err[i] = L_exc_err[i-1];
    }
    L_exc_err[0] = L_worst;

```

## 更新滤波器状态

```cpp
//
    // 得到语音合成结果  
    Syn_filt(Aq, &exc[i_subfr], &synth[i_subfr], L_SUBFR, mem_syn, 1);

    for (i = L_SUBFR-M, j = 0; i < L_SUBFR; i++, j++)   // M：LP滤波器借书
    {
      mem_err[j] = sub(speech[i_subfr+i], synth[i_subfr+i]);        // 合成误差
      temp       = extract_h(L_shl( L_mult(y1[i], gain_pit),  1) ); // (y1[i]*gain_pit)*2
      k          = extract_h(L_shl( L_mult(y2[i], gain_code), 2) ); // (y2[i]*gain_code)*4
      mem_w0[j]  = sub(xn[i], add(temp, k));    // xn[i] - (temp + k)
    }
```


# 存储更新

为下一帧做准备

```cpp
// Copy vector x[] to y[]
void Copy(
  Word16 x[],      /* (i)   : input vector   */
  Word16 y[],      /* (o)   : output vector  */
  Word16 L         /* (i)   : vector length  */
) 
```

```cpp
 /*--------------------------------------------------*
  * Update signal for next frame.                    *
  * -> shift to the left by L_FRAME:                 *
  *     speech[], wsp[] and  exc[]                   *
  *--------------------------------------------------*/
  Copy(&old_speech[L_FRAME], &old_speech[0], L_TOTAL-L_FRAME);
  Copy(&old_wsp[L_FRAME], &old_wsp[0], PIT_MAX);
  Copy(&old_exc[L_FRAME], &old_exc[0], PIT_MAX+L_INTERPOL);
```

# 最后

# 参考文献
1. [g729源码分析-2-共振锋感知加权]
2. [深蓝怒火的CSDN博客]
3. [JBoss - org.mobicents.media.server.impl.dsp.audio.g729.LD8KConstants]


[g729源码分析-2-共振锋感知加权]:http://blog.csdn.net/lsccsl/article/details/7449361
[深蓝怒火的CSDN博客]:http://blog.csdn.net/lsccsl/article/list/3
[JBoss - org.mobicents.media.server.impl.dsp.audio.g729.LD8KConstants]:http://grepcode.com/file/repository.jboss.org/nexus/content/repositories/releases/org.mobicents.servers.media/mms-impl/2.0.0.GA/org/mobicents/media/server/impl/dsp/audio/g729/LD8KConstants.java#LD8KConstants.0PI04